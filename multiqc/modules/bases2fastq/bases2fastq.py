from collections import defaultdict
import copy
from itertools import chain
import re
import json
import logging
import random
from typing import Any, Dict, List, Optional
import uuid
from pathlib import Path

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.types import LoadedFileDict
from multiqc.utils import mqc_colour

from multiqc.modules.bases2fastq.plot_runs import (
    plot_run_stats,
    tabulate_manifest_stats,
    tabulate_index_assignment_stats,
    tabulate_unassigned_index_stats,
    tabulate_run_stats,
    tabulate_project_stats,
    plot_base_quality_hist,
    plot_base_quality_by_cycle,
)
from multiqc.modules.bases2fastq.plot_samples import (
    tabulate_sample_stats,
    sequence_content_plot,
    plot_per_cycle_N_content,
    plot_adapter_content,
    plot_per_read_gc_hist,
)

log = logging.getLogger(__name__)


# Default minimum polony threshold - samples below this are skipped
DEFAULT_MIN_POLONIES = 10000


def _get_min_polonies() -> int:
    """
    Get the minimum polonies threshold from config or use default.

    Can be configured in multiqc_config.yaml:
        bases2fastq_config:
            min_polonies: 5000
    """
    cfg = getattr(config, "bases2fastq_config", {})
    if not isinstance(cfg, dict):
        return DEFAULT_MIN_POLONIES

    min_polonies = cfg.get("min_polonies", DEFAULT_MIN_POLONIES)
    try:
        min_polonies = int(min_polonies)
    except (ValueError, TypeError):
        log.warning(f"Invalid min_polonies value '{min_polonies}', using default {DEFAULT_MIN_POLONIES}")
        min_polonies = DEFAULT_MIN_POLONIES

    if min_polonies != DEFAULT_MIN_POLONIES:
        log.debug(f"Using custom min_polonies threshold: {min_polonies}")

    return min_polonies


class MultiqcModule(BaseMultiqcModule):
    """
    Bases2Fastq is Element Biosciences' secondary analysis software for demultiplexing
    sequencing data from AVITI systems and converting base calls into FASTQ files.

    The module parses the following output files from Bases2Fastq:

    - `RunStats.json`: Contains run-level and sample-level QC metrics
    - `RunManifest.json`: Contains sample sheet information including indexing and adapter settings
    - Project-level `RunStats.json`: Contains project-specific metrics when demultiplexing by project

    The module supports both run-level analysis (single run) and project-level analysis
    (aggregated metrics across projects), displaying metrics such as:

    - Polony counts and yields
    - Base quality distributions
    - Index assignment statistics
    - Per-sample sequence content and GC distribution
    - Adapter content analysis
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Bases2Fastq",
            anchor="bases2fastq",
            href="https://docs.elembio.io/docs/bases2fastq/introduction/",
            info="Demultiplexes and converts Element AVITI base calls into FASTQ files",
            doi="10.1038/s41587-023-01750-7",
        )

        # Get configurable minimum polonies threshold
        self.min_polonies = _get_min_polonies()

        # Initialize data structures
        self._init_data_structures()

        # Parse and validate input data
        summary_path = self._parse_and_validate_data()

        # Select data based on summary path and parse additional sources
        run_data, sample_data, samples_to_projects, manifest_data, index_assignment_data, unassigned_sequences = (
            self._select_data_by_summary_path(summary_path)
        )

        # Set up color schemes for groups and samples
        self._setup_colors(sample_data, samples_to_projects, summary_path)

        # Generate all plots and sections
        self._generate_plots(
            summary_path, run_data, sample_data, samples_to_projects,
            manifest_data, index_assignment_data, unassigned_sequences
        )

        # Write main data file at the very end after all sections are added
        self.write_data_file(sample_data, "bases2fastq")

    def _init_data_structures(self) -> None:
        """Initialize all data structures used by the module."""
        # File cache to avoid reading the same JSON files multiple times
        self._file_cache: Dict[str, Any] = {}

        # Run, project and sample level structures
        self.run_level_data: Dict[str, Any] = {}
        self.run_level_samples: Dict[str, Any] = {}
        self.run_level_samples_to_project: Dict[str, str] = {}
        self.project_level_data: Dict[str, Any] = {}
        self.project_level_samples: Dict[str, Any] = {}
        self.project_level_samples_to_project: Dict[str, str] = {}

        # Run and project groups
        self.group_dict: Dict[str, Any] = {}
        self.group_lookup_dict: Dict[str, Any] = {}
        self.project_lookup_dict: Dict[str, Any] = {}

        # Additional data structures
        self.b2f_sample_data: Dict[str, Any] = {}
        self.b2f_run_data: Dict[str, Any] = {}
        self.b2f_run_project_data: Dict[str, Any] = {}
        self.b2f_run_project_sample_data: Dict[str, Any] = {}
        self.missing_runs: set = set()
        self.sample_id_to_run: Dict[str, str] = {}

    def _read_json_file(self, file_path: Path) -> Optional[Dict[str, Any]]:
        """
        Read and parse a JSON file with caching.

        Args:
            file_path: Path to the JSON file

        Returns:
            Parsed JSON data or None if reading failed
        """
        cache_key = str(file_path.resolve())

        if cache_key in self._file_cache:
            return self._file_cache[cache_key]

        if not file_path.exists():
            log.error(
                f"{file_path.name} does not exist at {file_path}.\n"
                f"Please visit Elembio online documentation for more information - "
                f"https://docs.elembio.io/docs/bases2fastq/introduction/"
            )
            return None

        try:
            with open(file_path) as _infile:
                data = json.load(_infile)
                self._file_cache[cache_key] = data
                return data
        except (json.JSONDecodeError, OSError) as e:
            log.error(f"Error reading {file_path}: {e}")
            return None

    def _parse_and_validate_data(self) -> str:
        """
        Parse input data and validate that samples were found.

        Returns:
            summary_path: The determined summary path ('run_level', 'project_level', or 'combined_level')
        """
        # Check for available log files
        run_level_log_files = len(list(self.find_log_files("bases2fastq/run")))
        project_level_log_files = len(list(self.find_log_files("bases2fastq/project")))

        if run_level_log_files == 0 and project_level_log_files == 0:
            error_msg = "No run- or project-level log files found within the Bases2Fastq results."
            log.error(error_msg)
            raise ModuleNoSamplesFound(error_msg)

        # Parse data from available sources
        if run_level_log_files > 0:
            (self.run_level_data, self.run_level_samples, self.run_level_samples_to_project) = (
                self._parse_run_project_data("bases2fastq/run")
            )
        if project_level_log_files > 0:
            (self.project_level_data, self.project_level_samples, self.project_level_samples_to_project) = (
                self._parse_run_project_data("bases2fastq/project")
            )

        # Count samples
        num_run_level_samples = len(self.run_level_samples)
        num_project_level_samples = len(self.project_level_samples)

        # Ensure at least some data was found
        if all([
            len(self.run_level_data) == 0,
            num_run_level_samples == 0,
            len(self.project_level_data) == 0,
            num_project_level_samples == 0,
        ]):
            error_msg = "No run-, project- or sample-level data found"
            log.error(error_msg)
            raise ModuleNoSamplesFound(error_msg)

        # Determine summary path
        summary_path = self._determine_summary_path()

        # Log what was found
        log.info(f"Found {len(self.run_level_data)} run(s) within the Bases2Fastq results.")
        log.info(f"Found {len(self.project_level_data)} project(s) within the Bases2Fastq results.")
        if summary_path == "run_level":
            log.info(f"Found {num_run_level_samples} sample(s) within the Bases2Fastq results.")
        else:
            log.info(f"Found {num_project_level_samples} sample(s) within the Bases2Fastq results.")

        # Required call to confirm module is used
        self.add_software_version(None)

        # Warn if no data found
        if len(self.run_level_data) == 0 and len(self.project_level_data) == 0:
            log.warning("No run/project stats found!")
        if num_run_level_samples == 0 and num_project_level_samples == 0:
            log.warning("No sample stats found!")

        return summary_path

    def _determine_summary_path(self) -> str:
        """
        Determine which summary path to use based on available data.

        Returns:
            'run_level', 'project_level', or 'combined_level'
        """
        has_run_data = len(self.run_level_data) > 0
        has_project_data = len(self.project_level_data) > 0

        if has_run_data and not has_project_data:
            return "run_level"
        elif not has_run_data and has_project_data:
            return "project_level"
        elif has_run_data and has_project_data:
            return "combined_level"
        else:
            error_msg = "No run- or project-level data was retained. No report will be generated."
            log.error(error_msg)
            raise ModuleNoSamplesFound(error_msg)

    def _select_data_by_summary_path(self, summary_path: str):
        """
        Select the appropriate data sources based on the summary path.

        Returns:
            Tuple of (run_data, sample_data, samples_to_projects, manifest_data,
                index_assignment_data, unassigned_sequences)
        """
        if summary_path == "run_level":
            return (
                self.run_level_data,
                self.run_level_samples,
                self.run_level_samples_to_project,
                self._parse_run_manifest("bases2fastq/manifest"),
                self._parse_index_assignment("bases2fastq/manifest"),
                self._parse_run_unassigned_sequences("bases2fastq/run"),
            )
        elif summary_path == "project_level":
            return (
                self.project_level_data,
                self.project_level_samples,
                self.project_level_samples_to_project,
                self._parse_run_manifest_in_project("bases2fastq/project"),
                self._parse_index_assignment_in_project("bases2fastq/project"),
                {},  # No unassigned sequences for project level
            )
        elif summary_path == "combined_level":
            return (
                self.run_level_data,
                self.project_level_samples,
                self.project_level_samples_to_project,
                self._parse_run_manifest("bases2fastq/manifest"),
                self._parse_index_assignment("bases2fastq/manifest"),
                self._parse_run_unassigned_sequences("bases2fastq/run"),
            )
        else:
            error_msg = "No run- or project-level data was retained. No report will be generated."
            log.error(error_msg)
            raise ModuleNoSamplesFound(error_msg)

    def _setup_colors(self, sample_data: Dict, samples_to_projects: Dict, summary_path: str) -> None:
        """Set up color schemes for groups and samples."""
        # Create run and project groups
        run_groups: Dict[str, List] = defaultdict(list)
        project_groups: Dict[str, List] = defaultdict(list)
        in_project_sample_groups: Dict[str, List] = defaultdict(list)
        ind_sample_groups: Dict[str, List] = defaultdict(list)

        for sample in sample_data.keys():
            run_name, _ = sample.split("__")
            run_groups[run_name].append(sample)
            sample_project = samples_to_projects[sample]
            project_groups[sample_project].append(sample)
            ind_sample_groups[sample] = [sample]
            if summary_path == "project_level":
                in_project_sample_groups[sample].append(sample)

        merged_groups = {**run_groups, **project_groups, **in_project_sample_groups, **ind_sample_groups}

        # Build color palette
        self.color_getter = mqc_colour.mqc_colour_scale()
        self.palette = list(chain.from_iterable(
            self.color_getter.get_colours(hue)
            for hue in ["Set2", "Pastel1", "Accent", "Set1", "Set3", "Dark2", "Paired", "Pastel2"]
        ))

        # Add extra colors if needed
        if len(merged_groups) > len(self.palette):
            extra_colors = [
                f"#{random.randrange(0, 0xFFFFFF):06x}" for _ in range(len(self.palette), len(merged_groups))
            ]
            self.palette = self.palette + extra_colors

        # Assign colors to groups
        self.group_color = {
            group: color for group, color in zip(merged_groups.keys(), self.palette[:len(merged_groups)])
        }

        # Assign colors to samples
        self.sample_color: Dict[str, str] = {}
        for sample_name in samples_to_projects.keys():
            if summary_path == "project_level" or len(project_groups) == 1:
                sample_color = self.group_color[sample_name]
            else:
                sample_color = self.group_color[samples_to_projects[sample_name]]
            self.sample_color[sample_name] = sample_color

        # Copy group colors to run colors
        self.run_color = copy.deepcopy(self.group_color)
        self.palette = self.palette[len(merged_groups):]

    def _generate_plots(
        self,
        summary_path: str,
        run_data: Dict,
        sample_data: Dict,
        samples_to_projects: Dict,
        manifest_data: Dict,
        index_assignment_data: Dict,
        unassigned_sequences: Dict,
    ) -> None:
        """Generate all plots and add sections to the report."""
        # QC metrics table
        qc_metrics_function = (
            tabulate_run_stats if summary_path in ["run_level", "combined_level"] else tabulate_project_stats
        )
        self.add_run_plots(data=run_data, plot_functions=[qc_metrics_function])

        # Manifest stats
        self.add_run_plots(data=manifest_data, plot_functions=[tabulate_manifest_stats])

        # Index assignment stats
        self.add_run_plots(data=index_assignment_data, plot_functions=[tabulate_index_assignment_stats])

        # Unassigned sequences (only for run_level and combined_level)
        if summary_path in ["run_level", "combined_level"]:
            self.add_run_plots(data=unassigned_sequences, plot_functions=[tabulate_unassigned_index_stats])

        # Run-level plots
        self.add_run_plots(
            data=run_data,
            plot_functions=[plot_run_stats, plot_base_quality_hist, plot_base_quality_by_cycle],
        )

        # Sample-level plots
        self.add_sample_plots(
            data=sample_data,
            group_lookup=samples_to_projects,
            project_lookup=samples_to_projects,
        )

    def get_uuid(self):
        return str(uuid.uuid4()).replace("-", "").lower()

    def _extract_run_analysis_name(
        self,
        data: Dict[str, Any],
        source_info: str = "RunStats.json",
    ) -> Optional[str]:
        """
        Extract and validate run_analysis_name from data dict.

        Args:
            data: Dictionary containing RunName and AnalysisID keys
            source_info: Description of the data source for error messages

        Returns:
            The run_analysis_name (RunName-AnalysisID[0:4]) or None if extraction failed
        """
        run_name = data.get("RunName")
        analysis_id = data.get("AnalysisID")

        if not run_name or not analysis_id:
            log.error(
                f"Error with {source_info}. Either RunName or AnalysisID is absent.\n"
                f"RunName: {run_name}, AnalysisID: {analysis_id}\n"
                f"Please visit Elembio online documentation for more information - "
                f"https://docs.elembio.io/docs/bases2fastq/introduction/"
            )
            return None

        return f"{run_name}-{analysis_id[0:4]}"

    def _parse_run_project_data(self, data_source: str) -> List[Dict[str, Any]]:
        runs_global_data = {}
        runs_sample_data = {}
        sample_to_project = {}
        if data_source == "":
            return [runs_global_data, runs_sample_data, sample_to_project]

        for f in self.find_log_files(data_source):
            data = json.loads(f["f"])

            # Copy incomind data and reset samples to include only desired
            data_to_return = copy.deepcopy(data)
            data_to_return["SampleStats"] = []

            # get run + analysis
            run_name = data.get("RunName")
            run_analysis_name = self._extract_run_analysis_name(data, source_info=f"RunStats.json ({f['fn']})")
            if run_analysis_name is None:
                continue
            run_analysis_name = self.clean_s_name(run_analysis_name, f)

            # skip run if in user provider ignore list
            if self.is_ignore_sample(run_analysis_name):
                log.info(f"Skipping <{run_analysis_name}> because it is present in ignore list.")
                continue

            # Check run is present in the final dictionaries
            if run_analysis_name not in runs_global_data:
                runs_global_data[run_analysis_name] = data_to_return

            project = self.clean_s_name(data.get("Project", "DefaultProject"), f)

            # map sample UUIDs to run_analysis_name
            for sample_data in data["SampleStats"]:
                sample_id = sample_data["SampleID"]
                sample_name = sample_data["SampleName"]
                sample_data["RunName"] = run_name
                run_analysis_sample_name = "__".join([run_analysis_name, sample_name])

                num_polonies = sample_data["NumPolonies"]
                if num_polonies < self.min_polonies:
                    log.warning(
                        f"Skipping {run_analysis_sample_name} because it has "
                        f"<{self.min_polonies} assigned reads [n={num_polonies}]."
                    )
                    continue

                # skip run if in user provider ignore list
                if self.is_ignore_sample(sample_id) or self.is_ignore_sample(run_analysis_sample_name):
                    log.info(
                        f"Skipping <{sample_id}> ({run_analysis_sample_name}) because it is present in ignore list."
                    )
                    continue

                # If sample passes all checks add it back
                runs_sample_data[run_analysis_sample_name] = sample_data
                sample_to_project[run_analysis_sample_name] = project

            self.add_data_source(f=f, s_name=run_analysis_name, module="bases2fastq")

        return [runs_global_data, runs_sample_data, sample_to_project]

    def _parse_run_manifest(self, data_source: str) -> Dict[str, Any]:
        runs_manifest_data = {}

        if data_source == "":
            return runs_manifest_data

        for f in self.find_log_files(data_source):
            directory = f.get("root")
            if not directory:
                continue

            # Get RunName and RunID from RunStats.json
            run_stats_path = Path(directory) / "RunStats.json"
            run_stats = self._read_json_file(run_stats_path)
            if run_stats is None:
                continue

            run_analysis_name = self._extract_run_analysis_name(run_stats, source_info=str(run_stats_path))
            if run_analysis_name is None:
                continue

            run_manifest = json.loads(f["f"])
            if "Settings" not in run_manifest:
                log.warning(
                    f"<Settings> section not found in {directory}/RunManifest.json.\nSkipping RunManifest metrics."
                )
            else:
                for lane_data in run_manifest["Settings"]:
                    lane_id = lane_data.get("Lane")
                    if not lane_id:
                        log.error("<Lane> not found in Settings section of RunManifest. Skipping lanes.")
                        continue
                    lane_name = f"L{lane_id}"
                    run_lane = f"{run_analysis_name} | {lane_name}"
                    runs_manifest_data[run_lane] = {}

                    indices = []
                    indices_cycles = []
                    mask_pattern = re.compile(r"^I\d+Mask$")
                    matching_keys = [key for key in lane_data.keys() if mask_pattern.match(key)]
                    for key in matching_keys:
                        for mask_info in lane_data[key]:
                            if mask_info["Read"] not in indices:
                                indices.append(mask_info["Read"])
                            indices_cycles.append(str(len(mask_info["Cycles"])))
                    indexing = f"{' + '.join(indices_cycles)}<br>{' + '.join(indices)}"
                    runs_manifest_data[run_lane]["Indexing"] = indexing

                    runs_manifest_data[run_lane]["AdapterTrimType"] = lane_data.get("AdapterTrimType", "N/A")
                    runs_manifest_data[run_lane]["R1AdapterMinimumTrimmedLength"] = lane_data.get(
                        "R1AdapterMinimumTrimmedLength", "N/A"
                    )
                    runs_manifest_data[run_lane]["R2AdapterMinimumTrimmedLength"] = lane_data.get(
                        "R2AdapterMinimumTrimmedLength", "N/A"
                    )

            self.add_data_source(f=f, s_name=run_analysis_name, module="bases2fastq")

        return runs_manifest_data

    def _parse_run_manifest_in_project(self, data_source: str) -> Dict[str, Any]:
        project_manifest_data = {}

        if data_source == "":
            return project_manifest_data

        for f in self.find_log_files(data_source):
            directory = f.get("root")
            if not directory:
                continue

            # Get RunName and RunID from RunParameters.json
            run_manifest = Path(directory) / "../../RunManifest.json"
            project_stats = json.loads(f["f"])
            run_analysis_name = self._extract_run_analysis_name(
                project_stats, source_info=f"project RunStats.json ({f['fn']})"
            )
            if run_analysis_name is None:
                continue

            # skip run if in user provider ignore list
            if self.is_ignore_sample(run_analysis_name):
                log.info(f"Skipping <{run_analysis_name}> because it is present in ignore list.")
                continue

            run_manifest_data = self._read_json_file(run_manifest)
            if run_manifest_data is None:
                continue

            if "Settings" not in run_manifest_data:
                log.warning(f"<Settings> section not found in {run_manifest}.\nSkipping RunManifest metrics.")
            else:
                for lane_data in run_manifest_data["Settings"]:
                    lane_id = lane_data.get("Lane")
                    if not lane_id:
                        log.error("<Lane> not found in Settings section of RunManifest. Skipping lanes.")
                        continue
                    lane_name = f"L{lane_id}"
                    run_lane = f"{run_analysis_name} | {lane_name}"
                    project_manifest_data[run_lane] = {}

                    indices = []
                    indices_cycles = []
                    mask_pattern = re.compile(r"^I\d+Mask$")
                    matching_keys = [key for key in lane_data.keys() if mask_pattern.match(key)]
                    for key in matching_keys:
                        for mask_info in lane_data[key]:
                            if mask_info["Read"] not in indices:
                                indices.append(mask_info["Read"])
                            indices_cycles.append(str(len(mask_info["Cycles"])))
                    indexing = f"{' + '.join(indices_cycles)}<br>{' + '.join(indices)}"
                    project_manifest_data[run_lane]["Indexing"] = indexing

                    project_manifest_data[run_lane]["AdapterTrimType"] = lane_data.get("AdapterTrimType", "N/A")
                    project_manifest_data[run_lane]["R1AdapterMinimumTrimmedLength"] = lane_data.get(
                        "R1AdapterMinimumTrimmedLength", "N/A"
                    )
                    project_manifest_data[run_lane]["R2AdapterMinimumTrimmedLength"] = lane_data.get(
                        "R2AdapterMinimumTrimmedLength", "N/A"
                    )
            data_source_info: LoadedFileDict[Any] = {
                "fn": str(run_manifest.name),
                "root": str(run_manifest.parent),
                "sp_key": data_source,
                "s_name": str(run_manifest.with_suffix("").name),
                "f": run_manifest_data,
            }
            self.add_data_source(f=data_source_info, s_name=run_analysis_name, module="bases2fastq")

        return project_manifest_data

    def _parse_run_unassigned_sequences(self, data_source: str) -> Dict[str, Any]:
        run_unassigned_sequences = {}
        if data_source == "":
            return run_unassigned_sequences

        for f in self.find_log_files(data_source):
            data = json.loads(f["f"])

            # Get RunName and AnalysisID
            run_analysis_name = self._extract_run_analysis_name(data, source_info=f"RunStats.json ({f['fn']})")
            if run_analysis_name is None:
                continue
            run_analysis_name = self.clean_s_name(run_analysis_name, f)

            # skip run if in user provider ignore list
            if self.is_ignore_sample(run_analysis_name):
                log.info(f"Skipping <{run_analysis_name}> because it is present in ignore list.")
                continue

            # Get total polonies and build unassigned indices dictionary
            total_polonies = data.get("NumPoloniesBeforeTrimming", 0)
            if "Lanes" not in data:
                log.error(
                    f"Missing lane information in RunStats.json for run {run_analysis_name}."
                    f"Skipping building unassigned indices table."
                )
                continue
            index_number = 1
            for lane in data["Lanes"]:
                lane_id = lane.get("Lane")
                if lane_id:
                    lane_id = f"L{lane_id}"
                for sequence in lane.get("UnassignedSequences", []):
                    run_unassigned_sequences[index_number] = {
                        "Run Name": run_analysis_name,
                        "Lane": lane_id,
                        "I1": sequence["I1"],
                        "I2": sequence["I2"],
                        "Number of Polonies": sequence["Count"],
                        "% Polonies": float("nan"),
                    }
                    if total_polonies > 0:
                        run_unassigned_sequences[index_number]["% Polonies"] = round(
                            sequence["Count"] / total_polonies, 2
                        )
                    index_number += 1

        return run_unassigned_sequences

    def _parse_index_assignment(self, manifest_data_source: str) -> Dict[str, Any]:
        sample_to_index_assignment = {}

        if manifest_data_source == "":
            return sample_to_index_assignment

        for f in self.find_log_files(manifest_data_source):
            directory = f.get("root")
            if not directory:
                continue

            # Get RunName and RunID from RunStats.json
            run_stats_path = Path(directory) / "RunStats.json"
            run_stats = self._read_json_file(run_stats_path)
            if run_stats is None:
                continue

            total_polonies = 0

            # Get run name information
            run_analysis_name = self._extract_run_analysis_name(run_stats, source_info=str(run_stats_path))
            if run_analysis_name is None:
                continue

            # skip run if in user provider ignore list
            if self.is_ignore_sample(run_analysis_name):
                log.info(f"Skipping <{run_analysis_name}> because it is present in ignore list.")
                continue

            # Ensure sample stats are present
            if "SampleStats" not in run_stats:
                log.error(
                    f"Error, missing SampleStats in RunStats.json. Skipping index assignment metrics.\n"
                    f"Available keys: {list(run_stats.keys())}\n"
                    f"Please visit Elembio online documentation for more information - "
                    f"https://docs.elembio.io/docs/bases2fastq/introduction/"
                )
                continue

            # Extract per sample polony counts and overall total counts
            total_polonies = run_stats.get("NumPoloniesBeforeTrimming", 0)
            for sample_data in run_stats["SampleStats"]:
                sample_name = sample_data.get("SampleName")
                sample_id = None
                if run_analysis_name and sample_name:
                    sample_id = "__".join([run_analysis_name, sample_name])

                if "Occurrences" not in sample_data:
                    log.error(f"Missing data needed to extract index assignment for sample {sample_id}. Skipping.")
                    continue

                for occurrence in sample_data["Occurrences"]:
                    sample_expected_seq = occurrence.get("ExpectedSequence")
                    sample_counts = occurrence.get("NumPoloniesBeforeTrimming")
                    if any([element is None for element in [sample_expected_seq, sample_counts, sample_id]]):
                        log.error(
                            f"Missing data needed to extract index assignment for sample {sample_id}. Skipping."
                        )
                        continue
                    if run_analysis_name not in sample_to_index_assignment:
                        sample_to_index_assignment[run_analysis_name] = {}
                    if sample_expected_seq not in sample_to_index_assignment[run_analysis_name]:
                        sample_to_index_assignment[run_analysis_name][sample_expected_seq] = {
                            "SampleID": sample_id,
                            "SamplePolonyCounts": 0,
                            "PercentOfPolonies": float("nan"),
                            "Index1": "",
                            "Index2": "",
                        }
                    sample_to_index_assignment[run_analysis_name][sample_expected_seq]["SamplePolonyCounts"] += (
                        sample_counts
                    )

            for sample_data in sample_to_index_assignment[run_analysis_name].values():
                if total_polonies > 0:
                    sample_data["PercentOfPolonies"] = round(
                        sample_data["SamplePolonyCounts"] / total_polonies * 100, 2
                    )

            run_manifest = json.loads(f["f"])
            if "Samples" not in run_manifest:
                log.warning(
                    f"<Samples> section not found in {directory}/RunManifest.json.\n"
                    f"Skipping RunManifest sample index assignment metrics."
                )
            elif len(sample_to_index_assignment) == 0:
                log.warning("Index assignment data missing. Skipping creation of index assignment metrics.")
            else:
                for sample_data in run_manifest["Samples"]:
                    sample_name = sample_data.get("SampleName")
                    sample_id = None
                    if run_analysis_name is None or sample_name is None or "Indexes" not in sample_data:
                        continue
                    sample_id = "__".join([run_analysis_name, sample_name])
                    for index_data in sample_data["Indexes"]:
                        index_1 = index_data.get("Index1", "")
                        index_2 = index_data.get("Index2", "")
                        merged_indices = f"{index_1}{index_2}"
                        if merged_indices not in sample_to_index_assignment[run_analysis_name]:
                            log.error(f"Index assignment information not found for sample {sample_id}. Skipping.")
                            continue
                        if sample_id != sample_to_index_assignment[run_analysis_name][merged_indices]["SampleID"]:
                            log.error(
                                f"RunManifest SampleID <{sample_id}> does not match "
                                f"RunStats SampleID {sample_to_index_assignment[merged_indices]['SampleID']}."
                                "Skipping."
                            )
                            continue
                        sample_to_index_assignment[run_analysis_name][merged_indices]["Index1"] = index_1
                        sample_to_index_assignment[run_analysis_name][merged_indices]["Index2"] = index_2

        return sample_to_index_assignment

    def _parse_index_assignment_in_project(self, data_source: str) -> Dict[str, Any]:
        sample_to_index_assignment = {}

        if data_source == "":
            return sample_to_index_assignment

        for f in self.find_log_files(data_source):
            directory = f.get("root")
            if not directory:
                continue

            # Get RunManifest.json path for later use
            run_manifest = Path(directory) / "../../RunManifest.json"

            project_stats = json.loads(f["f"])
            project = self.clean_s_name(project_stats.get("Project", "DefaultProject"), f)

            run_analysis_name = self._extract_run_analysis_name(
                project_stats, source_info=f"project RunStats.json ({f['fn']})"
            )
            if run_analysis_name is None:
                continue

            # skip run if in user provider ignore list
            if self.is_ignore_sample(run_analysis_name):
                log.info(f"Skipping <{run_analysis_name}> because it is present in ignore list.")
                continue

            # Ensure sample stats are present
            if "SampleStats" not in project_stats:
                log.error(
                    f"Error, missing SampleStats in RunStats.json. Skipping index assignment metrics.\n"
                    f"Available keys: {list(project_stats.keys())}\n"
                    f"Please visit Elembio online documentation for more information - "
                    f"https://docs.elembio.io/docs/bases2fastq/introduction/"
                )
                continue

            # Extract per sample polony counts and overall total counts
            total_polonies = project_stats.get("NumPoloniesBeforeTrimming", 0)
            for sample_data in project_stats["SampleStats"]:
                sample_name = sample_data.get("SampleName")
                sample_id = None

                if run_analysis_name and sample_name:
                    sample_id = "__".join([run_analysis_name, sample_name])

                if "Occurrences" not in sample_data:
                    log.error(f"Missing data needed to extract index assignment for sample {sample_id}. Skipping.")
                    continue

                for occurrence in sample_data["Occurrences"]:
                    sample_expected_seq = occurrence.get("ExpectedSequence")
                    sample_counts = occurrence.get("NumPoloniesBeforeTrimming")
                    if any([element is None for element in [sample_expected_seq, sample_counts, sample_id]]):
                        log.error(f"Missing data needed to extract index assignment for sample {sample_id}. Skipping.")
                        continue
                    if run_analysis_name not in sample_to_index_assignment:
                        sample_to_index_assignment[run_analysis_name] = {}
                    if sample_expected_seq not in sample_to_index_assignment[run_analysis_name]:
                        sample_to_index_assignment[run_analysis_name][sample_expected_seq] = {
                            "SampleID": sample_id,
                            "Project": project,
                            "SamplePolonyCounts": 0,
                            "PercentOfPolonies": float("nan"),
                            "Index1": "",
                            "Index2": "",
                        }
                    sample_to_index_assignment[run_analysis_name][sample_expected_seq]["SamplePolonyCounts"] += (
                        sample_counts
                    )

            for sample_data in sample_to_index_assignment[run_analysis_name].values():
                if total_polonies > 0:
                    sample_data["PercentOfPolonies"] = round(
                        sample_data["SamplePolonyCounts"] / total_polonies * 100, 2
                    )

            run_manifest_data = self._read_json_file(run_manifest)
            if run_manifest_data is None:
                continue

            if "Samples" not in run_manifest_data:
                log.warning(
                    f"<Samples> section not found in {directory}/RunManifest.json.\n"
                    f"Skipping RunManifest sample index assignment metrics."
                )
            elif len(sample_to_index_assignment) == 0:
                log.warning("Index assignment data missing. Skipping creation of index assignment metrics.")
            else:
                for sample_data in run_manifest_data["Samples"]:
                    sample_name = sample_data.get("SampleName")
                    sample_id = None
                    if run_analysis_name is None or sample_name is None or "Indexes" not in sample_data:
                        continue
                    sample_id = "__".join([run_analysis_name, sample_name])
                    for index_data in sample_data["Indexes"]:
                        index_1 = index_data.get("Index1", "")
                        index_2 = index_data.get("Index2", "")
                        merged_indices = f"{index_1}{index_2}"
                        if merged_indices not in sample_to_index_assignment[run_analysis_name]:
                            continue
                        if sample_id != sample_to_index_assignment[run_analysis_name][merged_indices]["SampleID"]:
                            log.error(
                                f"RunManifest SampleID <{sample_id}> does not match "
                                f"RunStats SampleID {sample_to_index_assignment[merged_indices]['SampleID']}."
                                "Skipping."
                            )
                            continue
                        sample_to_index_assignment[run_analysis_name][merged_indices]["Index1"] = index_1
                        sample_to_index_assignment[run_analysis_name][merged_indices]["Index2"] = index_2

        return sample_to_index_assignment

    def add_run_plots(self, data, plot_functions):
        for func in plot_functions:
            plot_html, plot_name, anchor, description, helptext, plot_data = func(data, self.run_color)
            self.add_section(name=plot_name, plot=plot_html, anchor=anchor, description=description, helptext=helptext)
            self.write_data_file(plot_data, f"base2fastq:{plot_name}")

    def add_sample_plots(self, data, group_lookup, project_lookup):
        plot_functions = [
            tabulate_sample_stats,
            sequence_content_plot,
            plot_per_cycle_N_content,
            plot_adapter_content,
            plot_per_read_gc_hist,
        ]
        for func in plot_functions:
            plot_html, plot_name, anchor, description, helptext, plot_data = func(
                data, group_lookup, project_lookup, self.sample_color
            )
            self.add_section(name=plot_name, plot=plot_html, anchor=anchor, description=description, helptext=helptext)
            self.write_data_file(plot_data, f"base2fastq:{plot_name}")
