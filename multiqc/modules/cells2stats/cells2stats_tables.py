from multiqc.plots import table

from .queries import (
    get_batch_counts,
    get_batch_density,
    get_cell_count,
    get_median_cell_diameter,
    get_percent_assigned,
    get_percent_confluency,
    get_percent_mismatch,
    get_percent_nucleated_cells,
    get_total_counts,
    get_total_density,
    get_percent_assigned_spacer_polony,
    get_percent_mismatch_spacer_polony,
    get_spacer_cell_assignment_status,
    get_spacer_cell_metric_by_key
)


def merge_well_dictionaries(dict_list):
    """
    Merge metric dictionaries to build a table of well metrics
    Input is list of dictionaries of the form
    {
        "{run_name} {well_location}" :
        {
            "{metric_name}" : "{metric_value}",
            ...
        },
        ...
    }
    Output will be a merged dictionary of the form
    {
        "{run_name} {well_location}" : {
            "{metric_name}" : "{metric_value}",
            ...
        },
        ...
    }
    """
    result = {}
    for data_dict in dict_list:
        for key in data_dict:
            key_dict = result.setdefault(key, {})
            for sub_key, sub_value in data_dict[key].items():
                key_dict[sub_key] = sub_value
    return result


def tabulate_wells(c2s_run_data):
    """
    Generate a table of well metrics from the cells2stats report
    """
    plot_content = merge_well_dictionaries(
        [
            get_cell_count(c2s_run_data),
            get_percent_confluency(c2s_run_data),
            get_percent_nucleated_cells(c2s_run_data),
            get_median_cell_diameter(c2s_run_data),
            get_total_density(c2s_run_data),
            get_total_counts(c2s_run_data),
        ]
    )
    headers = {}
    headers["cell_count"] = {
        "title": "# Cells",
        "description": "The number of cells in the well",
        "min": 0,
        "scale": "GnBu",
        "format": "{d}",
    }
    headers["percent_confluency"] = {
        "title": "% Confluency",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["percent_nucleated_cells"] = {
        "title": "% Nucleated Cells",
        "description": "The percentage of cells with a segmented nucleus",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["median_cell_diameter"] = {
        "title": "Median Cell Diameter",
        "description": "Median cell diameter for cells in the well",
        "min": 0,
        "scale": "GnBu",
        "suffix": "um",
    }
    headers["total_density"] = {
        "title": "Assigned Counts K / mm2",
        "description": "Total density of assigned counts per mm2 of cell area across all batches",
        "min": 0,
        "scale": "GnBu",
    }
    headers["total_count"] = {
        "title": "Assigned Counts / Cell",
        "description": "Total average counts per cell across all batches",
        "min": 0,
        "scale": "GnBu",
    }

    pconfig = {
        "title": "cells2stats: Well QC metrics",
        "col1_header": "Run / Well",
        "id": "well_metrics_table",
        "ylab": "QC",
    }

    plot_name = "Well QC metrics table"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = "well_metrics"
    description = "Table of general well QC metrics"
    helptext = """Provides overall metrics summarizing the performance of each well"""
    return plot_html, plot_name, anchor, description, helptext, plot_content

def tabulate_spacer_wells(c2s_run_data, spacer_group_name):
    """
    Generate a table of well metrics from the cells2stats report
    """
    plot_content = merge_well_dictionaries(
        [
            get_percent_assigned_spacer_polony(c2s_run_data, spacer_group_name),
            get_percent_mismatch_spacer_polony(c2s_run_data, spacer_group_name),
            get_spacer_cell_metric_by_key(c2s_run_data, spacer_group_name, "PercentAssignedPureCells"),
            get_spacer_cell_metric_by_key(c2s_run_data, spacer_group_name, "PercentAssignedMixedCells"),
            get_spacer_cell_metric_by_key(c2s_run_data, spacer_group_name, "PercentUnassignedMixedCells"),
            get_spacer_cell_metric_by_key(c2s_run_data, spacer_group_name, "PercentUnassignedLowCountCells"),
            get_spacer_cell_metric_by_key(c2s_run_data, spacer_group_name, "AssignedCountsPerMM2", 1000),
            get_spacer_cell_metric_by_key(c2s_run_data, spacer_group_name, "MeanAssignedCountPerCell"),
            get_spacer_cell_metric_by_key(c2s_run_data, spacer_group_name, "MedianMaxSpacerCount"),
            get_spacer_cell_metric_by_key(c2s_run_data, spacer_group_name, "MeanUniqueSpacersPerCell"),
            get_spacer_cell_metric_by_key(c2s_run_data, spacer_group_name, "ExtraCellularRatio"),
            get_spacer_cell_metric_by_key(c2s_run_data, spacer_group_name, "PercentSpacerDropout"),
        ]
    )
    headers = {}
    headers["percent_assigned"] = {
        "title": "% Assigned",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["percent_mismatch"] = {
        "title": "% Mismatch",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["PercentAssignedPureCells"] = {
        "title": "% Assigned Pure Cells",
        "description": "Percentage of cells with a single spacer assigned",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["PercentAssignedMixedCells"] = {
        "title": "% Assigned Mixed Cells",
        "description": "Percentage of cells with multiple spacers assigned",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["PercentUnassignedMixedCells"] = {
        "title": "% Unassigned Mixed Cells",
        "description": "Percentage of cells with no spacer assigned, but multiple spacers in the well",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["PercentUnassignedLowCountCells"] = {
        "title": "% Unassigned Low Count Cells",
        "description": "Percentage of cells with no spacer assigned, but low counts in the well",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["AssignedCountsPerMM2"] = {
        "title": "Assigned Counts K / mm2",
        "description": "Total density of assigned spacer counts per mm2 of cell area",
        "min": 0,
        "scale": "GnBu",
        "suffix": "K",
    }
    headers["MeanAssignedCountPerCell"] = {
        "title": "Assigned Counts / Cell",
        "description": "Average assigned spacer counts per cell",
        "min": 0,
        "scale": "GnBu",
        "suffix": "",
    }
    headers["MedianMaxSpacerCount"] = {
        "title": "Median Max Spacer Count",
        "description": "Median maximum spacer count for cells in the well",
        "min": 0,
        "scale": "GnBu",
        "suffix": "",
    }
    headers["MeanUniqueSpacersPerCell"] = {
        "title": "Mean Unique Spacers / Cell",
        "description": "Mean number of unique spacers per cell in the well",
        "min": 0,
        "scale": "GnBu",
        "suffix": "",
    }
    headers["ExtraCellularRatio"] = {
        "title": "Extra Cellular Ratio",
        "description": "Ratio of extracellular spacer counts to total spacer counts in the well",
        "min": 0,
        "scale": "GnBu",
        "suffix": "",
    }
    headers["PercentSpacerDropout"] = {
        "title": "% Spacer Dropout",
        "description": "Percentage of cells with no spacer assigned, but spacer polonies present in the well",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }

    pconfig = {
        "title": f"cells2stats: {spacer_group_name} Well QC metrics",
        "col1_header": "Spacer Group / Well",
        "id": f"{spacer_group_name}_spacer_metrics_table",
        "ylab": "QC",
    }

    plot_name = f"{spacer_group_name} Well QC metrics table"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = f"{spacer_group_name}_spacer_well_metrics"
    description = "Table of general spacer well QC metrics"
    helptext = """Provides overall metrics summarizing the performance of each well"""
    return plot_html, plot_name, anchor, description, helptext, plot_content


def merge_batch_dictionaries(dict_list, metric_names):
    """
    Merge input dictionaries to build a table of batch metrics
    Input is list of dictionaries of the form
    {
        "{run_name} {well_location}" :
        {
            "{batch_name}" : "{batch_value}",
            ...
        },
        ...
    }

    Output will be a merged dictionary of the form
    {
        "{run name} {well_location} {batch_name}" : {
            "{metric_name}: "{metric_value}",
            ...
        },
        ...
    }
    """
    assert len(dict_list) == len(metric_names), "Input dictionary list and metric names must be the same length"
    result = {}
    for data_dict, metric_name in zip(dict_list, metric_names):
        for well_name in data_dict:
            well_dict = data_dict[well_name]
            for batch_name in well_dict:
                result.setdefault(f"{well_name} {batch_name}", {})[metric_name] = well_dict[batch_name]
    return result


def tabulate_batches(c2s_run_data):
    """
    Generate a table of batch metrics from the cells2stats report
    """
    plot_content = merge_batch_dictionaries(
        [
            get_batch_density(c2s_run_data),
            get_batch_counts(c2s_run_data),
            get_percent_assigned(c2s_run_data),
            get_percent_mismatch(c2s_run_data),
        ],
        ["batch_density", "batch_count", "percent_assigned", "percent_mismatch"],
    )

    headers = {}
    headers["batch_density"] = {
        "title": "Assigned Counts K / mm2",
        "description": "Assigned counts per mm2 for each batch in the well",
        "min": 0,
        "scale": "GnBu",
    }
    headers["batch_count"] = {
        "title": "Assigned Counts / Cell",
        "description": "Average assigned counts per cell for each batch in the well",
        "scale": "GnBu",
        "min": 0,
    }

    headers["percent_assigned"] = {
        "title": "% Assigned",
        "description": "Percent of polonies assigned to a barcode",
        "scale": "GnBu",
        "min": 0,
    }
    headers["percent_mismatch"] = {
        "title": "% Mismatch",
        "description": "Percent of assigned polonies assigned with a mismatch",
        "scale": "GnBu",
        "min": 0,
    }

    pconfig = {
        "title": "cells2stats: Batch QC metrics",
        "col1_header": "Run / Well / Batch",
        "id": "batch_metrics_table",
        "ylab": "QC",
    }

    plot_name = "Batch QC metrics table"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = "batch_metrics"
    description = "Table of general batch QC metrics"
    helptext = """Provides overall metrics summarizing the performance of each batch per well"""
    return plot_html, plot_name, anchor, description, helptext, plot_content
