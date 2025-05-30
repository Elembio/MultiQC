from .utils import json_decode_int, json_decode_float, is_nan, summarize_batch_names, find_entry

small_value = 0.00000001

def get_cell_count(c2s_run_data):
    """
    Get the cell count for each well in the run"""
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            result.setdefault(f"{run_name} {well_data['WellLocation']}", {})["cell_count"] = json_decode_int(
                well_data.get("CellCount", 0)
            )
    return result


def get_percent_confluency(c2s_run_data):
    """
    Get the percent confluency for each well in the run
    """
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            val = json_decode_float(well_data.get("PercentConfluency", float("nan")))
            if not is_nan(val):
                result.setdefault(f"{run_name} {well_data['WellLocation']}", {})["percent_confluency"] = val
    return result


def get_percent_nucleated_cells(c2s_run_data):
    """
    Get the percent nucleated cells for each well in the run
    """
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            val = json_decode_float(well_data.get("PercentNucleatedCells", float("nan")))
            if not is_nan(val):
                result.setdefault(f"{run_name} {well_data['WellLocation']}", {})["percent_nucleated_cells"] = val
    return result


def get_median_cell_diameter(c2s_run_data):
    """
    Get the median cell diameter for each well in the run
    """
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            val = json_decode_float(well_data.get("MedianCellDiameter", float("nan")))
            if not is_nan(val):
                result.setdefault(f"{run_name} {well_data['WellLocation']}", {})["median_cell_diameter"] = val
    return result


def get_total_density(c2s_run_data):
    """
    Get the total density for each well in the run
    """
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            well_location = well_data.get("WellLocation", "")
            if well_location != "":
                val = json_decode_float(well_data.get("AssignedCountsPerMM2", float("nan")))
                if not is_nan(val):
                    result.setdefault(f"{run_name} {well_location}", {})["total_density"] = val / 1000.0
    return result


def get_total_counts(c2s_run_data):
    """
    Get the average total counts per cell for each well in the run
    """
    batch_names = summarize_batch_names(c2s_run_data)
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            well_location = well_data.get("WellLocation", "")
            if well_location != "":
                total_count = 0
                for batch_name in batch_names:
                    batch_data = find_entry(well_data.get("Batches", []), "BatchName", batch_name, {})
                    val = json_decode_float(batch_data.get("MeanAssignedCountPerCell", float("nan")))
                    if not is_nan(val):
                        total_count += val
                result.setdefault(f"{run_name} {well_location}", {})["total_count"] = total_count
    return result


def get_batch_density(c2s_run_data):
    """
    Get the average density for each batch in each well in the run"""
    batch_names = summarize_batch_names(c2s_run_data)
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            well_location = well_data.get("WellLocation", "")
            if well_location != "":
                for batch_name in batch_names:
                    batch_data = find_entry(well_data.get("Batches", []), "BatchName", batch_name, {})
                    val = json_decode_float(batch_data.get("AssignedCountsPerMM2", float("nan")))
                    if not is_nan(val):
                        result.setdefault(f"{run_name} {well_location}", {})[batch_name] = val / 1000.0
    return result


def get_batch_counts(c2s_run_data):
    """
    Get the average counts per cell for each batch in each well in the run
    """
    batch_names = summarize_batch_names(c2s_run_data)
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            well_location = well_data.get("WellLocation", "")
            if well_location != "":
                for batch_name in batch_names:
                    batch_data = find_entry(well_data.get("Batches", []), "BatchName", batch_name, {})
                    val = json_decode_float(batch_data.get("MeanAssignedCountPerCell", float("nan")))
                    if not is_nan(val):
                        result.setdefault(f"{run_name} {well_location}", {})[batch_name] = val
    return result


def get_percent_assigned(c2s_run_data):
    """
    Get the percent assigned reads for each batch in each well in the run
    """
    batch_names = summarize_batch_names(c2s_run_data)
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for batch_name in batch_names:
            batch_entry = find_entry(run_data.get("DemuxStats", {}).get("Batches", []), "BatchName", batch_name, {})
            for well_data in batch_entry.get("Wells", []):
                well_location = well_data.get("WellLocation", "")
                if well_location != "":
                    val = json_decode_float(well_data.get("PercentAssignedReads", float("nan")))
                    if not is_nan(val):
                        result.setdefault(f"{run_name} {well_location}", {})[batch_name] = val
    return result


def get_percent_mismatch(c2s_run_data):
    """
    Get the percent mismatch for each batch in each well in the run
    """
    batch_names = summarize_batch_names(c2s_run_data)
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for batch_name in batch_names:
            batch_entry = find_entry(run_data.get("DemuxStats", {}).get("Batches", []), "BatchName", batch_name, {})
            for well_data in batch_entry.get("Wells", []):
                well_location = well_data.get("WellLocation", "")
                if well_location != "":
                    val = json_decode_float(well_data.get("PercentMismatch", float("nan")))
                    if not is_nan(val):
                        result.setdefault(f"{run_name} {well_location}", {})[batch_name] = val
    return result

def get_percent_assigned_spacer_polony(c2s_run_data, spacer_group_name):
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for batch_data in run_data.get("SpacerStats", {}).get("Batches", []):
            for well_data in batch_data.get("Wells", []):
                for spacer_group_data in well_data.get("SpacerGroups", []):
                    if spacer_group_data.get("GroupName", "") == spacer_group_name:
                        value = json_decode_float(spacer_group_data.get("PolonyAssignmentStats", {}).get("PercentAssigned", float("nan")))
                        if not is_nan(value) and value > 0:
                            result.setdefault(f"{run_name} {well_data['WellLocation']}", {})["percent_assigned"] = value
                        else:
                            result.setdefault(f"{run_name} {well_data['WellLocation']}", {})["percent_assigned"] = small_value
    return result

def get_percent_mismatch_spacer_polony(c2s_run_data, spacer_group_name):
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for batch_data in run_data.get("SpacerStats", {}).get("Batches", []):
            for well_data in batch_data.get("Wells", []):
                for spacer_group_data in well_data.get("SpacerGroups", []):
                    if spacer_group_data.get("GroupName", "") == spacer_group_name:
                        value = json_decode_float(spacer_group_data.get("PolonyAssignmentStats", {}).get("PercentMismatch", float("nan")))
                        if not is_nan(value) and value > 0:
                            result.setdefault(f"{run_name} {well_data['WellLocation']}", {})["percent_mismatch"] = value
                        else:
                            result.setdefault(f"{run_name} {well_data['WellLocation']}", {})["percent_mismatch"] = small_value

    return result



def get_spacer_cell_assignment_status(c2s_run_data, spacer_group_name):
    """
    spacer cell assignment"""
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for batch_data in run_data.get("SpacerStats", {}).get("Batches", []):
            for well_data in batch_data.get("Wells", []):
                for spacer_group_data in well_data.get("SpacerGroups", []):
                    if spacer_group_data.get("GroupName", "") == spacer_group_name:
                        cell_assignment_data = spacer_group_data.get("CellAssignmentStats", {})
                        for key in ["PercentAssignedPureCells", "PercentAssignedMixedCells", "PercentUnassignedMixedCells", "PercentUnassignedLowCountCells"]:
                            val = json_decode_float(cell_assignment_data.get(key, float("nan")))
                            if not is_nan(val) and val > 0:
                                result.setdefault(f"{run_name} {well_data['WellLocation']}", {})[key] = val
                            else:
                                result.setdefault(f"{run_name} {well_data['WellLocation']}", {})[key] = small_value
    return result

def get_spacer_cell_metric_by_key(c2s_run_data, spacer_group_name, key, factor = 1.0):
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for batch_data in run_data.get("SpacerStats", {}).get("Batches", []):
            for well_data in batch_data.get("Wells", []):
                for spacer_group_data in well_data.get("SpacerGroups", []):
                    if spacer_group_data.get("GroupName", "") == spacer_group_name:
                        value = json_decode_float(spacer_group_data.get("CellAssignmentStats", {}).get(key, float("nan")))
                        if not is_nan(value) and value > 0:
                            result.setdefault(f"{run_name} {well_data['WellLocation']}", {})[key] = value / factor
                        else:
                            result.setdefault(f"{run_name} {well_data['WellLocation']}", {})[key] = small_value
    return result



    # "SpacerStats": {
    #     "Batches": [
    #         {
    #             "BatchName": "B02",
    #             "SpacerGroups": [
    #                 {
    #                     "GroupName": "MCF-7",
    #                     "CellAssignmentStats": {
    #                         "PercentAssignedPureCells": 3.048121425997949,
    #                         "PercentAssignedMixedCells": 2.4133848205185053,
    #                         "PercentUnassignedMixedCells": 0.7312750776772874,
    #                         "PercentUnassignedLowCountCells": 1.0950264980399764,
    #                         "AssignedCountsPerMM2": 0.0003228037689836166,
    #                         "MeanAssignedCountPerCell": 0.767860387200877,
    #                         "MedianMaxSpacerCount": 9.0,
    #                         "MeanUniqueSpacersPerCell": 1.6691399824706743,
    #                         "ExtraCellularRatio": 0.0039840009036781095,
    #                         "PercentSpacerDropout": 5.877268798617113
    #                     },
    #                     "PolonyAssignmentStats": {
    #                         "PercentAssigned": 30.04947425306522,
    #                         "PercentMismatch": 84.41735732106471,
    #                         "UnassignedSequences": [
