from __future__ import annotations


def find_region_instance(cells_with_truth):
    for cell, truth in cells_with_truth:
        if truth:
            return dict(enumerate(cell.sample))
    return None


def component_instances(cells_with_truth, variables):
    out = []
    for cell, truth in cells_with_truth:
        if truth:
            out.append(
                {variable: sample for variable, sample in zip(variables, cell.sample, strict=True)}
            )
    return out
