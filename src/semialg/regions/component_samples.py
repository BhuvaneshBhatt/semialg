from __future__ import annotations


def component_sample_points(cells_with_truth, variables):
    return [tuple(cell.sample) for cell, truth in cells_with_truth if truth]
