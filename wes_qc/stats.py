"""
All functions used to calculate sample/variant statistics
"""

import hail as hl


def count_vars_per_cq(ht: hl.Table, cqs: list[str]) -> dict[str, int]:
    """
    Count variants per consequence type in one pass.
    Assumes ht.info.consequence is a string with consequences joined by '&'.
    """
    # Split consequences and explode into one row per consequence
    ht = ht.annotate(consequence_split=ht.info.consequence.split("&"))
    ht_exploded = ht.explode(ht.consequence_split)

    # Filter only consequences of interest
    ht_filtered = ht_exploded.filter(hl.literal(set(cqs)).contains(ht_exploded.consequence_split))

    # Group by consequence and count
    result_ht = ht_filtered.group_by(ht_filtered.consequence_split).aggregate(n=hl.agg.count())
    # Collect into Python dict
    result_dict = {row.consequence_split: row.n for row in result_ht.collect()}
    return result_dict

def print_variant_filter_stats(mt: hl.MatrixTable) -> None:
    """
    Print variant filter statistics in a table format.

    Expected labels in `mt.filters`:
        stringent_pass, stringent_fail
        medium_pass, medium_fail
        relaxed_pass, relaxed_fail
    """

    filter_levels = ["stringent", "medium", "relaxed"]
    filter_statuses = ["pass", "fail"]

    labels = [
        f"{level}_{status}"
        for level in filter_levels
        for status in filter_statuses
    ]

    # Count how many variants contain each filter label.
    label_counts = mt.aggregate_rows(
        hl.dict({
            label: hl.agg.count_where(mt.filters.contains(label))
            for label in labels
        })
    )

    print("\nVariant Filter Statistics:")
    print(f"{'Filter':<15} {'Pass':<12} {'Fail':<12} {'Total':<12}")
    print("-" * 51)

    for level in filter_levels:
        pass_count = label_counts.get(f"{level}_pass", 0)
        fail_count = label_counts.get(f"{level}_fail", 0)
        total = pass_count + fail_count

        print(f"{level:<15} {pass_count:<12} {fail_count:<12} {total:<12}")

    print()