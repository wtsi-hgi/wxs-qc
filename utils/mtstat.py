"""Print number of samples and variants in a Hail matrix table."""
import argparse

import hail as hl

def main():
    parser = argparse.ArgumentParser(description="Print stats for a Hail matrix table.")
    parser.add_argument("mt", help="Path to the Hail matrix table")
    parser.add_argument("-s", "--samples", action="store_true", help="Print all sample IDs instead of total counts")
    args = parser.parse_args()

    hl.init()

    mt = hl.read_matrix_table(args.mt)

    if args.samples:
        samples = mt.s.collect()
        for sample in samples:
            print(sample)
    else:
        variants, n_samples = mt.count()
        print(f"Samples:  {n_samples}")
        print(f"Variants: {variants}")


if __name__ == "__main__":
    main()
