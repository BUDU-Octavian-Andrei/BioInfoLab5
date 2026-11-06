"""Exercise 2: detect tandem repeats of 6-10 bases, built from a unit of 3-5 bases."""
from typing import NamedTuple

from sequence import load_sequence

MIN_UNIT = 3
MAX_UNIT = 5
MIN_TOTAL = 6
MAX_TOTAL = 10


class Repeat(NamedTuple):
    start: int  # 0-based position in the sequence
    unit: str
    copies: int

    @property
    def total(self):
        return len(self.unit) * self.copies

    @property
    def text(self):
        return self.unit * self.copies


def is_primitive(unit):
    """False if the unit is itself a repeat, like ATAT (AT twice) or AAA (A three times)."""
    return unit not in (unit + unit)[1:-1]


def find_tandem_repeats(sequence, min_total=MIN_TOTAL, max_total=MAX_TOTAL):
    """Finds each maximal run of adjacent copies of a unit, keeping those of allowed length.

    Runs are maximal, so ATGATGATG is one repeat of ATG x 3, and not also ATG x 2 twice
    or TGA x 2. Only whole copies count towards the length of a repeat.
    """
    repeats = []
    for unit_length in range(MIN_UNIT, MAX_UNIT + 1):
        start = 0
        while start + 2 * unit_length <= len(sequence):
            # Extend for as long as each base matches the one a unit length before it.
            end = start + unit_length
            while end < len(sequence) and sequence[end] == sequence[end - unit_length]:
                end += 1

            unit = sequence[start:start + unit_length]
            copies = (end - start) // unit_length
            if copies >= 2 and is_primitive(unit) and min_total <= copies * unit_length <= max_total:
                repeats.append(Repeat(start, unit, copies))

            # Any later run must start past this point, or the run would have reached further.
            start = end - unit_length + 1

    return sorted(repeats)


def main():
    description, sequence = load_sequence()
    print(f"Sequence: {description} ({len(sequence)} nt)\n")

    repeats = find_tandem_repeats(sequence)

    print(f"{'Position':>8}  {'Unit':<6}{'Copies':>6}  {'Length':>6}  Repeat")
    for repeat in repeats:
        # Positions are printed 1-based, as NCBI numbers them.
        print(f"{repeat.start + 1:>8}  {repeat.unit:<6}{repeat.copies:>6}  {repeat.total:>6}  {repeat.text}")

    print(f"\nFound {len(repeats)} repeats of {MIN_TOTAL}-{MAX_TOTAL} bases:")
    for unit_length in range(MIN_UNIT, MAX_UNIT + 1):
        count = sum(len(repeat.unit) == unit_length for repeat in repeats)
        print(f"  unit of {unit_length} bases: {count}")

    too_long = find_tandem_repeats(sequence, MAX_TOTAL + 1, len(sequence))
    print(f"Left out {len(too_long)} repeats longer than {MAX_TOTAL} bases.")


if __name__ == "__main__":
    main()
