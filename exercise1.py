"""Exercise 1: simulate gel electrophoresis of 10 random samples from a DNA sequence."""
import argparse
import math
import random
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle

from sequence import load_sequence

SAMPLE_COUNT = 10
MIN_SAMPLE_LENGTH = 100
MAX_SAMPLE_LENGTH = 3000

# Reference fragments, so the sizes of the sample bands can be read off the gel.
LADDER = [100, 200, 300, 500, 1000, 1500, 2000, 3000]

# Fragments of this size would be too large to enter the gel, so their speed is 0.
EXCLUSION_SIZE = 10_000
# The gel runs from 0 (wells) to 1 (bottom). At this time a 100 bp band is at 0.9.
RUN_TIME = 0.45
FRAMES = 60

OUTPUT_DIR = Path(__file__).parent / "output"
GEL_COLOR = "#11161c"
WELL_COLOR = "#4a5560"
BAND_COLOR = "#f2f7fb"
LABEL_COLOR = "#b8c4ce"
BAND_WIDTH = 0.7
BAND_HEIGHT = 0.012


def take_samples(sequence):
    """Cuts SAMPLE_COUNT random contiguous fragments out of the sequence."""
    # A sample cannot be longer than the sequence it is cut from.
    longest = min(MAX_SAMPLE_LENGTH, len(sequence))
    samples = []
    for _ in range(SAMPLE_COUNT):
        length = random.randint(MIN_SAMPLE_LENGTH, longest)
        start = random.randint(0, len(sequence) - length)
        samples.append(sequence[start:start + length])
    return samples


def migration_speed(size):
    """Longer fragments feel more friction: speed falls linearly with log10(size)."""
    return math.log10(EXCLUSION_SIZE / size)


def migration_distance(size, time):
    return migration_speed(size) * time


def animate_gel(sizes):
    """Builds the animated gel: lane 0 is the ladder, lanes 1-10 hold one sample each."""
    lanes = [(0, size, False) for size in LADDER]
    lanes += [(lane, size, True) for lane, size in enumerate(sizes, start=1)]

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_facecolor(GEL_COLOR)
    ax.set_xlim(-0.75, len(sizes) + 0.75)
    ax.set_ylim(1, -0.04)
    ax.set_xticks(range(len(sizes) + 1), ["Ladder"] + [str(i) for i in range(1, len(sizes) + 1)])
    ax.set_yticks([])
    ax.tick_params(axis="x", length=0)
    ax.set_xlabel("Lane")
    ax.set_ylabel("Migration direction (wells at the top)")

    for lane in range(len(sizes) + 1):
        ax.add_patch(Rectangle((lane - BAND_WIDTH / 2, -0.03), BAND_WIDTH, 0.02, color=WELL_COLOR))

    bands = []
    for lane, size, _ in lanes:
        band = Rectangle((lane - BAND_WIDTH / 2, 0), BAND_WIDTH, BAND_HEIGHT, color=BAND_COLOR)
        label = ax.text(lane, 0, f"{size} bp", color=LABEL_COLOR, fontsize=8, ha="center", va="top")
        ax.add_patch(band)
        bands.append((band, label, size))

    def update(frame):
        time = RUN_TIME * frame / (FRAMES - 1)
        for band, label, size in bands:
            distance = migration_distance(size, time)
            band.set_y(distance)
            label.set_y(distance + BAND_HEIGHT + 0.004)
        ax.set_title(f"Gel electrophoresis simulation: {100 * frame / (FRAMES - 1):.0f}% of run time")
        return []

    update(0)
    return fig, FuncAnimation(fig, update, frames=FRAMES, interval=60, repeat=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, help="repeat the same random samples")
    args = parser.parse_args()
    random.seed(args.seed)

    description, sequence = load_sequence()
    print(f"Sequence: {description} ({len(sequence)} nt)\n")

    samples = take_samples(sequence)
    sizes = [len(sample) for sample in samples]

    print("Lane  Size (bp)  Distance (fraction of gel)")
    for lane, size in enumerate(sizes, start=1):
        print(f"{lane:>4}  {size:>9}  {migration_distance(size, RUN_TIME):.3f}")

    print("\nOrder of the bands, furthest travelled first:")
    print(" > ".join(f"{size} bp" for size in sorted(sizes)))

    fig, animation = animate_gel(sizes)
    OUTPUT_DIR.mkdir(exist_ok=True)
    animation.save(OUTPUT_DIR / "gel.gif", writer="pillow", fps=15)
    fig.savefig(OUTPUT_DIR / "gel.png", dpi=150, bbox_inches="tight")
    print(f"\nSaved {OUTPUT_DIR / 'gel.png'} and {OUTPUT_DIR / 'gel.gif'}")

    if plt.get_backend().lower() != "agg":
        plt.show()


if __name__ == "__main__":
    main()
