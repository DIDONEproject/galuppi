"""Counts how often the violas double the bass in every aria of gp.csv and
gp_multi.csv, from the MusicXML sources, and writes viola_doubling.csv. A
viola note counts as doubling when a bass note sounds at the same moment at
the unison or the octave, that is, when the two share a pitch class; viola
notes sounding over a bass rest are left out of the denominator. The feature
extraction of galuppi.qmd describes each part on its own, so the relation
between the two parts is measured here from the scores themselves.

The output is skipped when already present, so delete the CSV to rebuild it.
Requires music21 and the MusicXML corpus, whose path is set in xml_dir below.

    python3 code/discussion/viola_doubling.py
"""

import bisect
import csv
import glob
import os
import re

from music21 import converter

# The MusicXML sources, outside the repository. The tracked
# viola_doubling.csv makes galuppi.qmd independent of them.
xml_dir = os.path.expanduser("~/Documents/Corpus_DIDONE/xml")

here = os.path.dirname(os.path.abspath(__file__))
out_file = os.path.join(here, "viola_doubling.csv")


def read_arias():
    """AriaId to composer, over the binary and the five-class corpora. The
    doubtful arias are kept, with Composer NA, so that they can be placed
    against the corpus."""

    arias = {}
    for name in ("gp.csv", "gp_multi.csv"):
        with open(os.path.join(here, "..", name)) as fh:
            for row in csv.DictReader(fh):
                arias[row["AriaId"]] = row["Composer"]
    return arias


def read_scores(arias):
    """One MusicXML file per AriaId, keyed by the number of the filename."""

    scores = {}
    for path in sorted(glob.glob(os.path.join(xml_dir, "*.xml"))):
        found = re.search(r"\[(\d+)\]\.xml$", os.path.basename(path))
        if found and found.group(1) in arias:
            scores.setdefault(found.group(1), path)
    missing = set(arias) - set(scores)
    if missing:
        raise SystemExit("no MusicXML source for AriaId "
                         + ", ".join(sorted(missing)))
    return scores


def part_named(score, prefix):
    """The first part whose name begins with prefix, or None."""

    for part in score.parts:
        if (part.partName or "").strip().lower().startswith(prefix):
            return part
    return None


def doubling(path):
    """Doubling counts of one aria, or None when a part is absent."""

    score = converter.parse(path)
    viola, bass = part_named(score, "viola"), part_named(score, "bass")
    if viola is None or bass is None:
        return None

    # Bass notes as (start, end, pitch classes), for lookup by onset
    events = []
    for note in bass.flatten().notes:
        onset = float(note.offset)
        events.append((onset, onset + float(note.quarterLength),
                       {p.pitchClass for p in note.pitches}))
    events.sort()
    starts = [e[0] for e in events]

    over_bass = doubled = over_rest = 0
    for note in viola.flatten().notes:
        onset = float(note.offset)
        first = bisect.bisect_right(starts, onset) - 1
        sounding = [e for e in events[max(0, first - 3):first + 1]
                    if e[0] <= onset < e[1]]
        if not sounding:
            over_rest += 1
            continue
        over_bass += 1
        classes = {p.pitchClass for p in note.pitches}
        if any(classes & e[2] for e in sounding):
            doubled += 1

    return dict(n_viola_over_bass=over_bass, n_doubled=doubled,
                n_viola_over_rest=over_rest,
                doubling_share=doubled / over_bass if over_bass else "")


def main():

    if os.path.exists(out_file):
        print(out_file, "already present, nothing to do")
        return

    arias = read_arias()
    scores = read_scores(arias)
    rows = []
    for done, (aria_id, path) in enumerate(sorted(scores.items()), start=1):
        counts = doubling(path)
        if counts is None:
            raise SystemExit("no viola or bass part for AriaId " + aria_id)
        counts.update(AriaId=aria_id, Composer=arias[aria_id])
        rows.append(counts)
        if done % 100 == 0:
            print(done, "of", len(scores), flush=True)

    columns = ["AriaId", "Composer", "n_viola_over_bass", "n_doubled",
               "n_viola_over_rest", "doubling_share"]
    with open(out_file, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

    # Shares by composer, as reported in the paper
    for composer in sorted({r["Composer"] for r in rows}):
        subset = [r for r in rows if r["Composer"] == composer]
        doubled = sum(r["n_doubled"] for r in subset)
        total = sum(r["n_viola_over_bass"] for r in subset)
        print(f"{composer:22s} {len(subset):4d} arias  "
              f"{100 * doubled / total:.1f}% doubled")


if __name__ == "__main__":
    main()
