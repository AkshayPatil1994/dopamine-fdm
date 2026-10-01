"""Example: Reynolds-stress budget statistics for a les_chan395 case.

Target case: ../../examples/les_chan395 (docs/Tools.md associates this case
with RSB stats; this script has nothing to run against until the case has
actually been simulated and stats/ populated in this checkout).
"""
from pathlib import Path

from dopamine_post.rsb import RSBStats

case = Path(__file__).resolve().parents[2] / "examples" / "les_chan395"

stats = RSBStats.read(case)
print(stats)

stats.plot()

stats.plot(last=20, out_dir="plots_last20")
