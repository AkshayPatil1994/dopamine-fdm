"""Example: solver diagnostics from a case's run.log.

Target case: ../../examples/les_chan395/run.log (a plausible run.log path
relative to a case directory; this script has nothing to run against until
that case has actually been simulated and run.log written in this
checkout).
"""
from pathlib import Path

from dopamine_post.runlog import RunLog

case = Path(__file__).resolve().parents[2] / "examples" / "les_chan395"

log = RunLog.read(case / "run.log")
print(log)

log.plot()

log.plot(variables=["Umean", "Umax", "divergence", "CFLc"])
