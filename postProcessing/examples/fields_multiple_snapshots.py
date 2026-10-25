"""Work through MANY snapshots: iterate, subset by step, time-average, and read the on-the-fly averages (*_tavg.<step>).

Run:   python fields_multiple_snapshots.py [case_dir]        (default: examples/les_chan395)
Writes mean_profile_vs_time.png and stats/ in the current directory.
"""
import numpy as np
import matplotlib.pyplot as plt
from _case import case_dir
from dopamine_post.fields import FieldSeries, plot_stats
from dopamine_post.tavg import read_tavg

case = case_dir("les_chan395")
series = FieldSeries(case)
print(len(series), "snapshots, steps", series.steps())

# every snapshot in turn (one is held in memory at a time); pick a subset with series.steps()[::2]
fig, ax = plt.subplots()
for step in series.steps()[-5:]:
    ym, U = series.read(step).profile("U")
    ax.plot(U, ym, label=f"step {step}")
ax.set_xlabel("<U>"); ax.set_ylabel("y"); ax.legend()
fig.savefig("mean_profile_vs_time.png", dpi=150)

# time- and (x,z)-averaged statistics over a step range (first, last), snapshot interval in steps
steps = series.steps()
stats = series.time_average((steps[0], steps[-1]), steps[1] - steps[0])
stats.save("stats")
plot_stats(stats, case_dir=case, save="stats.png")

# averages the solver accumulated during the run (time_average output): fields/<fileout>_tavg.<step>
# t_avg, tavg = read_tavg(f"{case}/fields/<fileout>_tavg.1000")      # dict of 3-D arrays: U, "u'u'", T, ...
