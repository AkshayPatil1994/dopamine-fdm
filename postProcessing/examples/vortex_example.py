from dopamine_post.fields import FieldSeries
from dopamine_post.vortex import q_criterion, q_volume_fraction, animate_q_isosurface

case = "."

series = FieldSeries(case)

snap = series.latest()
Q = q_criterion(snap)
print(f"step {snap.step}: Q in [{Q.min():.3g}, {Q.max():.3g}], "
      f"volume fraction Q>5: {q_volume_fraction(snap, 5.0):.4f}")

animate_q_isosurface(series, "q.mp4", qval=5.0, start=1000, stride=2)
