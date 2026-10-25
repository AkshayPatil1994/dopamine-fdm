"""ParaView time series of chunked snapshots (chunked_snapshot = .true. writes the box <fileout>Block.<step>).

Run:   python xmf_box.py <case_dir>
Writes <case>/paraview/box.xmf; fields are on the box's cell centres, staggered velocities take the lower face.
The time-average boxes (<fileout>_tavgBlock.<step>) are not part of the series.
CLI equivalent:  dopamine-post fields box-xmf --case-dir <case>
"""
import sys
from dopamine_post.fields import read_box_snapshot, write_box_xmf

if len(sys.argv) < 2:
    sys.exit(__doc__)
case = sys.argv[1]
write_box_xmf(case)
# one box on its own: names of the fields it holds and the cell range it covers (f1, f2)
# box = read_box_snapshot(f"{case}/fields/<fileout>Block.1000"); print(box["f1"], box["f2"], list(box["fields"]))
