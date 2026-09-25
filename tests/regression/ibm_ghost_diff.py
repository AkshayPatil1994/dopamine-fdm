#!/usr/bin/env python3
"""Compare the IBM ghost-cell lists (ibm_ghosts.rank* written under DOPAMINE_TRACE_DIR) of two runs by global index.

  ibm_ghost_diff.py <trace_dir_a> <trace_dir_b>

Reports, per velocity component, cells listed in only one run and cells listed by more than one rank (a ghost cell
shared between ranks gets its boundary condition applied more than once, each time from different halo data).
"""
import collections
import glob
import sys


def load(d):
    cells = collections.defaultdict(list)
    for f in glob.glob(d + '/ibm_ghosts.rank*'):
        rank = int(f.rsplit('rank', 1)[1])
        for line in open(f):
            c, i, j, k = line.split()
            cells[(c, int(i), int(j), int(k))].append(rank)
    return cells


a, b = load(sys.argv[1]), load(sys.argv[2])
only_a, only_b = sorted(set(a) - set(b)), sorted(set(b) - set(a))
dup_a = sorted(c for c, r in a.items() if len(r) > 1)
dup_b = sorted((c, r) for c, r in b.items() if len(r) > 1)
print(f'cells: A {len(a)}  B {len(b)};  only in A: {len(only_a)}  only in B: {len(only_b)};  '
      f'listed by several ranks: A {len(dup_a)}  B {len(dup_b)}')
for name, lst in (('only in A', only_a), ('only in B', only_b)):
    for c in lst[:6]:
        print(f'  {name}: {c}')
for c, r in dup_b[:6]:
    print(f'  shared in B: {c} ranks {r}')
sys.exit(1 if (only_a or only_b or dup_a or dup_b) else 0)
