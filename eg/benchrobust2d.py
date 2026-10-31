#!/usr/bin/env python3
"""Compare two builds using benchrobust2d_r; requires only Python's standard library.

Compile the helper against each build's headers and static libqhullstatic_r.a:
  cc -O2 -I src eg/benchrobust2d_r.c build/libqhullstatic_r.a -lm -o bench
Run: taskset -c 10 python3 eg/benchrobust2d.py /path/to/baseline /path/to/patched
The helper reports CPU time for qh_new_qhull plus cleanup, excluding input parsing.
Qhull license: COPYING.txt. Oleksandr Kulkov with AI assistance, 2026-10-31.
"""
import json
import math
import random
import statistics
import subprocess
import sys

rng = random.Random(71)
n = 100_000
cloud = [(rng.random(), rng.random()) for _ in range(n)]
circle = [(math.cos(2*math.pi*i/n), math.sin(2*math.pi*i/n)) for i in range(n)]
near_line = [(x, x+rng.uniform(-1e-14, 1e-14)) for x in (rng.random() for _ in range(n))]
order_rng = random.Random(811)
rows = []
for name, points in [('random', cloud), ('circle', circle), ('near_line', near_line)]:
    data = '2\n%d\n' % n + ''.join('%.17g %.17g\n' % p for p in points)
    times = [[], []]
    vertices = [None, None]
    for block in range(20):
        order = [0, 1]
        order_rng.shuffle(order)
        for i in order:
            result = subprocess.run([sys.argv[i+1], '3'], input=data, text=True,
                                    capture_output=True, check=True)
            elapsed, count = result.stdout.split()
            vertices[i] = int(count)
            if block:
                times[i].append(float(elapsed))
    row = dict(input=name, milliseconds=[statistics.median(t) for t in times],
               vertices=vertices, samples_ms=times)
    rows.append(row)
    print(json.dumps({k:v for k,v in row.items() if k != 'samples_ms'}), flush=True)
if len(sys.argv)>3:
    with open(sys.argv[3], 'w') as out:
        json.dump(rows, out, indent=2)
