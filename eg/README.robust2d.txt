Robust 2-d hull predicates

Scope
The reentrant library uses filtered exact predicates for binary64 2-d convex
hull construction with the default precision policy. The Quickhull construction
and public qhT layout are retained. Uncertain decisions are resolved locally;
this is a single native engine, without an external retry wrapper.

Delaunay, halfspace, higher-dimensional, joggled, rerun and explicit tolerance
overrides keep their existing path. Unsupported arithmetic configurations
use the existing path too. The new arithmetic assumes ordinary IEEE binary64
rounding; fast-math/flush-to-zero modes are not covered.

The membership argument below covers static construction. Incremental and
arbitrary-start nearest-facet queries are not included in that argument.
Normals, offsets and measures are floating-point outputs. Guard failures produce
Qhull precision errors; successful output is not guaranteed for every input.

Regression
Run qhull Qt Fx with:
2
4
0 0
100000000 99999999
-100000000 100000000
99999999 99999998

Upstream omits the last point and returns3 vertices. All4 are extreme. The native
regression checks all24 permutations, extreme-scale rectangles and collinear input.
Enable it with CMake's QHULL_ENABLE_TESTING option and run ctest.

Search-membership argument
Start with a nondegenerate triangle whose original-vertex centroid is certified
strictly inside its three supporting edges. A later flipped facet raises a
precision error. Selected insertion points have a certified positive owner edge;
they need not be exactly the furthest point.

For a convex polygon K and exterior insertion v, its visible edges form one chain.
The two new supporting edges join v to tangent endpoints u,w. A point p inside
both new halfplanes is v+a(u-v)+b(w-v), with a,b>=0. For its old visible owner
support F, F(v)>0 and F(u),F(w)<=0. If F(p)>=0, then
  0 <= F(p) <= F(v)(1-a-b),
so a+b<=1 and p belongs to the triangle(v,u,w), hence the expanded hull.
This also covers exact-zero owners retained for narrow-hull bookkeeping.

Without merges, directed search tests the two adjacent new edges. A permitted
collinear merge preserves a supporting halfplane and places its nonsimplicial
survivor on the new-facet list, selecting exhaustive new-facet search. Its shared
vertex must lie between the retained endpoints. Other merge/deletion paths fail
conservatively. Thus a negative-distance local maximum cannot hide an owned point
that lies outside the new hull. Horizon search only improves a positive best;
selecting an old facet requires a strict improvement and satisfies its requeue test.
This argument depends on the certified predicates and guarded construction; it
does not assert that arbitrary distance hill-climbing finds the globally nearest
facet. A separate numerical search margin is retained.

Benchmark reproduction
Compile eg/benchrobust2d_r.c against each version's headers and sources, with the
same flags. For example, from each tree (using the helper from this branch):
  cc -O2 -DNDEBUG -I src eg/benchrobust2d_r.c src/libqhull_r/*.c -lm -o bench
Then run eg/benchrobust2d.py with paths to the baseline and patched executables.
An optional third argument saves raw JSON samples. Pin the script to one CPU for
comparison. The helper times qh_new_qhull plus cleanup, excluding input parsing;
it performs one untimed warmup before each measured batch.

Modified by Oleksandr Kulkov with AI assistance, 2026-10-31.
Qhull license: COPYING.txt. Original source: https://www.qhull.org/
