[README_EN.md](https://github.com/user-attachments/files/33908037/README_EN.md)
# PLMA-RRT

A MATLAB implementation of 3D path planning, containing PLMA-RRT itself and seven comparison algorithms, together with the three scenarios used.

## Environment

After extracting the archive, launch MATLAB from the root directory.

## Running

```matlab
% In the PLMA-RRT-Code directory
[result, scene, options] = run_planner('PLMA', 'S1', 1);
result.success
result.path
result.wall_time

% Switch the algorithm and the scenario, and give the search budget explicitly
r = run_planner('BI', 'S2', 1, ...
    struct('max_iter', 30000, 'time_limit', 120));

% Run the environment only
r = run_planner('PLMA', 'S3', 1, ...
    struct('max_iter', 100, 'time_limit', 2));
```

`run_planner(algorithm, environment, seed, overrides)` runs a single case. The defaults are `PLMA`, `S1`, `1`, and an empty struct.

## Directory structure

```text
PLMA-RRT-Code/
  README.md
  run_planner.m              # single-case entry point
  planner_options.m          # budget and per-scenario scaled parameters
  environments/
    S1.json  S2.json  S3.json
    load_environment.m
  src/
    +planners/               # the eight algorithm entry points and their support functions
    +rr/                     # tree structure, nearest neighbors, rewiring, collision checking
    plma_repaired/           # PLMA core, steering control, adaptive step size
      plma_adaptive.m
      common/
      controller/
```

## Algorithms

| ID | Algorithm | Entry point (`src/+planners`) | Native stopping condition | Source |
|---|---|---|---|---|
| `GB` | RRT* with goal bias | `gb_karaman2011.m` | search budget | [Karaman & Frazzoli, 2011](https://doi.org/10.1177/0278364911406761) |
| `BI` | BI-RRT* | `bi_jordan2013.m` | search budget | [Jordan & Perez, 2013](https://hdl.handle.net/1721.1/79884) |
| `F` | F-RRT* | `frrt_liao2021.m` | first feasible path | [Liao et al., 2021](https://doi.org/10.1016/j.eswa.2021.115457) |
| `Q` | Q-RRT* | `qrrt_jeong2019.m` | search budget | [Jeong et al., 2019](https://doi.org/10.1016/j.eswa.2019.01.032) |
| `BAI` | BAI-RRT* | `bai_xu2025.m` | search budget, followed by post-processing | [Xu et al., 2025](https://doi.org/10.1016/j.oceaneng.2025.121511) |
| `IBIP` | IBI-P-RRT* | `ibip_fan2023.m` | search budget, followed by post-processing | [Fan et al., 2023](https://doi.org/10.1016/j.eswa.2022.119137) |
| `BIAPF` | BI-APF-RRT* | `biapf_huang2024.m` | first feasible path | [Huang et al., 2024](https://doi.org/10.3390/drones8120760) |
| `PLMA` | PLMA-RRT | `plma_repaired.m` | first connection | implemented in this repository |

The code for the seven comparison algorithms was reproduced from the corresponding original papers, and the completed 3D adaptation is provided.

## Scenarios

Coordinates and lengths are all in meters. `space` is the `[min, max]` range for each of x/y/z, and `sf` is the scale factor by which the algorithm parameters are scaled.

| Scenario | Workspace (x / y / z) | Start | Goal | Number of obstacles | `sf` |
|---|---|---|---|---|---|
| S1 | [0,200] / [0,200] / [0,100] | [0,0,0] | [200,200,50] | 50 | 1 |
| S2 | [0,20] / [0,20] / [0,20] | [0,0,0] | [20,20,15] | 21 | 0.1 |
| S3 | [0,20] / [0,20] / [0,20] | [0,0,0] | [20,20,15] | 8 | 0.1 |

Each obstacle row is a center-based description, distinguished by its number of entries:

- 4 entries: sphere, `[cx, cy, cz, r]`
- 5 entries: vertical cylinder, `[cx, cy, cz, r, h]`
- 6 entries: axis-aligned box, `[cx, cy, cz, dx, dy, dz]`

Collision checking tests the intersection between the complete line segment and the closed obstacle body.

Trade-offs made while adapting each algorithm:

- **GB** is RRT*.
- **BI** implements the connect, bound, contract, and prune heuristics from the original paper; the local-consistency repairs and the parameters I chose can all be seen in the source.
- **F / Q** share the geometry and interface in `fq_support.m`. The connection radius, the neighbor radius, and F's bisection distance are scaled per scenario in `planner_options.m`; these are the values I settled on for benchmarking.
- **BAI** was originally for ship planning and has been adapted to 3D geometry. `theta_max = pi`, with empty start and goal headings, so that the hull turning constraints are disabled.
- **IBI-P** is Fan et al.'s ImprovedBi-P-RRT*. The influence radius is `max(10*sf, max inscribed sphere radius of obstacles + sqrt(3)*step)`, and the repulsive gain is `(max workspace span/20)^2`.
- **BIAPF** is implemented following Algorithms 5–9 and Figure 3 of the original paper.

## Return values

| Field | Meaning |
|---|---|
| `success` | whether the algorithm returned a path |
| `path` | final path, one `[x,y,z]` per row; empty on failure |
| `raw_path` | the searched path before smoothing/post-processing |
| `treeA`, `treeB` | the retained search trees (`X` coordinates, `p` parent indices, `c` cost to root); `treeB` is empty for single-tree algorithms |
| `iterations` | number of outer iterations |
| `first_time` | the internally recorded time to first solution; `NaN` if no first solution was found |
| `elapsed` | the algorithm's internal timing |
| `wall_time` | the full call duration measured by the runner |
| `history` | the recorded solution-improvement history |
| `details` | each algorithm's own options, stopping reason, and internal information |
