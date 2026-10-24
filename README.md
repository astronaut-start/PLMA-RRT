[README.md](https://github.com/user-attachments/files/33907351/README.md)
# PLMA-RRT

三维路径规划的 MATLAB 实现，包含 PLMA-RRT 本体和七种对比算法，以及三个所使用的场景。

## 环境

解压之后在根目录启动 MATLAB。
## 运行

```matlab
% 在 PLMA-RRT-Code 目录下
[result, scene, options] = run_planner('PLMA', 'S1', 1);
result.success
result.path
result.wall_time

% 更换算法和更换场景，显式给搜索预算
r = run_planner('BI', 'S2', 1, ...
    struct('max_iter', 30000, 'time_limit', 120));

% 仅运行环境
r = run_planner('PLMA', 'S3', 1, ...
    struct('max_iter', 100, 'time_limit', 2));
```

`run_planner(algorithm, environment, seed, overrides)` 单次运行一个算例，默认为 `PLMA`、`S1`、`1`、空 struct。

## 目录结构

```text
PLMA-RRT-Code/
  README.md
  run_planner.m              # 单算例入口
  planner_options.m          # 预算和按场景缩放的参数
  environments/
    S1.json  S2.json  S3.json
    load_environment.m
  src/
    +planners/               # 八个算法入口和它们的支撑函数
    +rr/                     # 树结构、近邻、重连、碰撞检测
    plma_repaired/           # PLMA 核心、方向控制、自适应步长
      plma_adaptive.m
      common/
      controller/
```

## 算法

| ID | 算法 | 入口（`src/+planners`） | 原生停止条件 | 出处 |
|---|---|---|---|---|
| `GB` | 带目标偏置的 RRT* | `gb_karaman2011.m` | 搜索预算 | [Karaman & Frazzoli, 2011](https://doi.org/10.1177/0278364911406761) |
| `BI` | BI-RRT* | `bi_jordan2013.m` | 搜索预算 | [Jordan & Perez, 2013](https://hdl.handle.net/1721.1/79884) |
| `F` | F-RRT* | `frrt_liao2021.m` | 首条可行路径 | [Liao et al., 2021](https://doi.org/10.1016/j.eswa.2021.115457) |
| `Q` | Q-RRT* | `qrrt_jeong2019.m` | 搜索预算 | [Jeong et al., 2019](https://doi.org/10.1016/j.eswa.2019.01.032) |
| `BAI` | BAI-RRT* | `bai_xu2025.m` | 搜索预算，之后做后处理 | [Xu et al., 2025](https://doi.org/10.1016/j.oceaneng.2025.121511) |
| `IBIP` | IBI-P-RRT* | `ibip_fan2023.m` | 搜索预算，之后做后处理 | [Fan et al., 2023](https://doi.org/10.1016/j.eswa.2022.119137) |
| `BIAPF` | BI-APF-RRT* | `biapf_huang2024.m` | 首条可行路径 | [Huang et al., 2024](https://doi.org/10.3390/drones8120760) |
| `PLMA` | PLMA-RRT | `plma_repaired.m` | 首次连通 | 本仓库实现 |

七个对比算法代码来源通过对应原文复现，所提供的已完成三维适配。。



## 场景

坐标和长度单位均为米。`space` 是 x/y/z 各自的 `[min, max]`，`sf` 为缩放因子，算法参数会其缩放。
| 场景 | 工作空间 (x / y / z) | 起点 | 终点 | 障碍数 | `sf` |
|---|---|---|---|---|---|
| S1 | [0,200] / [0,200] / [0,100] | [0,0,0] | [200,200,50] | 50 | 1 |
| S2 | [0,20] / [0,20] / [0,20] | [0,0,0] | [20,20,15] | 21 | 0.1 |
| S3 | [0,20] / [0,20] / [0,20] | [0,0,0] | [20,20,15] | 8 | 0.1 |

障碍物每行是以中心为基准的描述，按长度区分：

- 4 元：球，`[cx, cy, cz, r]`
- 5 元：竖直圆柱，`[cx, cy, cz, r, h]`
- 6 元：轴对齐长方体，`[cx, cy, cz, dx, dy, dz]`

碰撞检测判定为完整线段和闭合障碍体的相交。



各算法适配时做的取舍：

- **GB** 即 RRT*。
- **BI** 实现了原文的连接、边界、收缩和剪枝启发式，局部一致性修复和我选用的参数在源码里都能看到。
- **F / Q** 共用 `fq_support.m` 的几何和接口。连接半径、近邻半径和 F 的二分距离都在 `planner_options.m` 里按场景缩放，这些是我为做基准测试定的。
- **BAI** 原文为船舶规划，进行了三维几何适配。`theta_max = pi`、起终点航向置空，防止船体转向约束了。
- **IBI-P** 是 Fan 等人的 ImprovedBi-P-RRT*，影响半径取 `max(10*sf, 最大障碍内切球半径 + sqrt(3)*step)`，斥力增益取 `(最大工作空间跨度/20)^2`。
- **BIAPF** 按原文算法 5–9 和图 3 实现。
## 返回值

| 字段 | 含义 |
|---|---|
| `success` | 算法是否返回了路径 |
| `path` | 最终路径，一行一个 `[x,y,z]`，失败时为空 |
| `raw_path` | 平滑/后处理之前的搜索路径 |
| `treeA`、`treeB` | 保留的搜索树（`X` 坐标、`p` 父节点索引、`c` 到根的代价），单树算法 `treeB` 为空 |
| `iterations` | 外层迭代数 |
| `first_time` | 内部记录的首解时间，未获得首解记为 `NaN` |
| `elapsed` | 算法内部计时 |
| `wall_time` | runner 量的完整调用耗时 |
| `history` | 记录的解改进过程 |
| `details` | 各算法自己的选项、停止原因和内部信息 |



