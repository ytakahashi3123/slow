# SLOW
Supersonic flow simulation code with low-learning cost


# Code description

`SLOW` solves compressible Navier-Stokes equations for unstructured mesh by a python script.

Note: the current implementation is **two-dimensional**. Three-dimensional support is not
completed yet (`meshdata.py` stops with a message for a 3D mesh).


# Governing equations

`SLOW` solves the compressible Navier--Stokes equations in conservative form,

$$
\frac{\partial \mathbf{Q}}{\partial t} + \nabla \cdot \left( \mathbf{F}_c - \mathbf{F}_v \right) = \mathbf{0},
\qquad
\mathbf{Q} = \left( \rho,\ \rho u,\ \rho v,\ \rho w,\ E \right)^{T}
$$

where the convective and viscous fluxes projected on a face normal $\mathbf{n}$ are

$$
\mathbf{F}_c \cdot \mathbf{n} =
\begin{pmatrix}
\rho U \\
\rho u U + p\, n_x \\
\rho v U + p\, n_y \\
\rho w U + p\, n_z \\
\rho H U
\end{pmatrix},
\qquad
\mathbf{F}_v \cdot \mathbf{n} =
\begin{pmatrix}
0 \\
\tau_{xj} n_j \\
\tau_{yj} n_j \\
\tau_{zj} n_j \\
\left( \tau_{ij} u_i + \lambda \dfrac{\partial T}{\partial x_j} \right) n_j
\end{pmatrix},
\qquad U = \mathbf{u} \cdot \mathbf{n} .
$$

The viscous stress tensor uses the Stokes hypothesis,

$$
\tau_{ij} = \mu \left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}
          - \frac{2}{3} \frac{\partial u_k}{\partial x_k} \delta_{ij} \right).
$$

The system is closed by the ideal gas law,

$$
p = \rho R T, \qquad
E = \rho \left( C_v T + \frac{1}{2} \lvert \mathbf{u} \rvert ^2 \right), \qquad
H = \frac{E + p}{\rho}, \qquad
c = \sqrt{\frac{\gamma p}{\rho}},
$$

and the transport coefficients follow Sutherland's law with a constant Prandtl number,

$$
\mu = \mu_0 \left( \frac{T}{T_0} \right)^{3/2} \frac{T_0 + C_s}{T + C_s},
\qquad
\lambda = \frac{\mu\, C_p}{Pr}.
$$

Laminar flow only; no turbulence model, source term, or chemistry is included.


# Numerical schemes

## Finite volume discretisation

A cell-centred finite volume method on an unstructured mesh. Integrating over cell $i$,

$$
V_i \frac{d \mathbf{Q}_i}{d t}
  = - \sum_{f \in \partial i} \left( \mathbf{F}_c - \mathbf{F}_v \right) \cdot \mathbf{n}_f\, S_f
  \equiv - \mathbf{R}_i ,
$$

where $S_f$ is the face area and $\mathbf{n}_f$ the outward unit normal of cell $i$.

## Gradient and reconstruction

Gradients of the primitive variables $\phi = (\rho, u, v, w, T, p)$ are obtained by the
Green--Gauss theorem, with the face value interpolated by inverse distance
($d_i$, $d_j$ are the distances from the face centre to the two cell centres),

$$
\nabla \phi_i = \frac{1}{V_i} \sum_{f \in \partial i} \phi_f\, \mathbf{n}_f\, S_f ,
\qquad
\phi_f = \frac{d_j\, \phi_i + d_i\, \phi_j}{d_i + d_j} .
$$

The left and right face states are extrapolated from the cell centres (MUSCL),

$$
\phi_L = \phi_i + \varepsilon\, \Phi_i\, d_i \left( \nabla \phi_i \cdot \mathbf{n} \right),
\qquad
\phi_R = \phi_j - \varepsilon\, \Phi_j\, d_j \left( \nabla \phi_j \cdot \mathbf{n} \right),
$$

where $\mathbf{n}$ points from cell $i$ (left) to cell $j$ (right) and
$\varepsilon \in [0, 1]$ (`eps_muscl`; $\varepsilon = 0$ gives first order). The slope
limiter $\Phi_i \in [0, 1]$ is of Barth--Jespersen type: for each face of the cell,

$$
\Phi_i = \min_{f \in \partial i}
\min \left[ 1, \ \max \left( 0, \ \frac{\Delta_i^{\pm}}{\Delta_f} \right) \right],
\qquad
\Delta_f = \left( \nabla \phi_i \cdot \mathbf{n}_f \right) \left( d_i + d_j \right),
$$

where $\Delta_i^{+} = \phi_i^{\max} - \phi_i$ if $\Delta_f \ge 0$ and
$\Delta_i^{-} = \phi_i^{\min} - \phi_i$ otherwise, the extrema being taken over the
face neighbours of cell $i$ and the cell itself.

## Convective flux

Two upwind schemes are available (`advection_scheme`). Both split the flux into a mass
flux and a pressure term,

$$
\mathbf{F}_c \cdot \mathbf{n} =
\dot{m}^{+} \boldsymbol{\psi}_L + \dot{m}^{-} \boldsymbol{\psi}_R + \tilde{p}\, \mathbf{n}_p,
\qquad
\boldsymbol{\psi} = \left( 1,\ u,\ v,\ w,\ H \right)^{T},
\qquad
\dot{m}^{\pm} = \frac{\dot{m} \pm \lvert \dot{m} \rvert}{2},
$$

where $\mathbf{n}_p = (0, n_x, n_y, n_z, 0)^T$ carries the pressure.

**SLAU2** (default) is an all-speed scheme whose mass flux is free of a cut-off Mach number:

$$
\dot{m} = \frac{1}{2} \left[ \rho_L U^{+} + \rho_R U^{-}
        - \chi\, \frac{p_R - p_L}{\bar{c}} \right],
\qquad \bar{c} = \frac{c_L + c_R}{2},
$$

$$
U^{+} = U_L + \left( 1 - g \right) \overline{\lvert U \rvert} + g \lvert U_L \rvert,
\qquad
U^{-} = U_R - \left( 1 - g \right) \overline{\lvert U \rvert} - g \lvert U_R \rvert,
$$

$$
\overline{\lvert U \rvert} = \frac{\rho_L \lvert U_L \rvert + \rho_R \lvert U_R \rvert}{\rho_L + \rho_R},
\qquad
g = - \max \left[ \min \left( M_L, 0 \right), -1 \right] \cdot
      \min \left[ \max \left( M_R, 0 \right), 1 \right],
$$

$$
\chi = \left( 1 - \bar{M} \right)^2,
\qquad
\bar{M} = \min \left( 1, \ \frac{1}{\bar{c}}
          \sqrt{ \frac{\lvert \mathbf{u}_L \rvert^2 + \lvert \mathbf{u}_R \rvert^2}{2} } \right).
$$

The pressure flux uses the van Leer type splitting functions $\beta^{\pm}$ evaluated at
$M_L$ and $M_R$,

$$
\beta^{\pm} =
\begin{cases}
\dfrac{1}{4} \left( 2 \mp M \right) \left( M \pm 1 \right)^2, & \lvert M \rvert < 1 \\
\dfrac{1}{2} \left( 1 \pm \mathrm{sign}\, M \right), & \lvert M \rvert \ge 1
\end{cases}
$$

$$
\tilde{p} = \frac{p_L + p_R}{2}
          + \frac{\beta^{+} - \beta^{-}}{2} \left( p_L - p_R \right)
          + \sqrt{ \frac{\lvert \mathbf{u}_L \rvert^2 + \lvert \mathbf{u}_R \rvert^2}{2} }
            \left( \beta^{+} + \beta^{-} - 1 \right) \frac{\rho_L + \rho_R}{2} \bar{c} .
$$

**Haenel** is the van Leer flux vector splitting in Hänel's form, in which each side uses
its own speed of sound,

$$
U^{+} = \frac{\left( U_L + c_L \right)^2}{4 c_L},
\quad
p^{+} = \frac{p_L}{4} \left( M_L + 1 \right)^2 \left( 2 - M_L \right)
\qquad \left( \lvert M_L \rvert \le 1 \right),
$$

$$
U^{-} = - \frac{\left( U_R - c_R \right)^2}{4 c_R},
\quad
p^{-} = \frac{p_R}{4} \left( M_R - 1 \right)^2 \left( 2 + M_R \right)
\qquad \left( \lvert M_R \rvert \le 1 \right),
$$

reducing to the fully upwind values for supersonic normal Mach numbers, with
$\dot{m} = \rho_L U^{+} + \rho_R U^{-}$ and $\tilde{p} = p^{+} + p^{-}$.

## Viscous flux

The face gradients are the arithmetic mean of the two adjacent cell gradients, and the
transport coefficients are averaged in the same way. No non-orthogonality correction is
applied.

## Time integration

The pseudo time step of each cell follows from the spectral radius of the flux Jacobian
including a viscous contribution,

$$
\lambda_f = \lvert U \rvert + c + \frac{2 \mu}{\rho\, d},
\qquad
\Delta \tau_i = \mathrm{CFL} \cdot \frac{V_i}{\max_{f \in \partial i} \left( \lambda_f S_f \right)} ,
$$

used either per cell (local) or as a single global minimum.

For unsteady flow, the physical time derivative is discretised by a backward difference
(BDF1 or BDF2) and the resulting nonlinear system is driven to convergence in pseudo time,

$$
V_i \frac{\Delta \mathbf{Q}_i}{\Delta \tau_i}
+ V_i \frac{3 \mathbf{Q}_i^{n+1} - 4 \mathbf{Q}_i^{n} + \mathbf{Q}_i^{n-1}}{2 \Delta t}
= - \mathbf{R}_i \left( \mathbf{Q}^{n+1} \right).
$$

Because BDF2 needs two stored levels, the very first step after a fresh start falls back
to BDF1; using BDF2 there would leave an $O(\Delta t)$ error that degrades the whole run
to first order.

The linear system of each inner iteration is solved by LU--SGS, which factorises the
implicit operator without ever forming it as a matrix,

$$
\left( \mathbf{D} + \mathbf{L} \right) \mathbf{D}^{-1}
\left( \mathbf{D} + \mathbf{U} \right) \Delta \mathbf{Q} = - \mathbf{R} - \mathbf{R}_{\mathrm{unst}},
$$

with the diagonal and the off-diagonal blocks

$$
\mathbf{D}_i = \left( \frac{V_i}{\Delta \tau_i} + \frac{3 V_i}{2 \Delta t} \right) \mathbf{I}
             + \theta \sum_{f \in \partial i} \mathbf{A}^{+}_f S_f ,
\qquad
\mathbf{A}^{\mp}_{f} = \frac{1}{2} \left( \mathbf{A} \left( \mathbf{Q}_j, \mathbf{n}_f \right)
                        \mp \mathbf{\Gamma}_f \right),
$$

where $\mathbf{A}(\mathbf{Q}, \mathbf{n}) = \partial (\mathbf{F}_c \cdot \mathbf{n}) / \partial \mathbf{Q}$
is the convective flux Jacobian, $\mathbf{n}_f$ is the outward normal of cell $i$, and
$\theta = \beta/2$ for steady flow ($\beta$ = `lusgs_beta`, over-relaxation) or $1/2$ for
unsteady flow. The dissipation $\mathbf{\Gamma}_f$ is selected by `kind_lusgs_dissipation`:

$$
\mathbf{\Gamma}_f = \lambda_f \mathbf{I}
\quad (\texttt{scalar}),
\qquad
\mathbf{\Gamma}_f = \lvert \mathbf{A}_{\mathrm{Roe}} \rvert = \mathbf{P} \lvert \mathbf{\Lambda} \rvert \mathbf{P}^{-1}
\quad (\texttt{matrix}),
$$

with $\mathbf{P}$ the matrix of right eigenvectors of $\mathbf{A}$,
$\mathbf{\Lambda} = \mathrm{diag}(U, U, U, U + c, U - c)$, and the Roe average

$$
\mathcal{R} = \sqrt{\frac{\rho_R}{\rho_L}},
\qquad
\tilde{\mathbf{u}} = \frac{\mathcal{R} \mathbf{u}_R + \mathbf{u}_L}{\mathcal{R} + 1},
\qquad
\tilde{H} = \frac{\mathcal{R} H_R + H_L}{\mathcal{R} + 1},
\qquad
\tilde{c} = \sqrt{ \left( \gamma - 1 \right)
            \left( \tilde{H} - \frac{1}{2} \lvert \tilde{\mathbf{u}} \rvert ^2 \right) } .
$$

The forward sweep solves $(\mathbf{D} + \mathbf{L}) \Delta \mathbf{Q}^{*} = \mathbf{b}$ and the
backward sweep $(\mathbf{D} + \mathbf{U}) \Delta \mathbf{Q} = \mathbf{D} \Delta \mathbf{Q}^{*}$.
See [Implicit operator (LU-SGS)](#implicit-operator-lu-sgs) for the practical differences
between the two dissipations.

An explicit Euler scheme (`explicit_euler`) is also available for reference.


# Installation

```console
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

Development (test) dependencies:

```console
pip install -r requirements-dev.txt
```

Runtime dependencies alone are listed in `requirements.txt`.


# How to start calculation

## Simulation

After `pip install -e .`:

```console
slow
```

Without installing, `src/` only needs to be on the module search path:

```console
PYTHONPATH=src python3 -m slow
```

Tutorial cases live in `tutorial/`. Each of them has a `run_slow.sh` that works either way:

```console
cd tutorial/work_nozzle_prism
./run_slow.sh
```

ケースの一覧と、`numba` を入れて長い非定常計算を回す手順は
[tutorial/README.md](tutorial/README.md) にある。


## Configuration file

Numerical simulation by `SLOW` is controlled by the configuration file: `config.yml`.
A reference copy with every key documented is `src/config.yml`.
Another file name can be given by `-file`:

```console
slow -file my_config.yml
```


# Tests

```console
pytest
```

- `tests/test_flux_jacobian.py` LU-SGS の流束ヤコビアンを厳密な微分と突き合わせる
- `tests/test_eigenvalue.py` LU-SGS が用いる最大固有値の定義を固定する
- `tests/test_roe_dissipation.py` 行列散逸 `|A_Roe|=P|Lambda|P^-1` を検証する
- `tests/test_lusgs_options.py` LU-SGS の設定（散逸の種類・時間項）と時間刻みのガードを固定する
- `tests/test_lusgs_sweep.py` LU-SGS のスイープを、密行列に組み直した `D`, `L`, `U` と突き合わせる
- `tests/test_gradient_boundary.py` Green-Gauss の境界寄与を線形場の厳密再現で検証する
- `tests/test_convergence_check.py` 内側・外側ループの収束判定を固定する
- `tests/test_rhs_snapshot.py` メッシュ読み込みから残差までを通し、基準値との一致を確認する

`tests/test_rhs_snapshot.py` は残差の値を意図的に変えたときだけ基準値を作り直す
（手順はファイル先頭のコメントを参照）。

数値部分の検証には、複素ステップ微分による解析的恒等式、陰解演算子を密行列に
組み直しての突き合わせ、時間精度の次数測定、前後のビット比較を使っている。
新しい数値を入れるとき、あるいは「挙動は変えていない」と主張するときは、
このいずれかの型で裏を取る。


## Example

![Temperature distribution around sphere at supersonic flow.\label{fig:temperature}](figure/sphere.1000.png)


## Implicit operator (LU-SGS)

Time integration by `implicit_lusgs` builds its implicit operator as
`D + L + U` with `D` on the cell and `L`, `U` from the faces. The dissipation
used for the off-diagonal part is selected by `kind_lusgs_dissipation`:

| | `scalar` (default) | `matrix` |
|---|---|---|
| 散逸の型 | スカラー（局所 Lax--Friedrichs 型） | 行列（Roe 型） |
| 非対角項 | `0.5*(A - lambda*I)*S` | `0.5*(A - \|A_Roe\|)*S` |
| `lambda` / `\|A_Roe\|` | `max(lambda_a, lambda_b)` | `P\|Lambda\|P^-1`（Roe 平均） |
| 対角 `D` | セルごとのスカラー | セルごとの 5x5 ブロック |

`matrix` は内部反復あたりのコストが `scalar` の約 1.6 倍だが、内部反復の減衰が
ノズル格子（7734 要素, CFL 2.5）で 1.4--3.3 倍速くなるため差し引きで有利になる。
陰解演算子が違っても収束後の解は一致する。項目を省略した場合は `scalar` になる。

**`matrix` を定常計算で使うときは残差の履歴 (`output_result/history.csv`) を
確認すること。** 格子と `courant_number` の組み合わせによっては、残差が下がった
あと再上昇する（chimera 格子は CFL 5.0 で反復 78 から再上昇、ノズル格子は
CFL 2.5 で 82 反復で負の状態になって停止）。一方 wedge 格子は CFL 2.5・既定の
`lusgs_beta` でも単調で `scalar` より速い。

原因は対角重みの不足なので、**`courant_number` を下げる**か **`lusgs_beta` を
上げる**かのどちらでも解消する。chimera 格子は CFL を 5.0 から 2.5 に下げるだけで
168 反復・残差 5.1e1 まで収束し、CFL 5.0 のまま `lusgs_beta` を 2.0 にした場合
(230 反復、5.9e1) より良い。ノズル格子も CFL 1.0 なら 300 反復を単調に完走する
（ただしこの格子は `scalar` でも定常では 3e6 程度で停滞するので、定常向きの設定
ではない）。`scalar` は同じ chimera 格子で CFL 20 まで再上昇しないので、両者の
差は「許容できる CFL の幅」として現れる。

定常では内部反復が 1 回に固定されるので、1 回の掃引の質がそのまま解の質になり、
対角の重みが足りないと発散する。`scalar` の対角は `lambda=|u.n|+c+2mu/(rho*d)` を使うため
すべての固有値の絶対値を過大に見積もっており、その余剰が 1 回の掃引を安定化させている。
`matrix` の対角は正確な `A+` なので、その余剰を `lusgs_beta` で補うか、
`courant_number` を下げて時間項 `V/dt` に補わせる必要がある。
非定常計算では対角が時間項に支配され内部反復も複数回あるため影響しない。
tutorial のうち定常なのは `work_wedge` だけで、そこでは `matrix` が有利。


## Performance

面ループ・セルループは `*_kernel.py` に切り出してあり、`numba` が入っていれば
コンパイルされる。`numba` は**任意依存**で、無くても動く（遅いだけ）。

```console
pip install -e ".[fast]"      # numba を入れる
SLOW_DISABLE_NUMBA=1 slow     # 明示的に切る（デバッグ時）
```

実測（ノズル格子 7,325 セル）: 内側反復 1 回が 1,254 ms --> 4.9 ms、
メッシュ生成が 0.47 --> 0.22 秒、外側 5 反復の実行全体で 166.8 秒 --> 1.9 秒。
結果は変わらない（差は 1e-15 程度）。
`numba` を入れない場合は 180.3 秒で、導入前より 8% 遅い。


## Requirements

- python (version >= 3.9)
- numpy (version >= 1.22)
- pyyaml (version >= 6.0)
- gmsh (version >= 4.9.5)

Optional:

- numba (version >= 0.59) 面ループ・セルループのコンパイルに使う


# Contact:

Yusuke Takahashi, Hokkaido University

ytakahashi@eng.hokudai.ac.jp


# References
