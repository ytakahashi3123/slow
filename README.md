# SLOW
Supersonic flow simulation code with low-learning cost


# Code description

`SLOW` solves compressible Navier-Stokes equations for unstructured mesh by a python script.

Note: the current implementation is **two-dimensional**. Three-dimensional support is not
completed yet (`meshdata.py` stops with a message for a 3D mesh).


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
