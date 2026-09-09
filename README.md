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
- `tests/test_lusgs_options.py` LU-SGS の設定（散逸の種類・時間項）を固定する
- `tests/test_rhs_snapshot.py` メッシュ読み込みから残差までを通し、基準値との一致を確認する

`tests/test_rhs_snapshot.py` は残差の値を意図的に変えたときだけ基準値を作り直す
（手順はファイル先頭のコメントを参照）。


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


## Requirements

- python (version >= 3.9)
- numpy (version >= 1.22)
- pyyaml (version >= 6.0)
- gmsh (version >= 4.9.5)


# Contact:

Yusuke Takahashi, Hokkaido University

ytakahashi@eng.hokudai.ac.jp


# References
