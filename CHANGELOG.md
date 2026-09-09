# Change Log
All notable changes to this project will be documented in this file.

## [Unreleased]
### Added
- `tests/` テストを新設（`pytest` で実行）
  - `tests/test_flux_jacobian.py` LU-SGS の流束ヤコビアンを厳密な微分（複素ステップ）と突き合わせる
  - `tests/test_eigenvalue.py` LU-SGS が用いる最大固有値の定義を固定する
  - `tests/test_roe_dissipation.py` 行列散逸を検証する。固有ベクトル行列 P について `P*Lambda*P^-1 == A` を確認し、
    法線方向超音速で非対角ブロック `0.5*(A-|A|)` が厳密に 0 になること（上流化）も確かめる
  - `tests/test_lusgs_options.py` 散逸の種類の選択と時間項の組み立てを固定する
  - `tests/test_rhs_snapshot.py` メッシュ読み込みから残差までを通し、基準値 (`tests/data/rhs_snapshot.npz`) との一致を確認する
- `pyproject.toml` `pip install -e .` でインストールでき、`slow` コマンドと `python3 -m slow` が使えるようになった
- `requirements.txt`, `requirements-dev.txt` 実行時／開発時の依存関係
- LU-SGS の陰解演算子に Roe 型の行列散逸を選べるようにした（`config.yml` の `kind_lusgs_dissipation`）
  - `'scalar'`（既定・従来）: `0.5*(A-lambda*I)*S`。対角はセルごとのスカラー
  - `'matrix'`: `0.5*(A-|A_Roe|)*S`。Roe 平均による行列散逸で、対角も 5x5 ブロックになる。
    項目を省略すると `'scalar'` になるので既存の `config.yml` はそのまま使える。
    陰解演算子が違っても収束後の解は一致する（デバッグ格子で相対差 1e-8 以下を確認）
  - `src/slow/time_integration/roe_dissipation.py` `|A_Roe|=P|Lambda|P^-1`。P は圧縮性 Euler 方程式の右固有ベクトル行列（3 次元形式）
  - 実測: 内部反復あたりのコストは scalar の 1.6 倍。内部反復の減衰は
    ノズル格子 (7734 要素, CFL 2.5) で 1.4--3.3 倍、デバッグ格子 (25 セル) で約 3 倍改善する
    （いずれも流れ場が発達するほど効果が大きい）。差し引きで有利

### Changed
- `src/` 以下を `src/slow/` パッケージ配下へ移動。`import` は `slow.` 起点の絶対 import に統一
  - `src/slow.py` --> `src/slow/cli.py`（`main()` を関数として呼べる形にした）
- バージョンを `src/slow/__init__.py` の `__version__` に一元化。`pyproject.toml` はこれを読む
- `tutorial/*/run_slow.sh` `python3.9` の決め打ちをやめ、インストール有無のどちらでも動くようにした
- `src/slow/time_integration/flux_jacobian.py` `lusgs_sweep.py` 内の入れ子関数 `jacobian_routine` を独立モジュールへ切り出し、外部から検証可能にした（挙動は不変）
- `src/slow/time_integration/eigenvalue.py` lambda の式を 1 か所に集約し、`lusgs_diagonal.py` と `lusgs_sweep.py` の双方から呼ぶようにした
- LU-SGS のスカラー散逸の lambda を面の左右セルの `max(lambda_a, lambda_b)` に変更し、対角項と非対角項で同一の値を使うようにした。
  日常的な条件では減衰率が 4--8% 悪化するが、現状が破綻する条件（ノズル CFL 50/200 の発達した流れ場）で 1.9--2.4 倍改善する。
  なお面上の「平均」状態にそろえる案は全条件で 0--3% 悪化したため採用していない（衝撃波近傍で散逸が不足する）
- `src/slow/time_integration/lusgs.py` 散逸の種類の選択と、時間項（定常／1 次・2 次後退差分）の組み立てを集約。
  時間項は `lusgs_diagonal.py` のスカラー版とブロック版が共有する
- `src/slow/time_integration/lusgs_sweep.py` `vecz` の 0 決め打ちをやめ法線の z 成分を用いるようにした
  （`area_vec[3]` は現状常に 0 なので数値は不変。3D 化に備えた整理）

### Fixed
- `src/slow/time_integration/lusgs_sweep.py` 流束ヤコビアンの圧力の密度微分 `prho` の符号が反転していた（正: `dp/drho=(gamma-1)*q`）。
  誤差は密度列 `J[:,0]` に限られ、大きさが動圧に比例するため高マッハ数で顕著になる。
  残差は厳密なので収束後の解は変わらないが、内部反復の収束が遅くなっていた（超音速ノズル計算で減衰率が 1.1--1.3 倍改善）

### Removed
- `src/slow/general/general.py` 未使用の `import matplotlib` を削除（未導入環境ではこれが原因で起動できなかった）
- `requirements/requirements` `requirements.txt` と `pyproject.toml` に置き換えたため削除

## [0.5.0] - 2024-01-31
### Added
- None

### Changed
- `src/gradient/gradient.py:get_slopelimiter` Primitive variables iterationa are rewritten for efficient computation. 
- `src/meshdata/meshdata.py:` オーバーラップセル判定の高速化

### Fixed
- `src/boundary/boundary.py` Key error in config is fixed.

### Removed
- None
