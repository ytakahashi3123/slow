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
  - `tests/test_lusgs_sweep.py` LU-SGS のスイープを検証する。人工格子に対して `D`, `L`, `U` を密行列として
    独立に組み立て（非対角ブロックの `A` は複素ステップ微分で作る）、スイープの結果が
    `(D+L)*D^-1*(D+U)*dq = b` を厳密に満たすことを確認する。
    セル番号・法線の向き・掃引順のいずれかを取り違えると成立しない。
    併せて掃引が黙って依存している面の並び（`face2cell_inner[0,n] < face2cell_inner[1,n]` かつ
    第 0 行が単調非減少）と、各セルの外向き法線の面積重み和が 0 になること（境界面の向きを含む）を
    凍結した 25 セル格子で固定する
  - `tests/test_convergence_check.py` 内側・外側ループの収束判定を固定する。
    相対／絶対の判定が両者で同じ向きであることも確認する
  - `tests/test_gradient_boundary.py` Green-Gauss の境界寄与を検証する。境界面が内部面より
    多い最小の格子（正方形セル 2 個、内部面 1 枚・境界面 6 枚）を組み、境界ループが動くことと、
    線形場に対して勾配が厳密に再現されることを固定する。境界の重みが距離に依らないことも確認する
  - `tests/test_rhs_snapshot.py` メッシュ読み込みから残差までを通し、基準値 (`tests/data/rhs_snapshot.npz`) との一致を確認する
- `docs/performance.md` 実行速度の実測と高速化の方針をまとめた。ルーチンごとの内訳、
  「この粒度では numpy は何も買っていない」ことの実測、採った方針（スカラーループのまま
  numba でコンパイル）と採らなかった方針（面方向のベクトル化）の比較、進捗表
- `src/slow/general/jit.py` 面ループ・セルループを numba でコンパイルする `kernel` デコレータ。
  numba は**任意依存**で、無ければ恒等デコレータになるので挙動は変わらず速度だけが変わる
  （`pip install -e ".[fast]"` で入る。`SLOW_DISABLE_NUMBA=1` で明示的に切れる）
- `src/slow/general/thermodynamics.py` カーネルから呼べる熱力学関係（音速・全エンタルピー・
  最大固有値）。`orbital` の同名メソッドおよび `time_integration/eigenvalue.py` と
  ビット一致することを `tests/test_thermodynamics.py` で固定している（各 20,000 サンプル）
- `src/slow/rhs/viscous_kernel.py` 粘性流束（応力テンソル・熱流束）と面ループを切り出した。
  入れ子関数がクロージャでやりとりしていたのを独立した関数にした
  - 実測: `flux_viscous` が **147 ms --> 0.32 ms（458 倍）**。**ビット一致**
- `src/slow/time_integration/time_integration_kernel.py` 特性時間の面ループと、
  保存変数から原始変数への変換のセルループ
  - 実測: `get_characteristic_time` が **47 ms --> 0.11 ms（424 倍）**、
    `update_primitive` が **13 ms --> 1.80 ms（7.1 倍）**。どちらも**ビット一致**
  - `update_primitive` は負値検査のセルループが Python 側に残っている（1.8 ms）
- `src/slow/time_integration/lusgs_diagonal_kernel.py` LU-SGS 対角の面ループと
  セルごとの時間項（スカラー散逸）
  - 実測: `get_diagonal` が **77 ms --> 0.17 ms（453 倍）**。**ビット一致**
- **内側反復の合計が 1,254 ms --> 4.9 ms（256 倍）**（ノズル格子 7,325 セル）
- `src/slow/rhs/advection_kernel.py` 移流流束（SLAU2 / Haenel）と面ループを切り出した。
  もとは入れ子関数がクロージャで変数をやりとりしていたのを、引数と戻り値を明示した
  独立した関数にした（`get_flux_slau2` / `get_flux_haenel` / `accumulate_advection`）。
  スキームの選択はカーネル内で文字列分岐させないよう整数の識別子で渡す
  - 実測: `slau2` が **296 ms --> 0.85 ms（348 倍）**、`haenel` が **238 ms --> 0.72 ms（330 倍）**
  - スキーム内部の二乗を `x*x` にしたため結果が 1.2e-15 だけ変わる。
    `tests/test_rhs_snapshot.py` (`rtol=1e-12`) は両スキームで通る。
    デバッグ格子の定常を収束まで流すと **収束反復数は同一**（`slau2` 205、`haenel` 489、
    `eps_muscl: 0.5` で 207）、収束解の差は 8e-16--1.5e-15
  - numba を有効にした版と素の Python 版はビット一致する
- `src/slow/time_integration/lusgs_sweep_kernel.py` LU-SGS の前進・後退スイープ
  （スカラー散逸）の面ループを素の関数として切り出した。`flux_jacobian.set_flux_jacobian` も
  カーネルから呼べるようにした（コンパイル版と素の Python 版がビット一致することを確認済み）
  - 実測: `sweep_jacobian` が **356 ms --> 0.70 ms（508 倍）**
  - 5x5 の行列ベクトル積を `np.dot` から明示的な和に置き換えたため、結果が 8e-16 だけ変わる
    （`np.dot` は 5 要素でも BLAS を使い、どの明示的な和とも一致しない）。
    スイープは残差に影響しないので、収束の振る舞いで検証した:
    デバッグ格子の定常は **205 反復で同一**（収束解の差 1.2e-15、最終残差は 7 桁一致）、
    非定常 BDF2 は **内部反復の総数 40 で同一**（`Delta Q` の差 1.4e-15）
- `src/slow/gradient/gradient_kernel.py` に minmod 制限関数と隣接セル最大・最小の
  面ループも切り出した
  - 実測: `get_slopelimiter` が **259 ms --> 0.94 ms（276 倍）**。**ビット一致**
  - numba を切った経路ではこの部分が 261 --> 320 ms と約 23% 遅くなる。
    元は 6 要素スライスに対する `np.maximum` / `np.minimum` で、そこは numpy が
    有利だった箇所（1 回の呼び出しで複数要素を畳み込める）。結果はビット一致する
- `src/slow/gradient/gradient_kernel.py` Green-Gauss 勾配の面ループを、配列とスカラーだけを
  引数に取る素の関数として切り出した。ループの形は元のままで、原始変数についての内側ループを
  明示的に書いてある（6 要素の numpy スライスはスカラー演算より遅く、numba も掛けられない）
  - 実測（ノズル格子 7,325 セル）: `get_gradient` が **122 ms --> 0.24 ms（508 倍）**。
    numba を切ると 122 ms で従来どおり。**どちらもビット一致**するので、
    スナップショット回帰もビット比較もそのまま使える
- `docs/verification.md` 数値部分の検証方法をまとめた。実際に使って有効だった 7 つの型
  （複素ステップ微分による解析的恒等式、陰解演算子を密行列に組み直しての突き合わせ、
  故意にバグを入れてテストの検出力を確かめる、時間精度の次数測定、前後のビット比較、
  幾何・離散化の恒等式、異常系の挙動確認）と、それぞれで実際に見つかった不具合、
  および踏んだ落とし穴を記録している。検証していない範囲も明記した。
  `README.md` から参照する。開発者向けメモには「変更した箇所に応じてどの方法を使うか」の
  対応表と、実務上の要点（複素ステップ微分、密行列との突き合わせ、変異テストの落とし穴、
  次数の推移の見方、ビット比較の注意）を独立した節として置いた
- `pyproject.toml` `pip install -e .` でインストールでき、`slow` コマンドと `python3 -m slow` が使えるようになった
- `src/slow/general/history.py` 残差の履歴を CSV で出力するようにした（`output_result/history.csv`）
  - 列は `iteration, iteration_inner, residual_{rho,momx,momy,momz,energy},
    deltaq_{...}`。内部反復ごとに 1 行を追記する。値は `%.17g` で書くので倍精度が丸まらない
  - `config.yml` の `post_process` に `flag_output_history` / `filename_output_history` を追加。
    どちらも省略可で、省略時は `True` / `history.csv`（既存の `config.yml` はそのまま使える）
  - リスタート計算では追記し、見出し行を重ねて書かない。`explicit_euler` は delta Q を
    持たないので `deltaq_*` は空欄になる
  - 行ごとに flush するので、異常終了しても書けた分は残る
  - 検証: ログの `Residuals` / `Delta Q` の値と CSV が完全一致することを確認。
    既存の `log_slow` / `output_restart` / VTK は 4 通りの設定でビット一致
- `src/slow/general/logging_setup.py` ログ出力を `logging` に統一した（`print` を全廃）
  - `configure_logging(level, stream, filename, fmt)` で初期化する。`cli.py:main()` の先頭で呼ぶ
  - 既定は標準出力へ素のメッセージのみ（書式 `%(message)s`）。レベル名や時刻を付けないのは、
    `run_slow.sh` のリダイレクトで得られる `log_slow` の見た目を従来と変えないため
  - 各モジュールは `logger = logging.getLogger(__name__)` を使う。`slow` 以下に集約されるので
    `logging.getLogger('slow').setLevel(...)` で一括制御できる。ライブラリとして組み込むときに
    アプリ側の設定へ干渉しないよう `propagate = False` にしてある
  - 検証: デバッグ格子 5 反復の `log_slow` から gmsh 自身の `Info` 行を除いた
    ソルバ側の出力 174 行が変更前と完全一致。`output_restart` と VTK もビット一致。
    gmsh の `Info` 行との並び順のみ変わる（従来は `print` がブロックバッファされて
    Python 側の出力が後回しになっていた。内容は同一）
- `ruff` を linter として導入した（`pyproject.toml` の `[tool.ruff]`、`requirements-dev.txt`）。
  実バグや取り違えを拾う規則だけを選び (`E4,E7,E9,F,B,PLE`)、体裁を整えるだけの規則は入れていない。
  `ruff format` は既存コードには掛けない: このコードベースは代入の桁揃えと行継続を手で整えており、
  それが読みやすさの一部なので、一括再整形すると読めない差分になる。
  `F841`（未使用ローカル変数）は、関数先頭で dict から一括展開する定型と、
  3 次元化のために残しているプレースホルダ（`wvel`, `vecz`, `flux_tmp4` など）を消させないため除外している
  - 実行: `.venv/bin/ruff check src/ tests/`
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
- セルごとに同じ式を当てるだけのループを、素の numpy の配列演算に書き換えた。
  面ループと違ってカーネルにする必要が無く、1 行で書けて読みやすい
  - `time_integration/update.py` --> `var_conserv += var_dq`（6.8 ms --> ~0）
  - `time_integration/explicit_euler.py` --> `var_conserv -= var_rhs*var_dt/volume`
  - `flowfield.set_transport_coefficients` Sutherland 則（4.2 ms --> ~0）
  - `time_integration.set_timestep`（1.3 ms --> ~0。global は `np.min`）
  - `update_primitive` の負値検査（2.0 ms --> ~0）。`np.any` で判定し、
    `np.argmax` で最初のセルを示す。メッセージも何が負なのか分かる形にし、
    終了コードを `sys.exit(1)` にした（従来の `exit()` は 0 で成功に見えていた）
  - Sutherland 則だけは結果が最大 5.5e-16（4.4% のセル）変わる。numpy の配列に対する
    `**1.5` が SIMD の `pow` を使い、スカラーの `pow` と 1 ULP 違うため。
    デバッグ格子での実際の影響は 1e-18--1e-20 で、定常を収束まで流しても
    **反復数は同一（205 回）**、収束解の差は 1.5e-18。
    他の 4 箇所はビット一致（`explicit_euler` と global 時間刻みは差が厳密に 0）
- `src/slow/time_integration/lusgs.py` 時間項の設定を 1 回だけ解決する
  `get_time_setting` を切り出した。従来は `get_time_term` がセルごとに config の辞書を
  4 回引いており、7,325 セルでは対角の計算時間の大半を占めていた。
  併せて時間刻みの正値性の確認もセルごとから 1 回（`check_timestep_array`）に変えた。
  `get_time_term` の signature と挙動は変えていない（行列散逸版とテストがそのまま使う）。
  カーネル側の式が `get_time_term` と一致することは `tests/test_lusgs_options.py` で固定
- 数値計算で使う二乗を `x**2` から `x*x` にそろえた（`orbital.get_enthalpy` /
  `get_total_energy` / `get_primitive`、`rhs/advection.py` の `q2_l` / `q2_r`）。
  `x**2` と `x*x` は倍精度で 0.09% の値について 1 ULP 違い、numba は `x**2` を
  再現できない（`x**2` / `math.pow` / `np.power` のいずれも乗算になる）ため、
  カーネルとその外側で同じ式が食い違う状態を避けるための変更
  - 一度だけ結果が変わる。実測: `tests/test_rhs_snapshot.py` (`rtol=1e-12`) は通り
    基準値の作り直しは不要。デバッグ格子の定常・`haenel` はビット一致、
    非定常 BDF2 で `restart.dat` 最大相対差 1.5e-15、ノズル格子で `Residuals` 1.1e-16。
    **定常を収束まで流すと反復数は同一（205 回）で、収束解の差は 9.4e-16**
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
- `src/slow/gradient/gradient.py:get_gradient` の境界ループが内部面用の `length`
  （`2 x num_face_inner`）を境界面番号で引いていた（正は `length_boundary`）。
  重みは `dl_n = dl_s*vcell` なので `dl_s` が約分され結果は正しかったが、
  `num_face_boundary > num_face_inner` の格子では `IndexError` で起動できなかった。
  チュートリアルの格子はどれも内部面のほうが多く、スナップショット回帰では踏めていなかった
  - 修正は挙動不変。`tests/test_rhs_snapshot.py` (`rtol=1e-12`) が通り、
    デバッグ格子 8 反復を 6 通りの設定で流して `log_slow` / `output_restart` / VTK /
    `history.csv` がビット一致することを確認済み
- 時間刻みによる 0 除算にガードを入れた。`courant_number: 0` や `timestep_outer: 0` のような
  設定で `volume/var_dt` が `inf` になり、対角が `inf`、`delta Q` が 0 になって
  「反復しても解が動かない」計算が**警告も出ずに正常終了していた**（`courant_number: 0` では
  ログに `inf` すら現れない）。`lusgs.check_timestep_positive` で `var_dt` と `timestep_outer` を、
  `time_integration.set_timestep` で全セルの `var_dt` を確かめ、原因を示して
  `sys.exit(1)` で止める（従来の `exit()` は終了コード 0 になり、呼び出し側から失敗を検知できない）
- `src/slow/time_integration/time_integration.py` 未知の `kind_time_scheme` を指定したときの
  エラー処理が `kind_time_shceme`（綴り間違い）を参照しており、意図したメッセージではなく
  `NameError` で落ちていた
- `src/slow/orbital/orbital.py:write_gmsh2vtk` 欠けている `kwargs` を拾う `try` が裸の `except` で、
  無関係な例外まで飲み込んで出力を黙って落とす可能性があった。`except KeyError` に絞った
- 定常計算の外側ループの収束判定を修正した。2 つの不具合があった
  - `src/slow/orbital/orbital.py:check_convergence_outer` の相対／絶対の分岐が逆だった。
    `flag_convergence_relative_outerloop: True` のときに残差そのものを、`False` のときに
    初期残差との比を見ていた（`check_convergence_inner` は正しい向き）。
    エネルギーの残差は 1e9 の大きさなので、既定の設定（相対・`1e-8`）では実質いつまでも判定が成立しなかった
  - `src/slow/cli.py` 判定結果 `flag_converged_outer` がどこにも使われておらず、外側ループに `break` が無かった。
    このため `criterion_convergence_outerloop` は指定しても効かず、常に `iteration_maximum` まで回っていた
  - 実測（デバッグ格子 25 セル、定常、`iteration_maximum: 60`）: 既定の `1e-8` では 60 反復で収束条件に届かず、
    修正前後で `restart.dat` はバイト一致。基準を `1e-2` に緩めると 16 反復で停止する（従来は 60 反復完走）
- 2 次後退差分 (BDF2) の立ち上げを修正した。初期条件からの計算では過去の解が `Q^0` の 1 段しかなく、
  `var_conserv_prev[0] = var_conserv_prev[1] = Q^0` のまま BDF2 の式を使っていた。このとき
  `(1.5Q^1-2Q^0+0.5Q^0)/dt = 1.5*(Q^1-Q^0)/dt` となり、実効的に刻み幅 `dt/1.5` の 1 次後退差分を解いていた。
  この `O(dt)` の誤差は最初の 1 ステップだけでも最後まで残るため、全体の時間精度が 2 次から 1 次に落ちていた。
  過去の解の段数 `num_conserv_prev_level` を持ち回し、段数が足りないステップは 1 次後退差分で立ち上げるようにした
  （リスタート時は `unsteady.dat` があれば 2 段揃うので最初のステップから BDF2 になる）
  - 実測（デバッグ格子 25 セル、`T=2.0e-3 s` まで積分、内部反復は相対 `1e-13` まで収束させ、
    `dt=2.5e-6` の解を基準とした相対 L2 誤差）:

    | `dt` [s] | 修正前 誤差 | 修正前 次数 | 修正後 誤差 | 修正後 次数 |
    | --- | --- | --- | --- | --- |
    | 2.0e-4 | 1.742e-2 |      | 9.553e-3 |      |
    | 1.0e-4 | 5.387e-3 | 1.69 | 1.840e-3 | 2.38 |
    | 5.0e-5 | 2.083e-3 | 1.37 | 3.890e-4 | 2.24 |
    | 2.5e-5 | 9.006e-4 | 1.21 | 9.149e-5 | 2.09 |

    修正前は次数が 1 に向かって劣化する。修正後は 2 に収束し、最細の `dt` で誤差が 9.8 倍小さい。
    `scalar` / `matrix` の両散逸で同じ結果になる
  - 影響範囲: `kind_steady_mode: unsteady` かつ `kind_backward_difference: 2nd_backward_diff` のみ。
    定常計算、`1st_backward_diff`、`explicit_euler` はデバッグ格子でビット一致を確認済み
- `src/slow/time_integration/lusgs_sweep.py` 流束ヤコビアンの圧力の密度微分 `prho` の符号が反転していた（正: `dp/drho=(gamma-1)*q`）。
  誤差は密度列 `J[:,0]` に限られ、大きさが動圧に比例するため高マッハ数で顕著になる。
  残差は厳密なので収束後の解は変わらないが、内部反復の収束が遅くなっていた（超音速ノズル計算で減衰率が 1.1--1.3 倍改善）

### Removed
- `src/slow/time_integration/eigenvalue.py` 最大固有値の式を
  `src/slow/general/thermodynamics.py` に一本化したため削除。
  カーネルから呼べる必要があり（`orbital.get_speedofsound` を呼んでいたため numba に
  掛けられなかった）、同じ式を 2 箇所に置く状態を避けた。
  `tests/test_eigenvalue.py` の参照先を差し替えてある
- デッドコードを削除した。いずれも実行経路から外れており、削除前後でデバッグ格子の
  `log_slow` / `output_restart` / VTK がビット一致することを 6 通りの設定で確認済み
  （scalar/matrix、定常/非定常、explicit_euler、slau2/haenel、GG+制限関数なし）
  - `src/slow/geometry/` (684 行) どこからも import されておらず、`cli.py` の呼び出しは
    コメントアウト済み。参照する `config['geometry']` セクションはどの `config.yml` にも無く、
    実行することすらできない状態だった。`cli.py` のコメントアウトされた呼び出しも削除
  - `src/slow/pending/` (240 行)、`src/slow/rhs/advection_para_pending.py` (371 行)
  - `src/slow/orbital/orbital.py:parallel_execution_decorated` と `import concurrent.futures`
    （参照は `advection.py` のコメントアウト行だけ。GIL のためスレッド並列は効かないので採らない方針）
  - `src/slow/rhs/advection.py:rhs_advection_inner` 定義されているが呼ばれていない並列化の残骸
  - 未使用の `import numpy` 4 箇所（`cli.py`, `general.py`, `explicit_euler.py`, `update.py`）と
    未使用の `import orbital` 2 箇所、`meshdata.py` の未使用の空リスト
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
