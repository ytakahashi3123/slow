# Tutorial cases

各ディレクトリに `config.yml`、格子、`run_slow.sh` が入っている。
`run_slow.sh` は `slow` をインストールしていなくても動く（`src/` を `PYTHONPATH` に足す）。

```console
cd tutorial/work_nozzle_prism
./run_slow.sh                    # ログは log_slow、結果は output_result/
```

| ケース | 内容 | セル数 | 時間積分 |
|---|---|---|---|
| `work_debug_fordebug` | キメラ格子の最小ケース。動作確認用 | 25 | 非定常 BDF2 |
| `work_nozzle_prism` | 超音速ノズル内部流れ | 7,325 | 非定常 BDF2 |
| `work_sphere_chimera` | 球まわりの超音速流れ、キメラ格子 | 9,359 | 非定常 BDF2 |
| `work_wedge` | くさびの斜め衝撃波 | 20,000 | 定常 |
| `work_nozzle_prism_fast` | ノズルを 2,000 反復まで回す | 7,325 | 非定常 BDF2 |
| `work_sphere_chimera_fast` | 球を 2,000 反復まで回す | 9,359 | 非定常 BDF2 |

## Long runs with numba

`*_fast` の 2 ケースは、`numba` を入れた状態で長時間の非定常計算を回すためのもの。
物理設定は元のケース (`work_nozzle_prism` / `work_sphere_chimera`) と同一で、
**反復数 (500 --> 2,000) と出力頻度だけ**が違う。格子は元のケースのものを参照する。

`numba` は**任意依存**で、選ぶ設定項目は無い。**入っていれば面ループ・セルループが
自動でコンパイルされる**。入っていなければ素の Python で同じ計算をする（結果は同じ、遅いだけ）。

```console
pip install -e ".[fast]"         # numba を入れる
cd tutorial/work_nozzle_prism_fast
./run_slow.sh                    # 実行時間と、numba が効いているかを表示する
```

効いているかはログの先頭行でわかる。

```
--numerical kernels: compiled by numba          # 効いている
--numerical kernels: plain python               # 入っていない、または SLOW_DISABLE_NUMBA=1
```

### 実測

同一マシン、2,000 反復（内部反復は最大 25）。

| ケース | numba あり | numba なし | 倍率 |
|---|---|---|---|
| `work_nozzle_prism_fast` | **282 秒** | 約 25 時間 | 約 320 倍 |
| `work_sphere_chimera_fast` | **362 秒** | 約 34 時間 | 約 340 倍 |

「numba なし」は 5 反復の実測（ノズル 225.8 秒、球 306.2 秒）からの換算で、
通しでは測っていない。5 反復の実行時間そのものの比は 125 倍 / 144 倍だが、
これには起動と格子生成の約 1 秒が含まれるため、反復あたりの比はもっと大きい。

初回だけ numba のコンパイルに数秒かかる。`cache=True` なので 2 回目以降は
`__pycache__` から読むだけになる。

### 出力頻度

元のケースは毎反復 VTK を書く設定になっている。ノズルでは 1 回の書き出しが
約 0.13 秒で、**1 反復の計算（約 0.14 秒）と同程度**かかる。長く回すときは
出力を間引かないと、計算ではなく書き出しが実行時間を決めてしまう。
`*_fast` では結果を 50 反復ごと、リスタートを 500 反復ごとにしてある。

それでも 2,000 反復で `output_result/` は 60 MB 前後になる
（`history.csv` は内部反復ごとに 1 行なので、これだけで 8 MB 程度）。
`log_slow`、`output_result/`、`output_restart/` は `.gitignore` に入れてある。
