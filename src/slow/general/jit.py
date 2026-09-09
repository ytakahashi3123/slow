#!/usr/bin/env python3

# Program to compile the face and cell loops when numba is available

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/10

import logging
import os

logger = logging.getLogger(__name__)


# numba は任意依存にしてある。無くても動く（遅いだけ）ので、
# 実行時の requirements には入れず `pip install -e ".[fast]"` で入れる。
# 環境変数 SLOW_DISABLE_NUMBA=1 で明示的に切れる。njit の中はデバッグしづらいので、
# 素の Python で追いたいときに使う。
ENV_DISABLE = 'SLOW_DISABLE_NUMBA'


def _get_njit():
  # numba が使えるなら njit を返す。使えないなら None

  if os.environ.get(ENV_DISABLE, '') not in ('', '0'):
    return None

  try:
    from numba import njit
  except ImportError:
    return None

  return njit


_njit = _get_njit()

# 呼び出し側が状況を表示できるようにしておく
FLAG_AVAILABLE = _njit is not None


def kernel(func):
  """
  Compile a numerical kernel with numba, or leave it as plain Python

  面ループ・セルループの本体をこのデコレータで包む。numba が使えないときは
  関数をそのまま返すので、挙動は変わらず速度だけが変わる。

  カーネルは numpy 配列とスカラーだけを引数に取ること（numba は Python の dict を
  受け取れない）。値の取り出しは呼び出し側のメソッドで行う。

  cache=True にしてあるので、コンパイル結果は __pycache__ に残り 2 回目以降は
  ほぼ待たない。カーネルを書き換えると作り直しになる。
  """

  if _njit is None:
    return func

  try:
    return _njit(cache=True)(func)
  except RuntimeError:
    # ソースファイルが特定できない場所（対話環境など）ではキャッシュを使えない
    return _njit(cache=False)(func)


def describe():
  # Message telling whether the kernels are compiled

  if FLAG_AVAILABLE:
    return 'numerical kernels: compiled by numba'

  if os.environ.get(ENV_DISABLE, '') not in ('', '0'):
    return f'numerical kernels: plain python ({ENV_DISABLE} is set)'

  return 'numerical kernels: plain python (numba is not installed)'
