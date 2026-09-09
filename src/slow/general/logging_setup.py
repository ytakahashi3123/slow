#!/usr/bin/env python3

# Program to configure logging for the solver

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

import logging
import sys


# ソルバ全体のログはこの名前の下にまとめる（各モジュールは
# logging.getLogger(__name__) を使うので slow.* が全て子になる）
LOGGER_NAME = 'slow'

# 既定の書式は素のメッセージだけ。従来の print による log_slow と
# 1 バイトも変えないための選択で、レベル名や時刻は付けない
DEFAULT_FORMAT = '%(message)s'


def configure_logging(level=logging.INFO, stream=None, filename=None, fmt=DEFAULT_FORMAT):
  """
  Configure the logger used by every module of the solver

  ソルバは run_slow.sh から標準出力をリダイレクトして使う運用なので、
  既定では標準出力へ素のメッセージだけを出す。従来の print と同じ見た目になる。

  --level: 出力するレベル。ログを絞りたいときは logging.WARNING などを渡す
  --stream: 出力先ストリーム。既定は sys.stdout
  --filename: 指定するとファイルにも書き出す（標準出力への出力は残る）
  --fmt: logging の書式文字列

  main() から複数回呼ばれても handler が重複しないよう、毎回作り直す。
  """

  logger = logging.getLogger(LOGGER_NAME)
  logger.setLevel(level)
  # 呼び出し元のアプリ側の設定に干渉しないよう、上位へは伝播させない
  logger.propagate = False

  for handler_tmp in list(logger.handlers):
    logger.removeHandler(handler_tmp)
    handler_tmp.close()

  formatter = logging.Formatter(fmt)

  handler = logging.StreamHandler(sys.stdout if stream is None else stream)
  handler.setFormatter(formatter)
  logger.addHandler(handler)

  if filename is not None:
    handler_file = logging.FileHandler(filename, mode='w')
    handler_file.setFormatter(formatter)
    logger.addHandler(handler_file)

  return logger
