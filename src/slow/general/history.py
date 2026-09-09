#!/usr/bin/env python3

# Program to record the residual history in a machine-readable file (CSV)

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

import logging
import os

logger = logging.getLogger(__name__)


# 設定を省略した場合の既定値（既存の config.yml をそのまま使えるように）
FLAG_OUTPUT_DEFAULT = True
FILENAME_DEFAULT    = 'history.csv'

# 保存変数の並びに対応する列名。conservative variables の順序に合わせる
NAME_CONSERV = ('rho', 'momx', 'momy', 'momz', 'energy')


def get_history_setting(config):
  """
  Settings of the residual history, given by the control file

  --戻り値: (出力するか, ファイルのパス)
  """

  setting = config['post_process']

  flag_output = bool( setting.get('flag_output_history', FLAG_OUTPUT_DEFAULT) )
  filename    = str( setting.get('filename_output_history', FILENAME_DEFAULT) )
  directory   = setting['directory_output']

  return flag_output, os.path.join(directory, filename)


def open_history(config, dimension_dict):
  """
  Open the residual-history file and write its header

  ログは人が読む形でしか残らないので、収束の履歴を機械的に扱えるよう
  CSV でも書き出す。内部反復ごとに 1 行を追加する。

  リスタート計算では追記し、見出し行を重ねて書かない。
  出力しない設定のときは None を返し、呼び出し側はそのまま渡せばよい。
  """

  flag_output, filename = get_history_setting(config)
  if not flag_output :
    return None

  # 初期条件からの計算では作り直し、リスタートでは追記する
  flag_initial = config['computational_setup']['flag_initial']
  flag_header  = flag_initial or not os.path.exists(filename)

  logger.info('Writing residual history: %s', filename)

  # 途中で異常終了しても書けた分が残るよう、行ごとに flush する
  file = open(filename, 'w' if flag_initial else 'a', buffering=1)

  if flag_header :
    num_conserv = dimension_dict['num_conservative']
    columns = ['iteration', 'iteration_inner']
    columns += [ 'residual_'+NAME_CONSERV[m] for m in range(0,num_conserv) ]
    columns += [ 'deltaq_'+NAME_CONSERV[m]   for m in range(0,num_conserv) ]
    file.write( ','.join(columns) + '\n' )

  return file


def write_history(file, iteration, iteration_inner, sum_rhs, sum_dq):
  """
  Append one row of the residual history

  --sum_rhs: display_residual が返す残差の二乗和（保存変数ごと）
  --sum_dq:  display_deltaq が返す delta Q の二乗和（保存変数ごと）。
             explicit_euler では delta Q を持たないので None でよい
  """

  if file is None :
    return

  values = [ str(iteration), str(iteration_inner) ]
  # 倍精度を丸めずに残す（%.17g は往復して同じ値に戻る）
  values += [ '%.17g'%value for value in sum_rhs ]
  if sum_dq is None :
    values += [ '' for _value in sum_rhs ]
  else :
    values += [ '%.17g'%value for value in sum_dq ]

  file.write( ','.join(values) + '\n' )

  return


def close_history(file):
  # Close the residual-history file

  if file is None :
    return

  file.close()

  return
