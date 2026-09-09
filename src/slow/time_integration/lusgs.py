#!/usr/bin/env python3

# Program to provide settings shared by the LU-SGS routines

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

import logging
import sys

logger = logging.getLogger(__name__)

# Kind of dissipation used for the LU-SGS implicit operator
# --'scalar': 0.5*(A-lambda*I)*S。lambda=max(lambda_a,lambda_b) のスカラー散逸
# --'matrix': 0.5*(A-|A_Roe|)*S。Roe 平均による行列散逸。
#             対角もブロック行列になる
# 2 次後退差分 (BDF2) が必要とする過去の解の段数: Q^(n) と Q^(n-1)
NUM_PREV_LEVEL_REQUIRED_BDF2 = 2

KIND_DISSIPATION_SCALAR = 'scalar'
KIND_DISSIPATION_MATRIX = 'matrix'
AVAIL_KIND_DISSIPATION  = (KIND_DISSIPATION_SCALAR, KIND_DISSIPATION_MATRIX)


def get_kind_dissipation(config):
  """
  Kind of dissipation for the LU-SGS implicit operator, given by the control file

  設定が無い場合は従来の挙動である 'scalar' とする（既存の config.yml をそのまま使えるように）。
  """

  kind_dissipation = str( config['time_integration'].get('kind_lusgs_dissipation',
                                                          KIND_DISSIPATION_SCALAR) )

  if kind_dissipation not in AVAIL_KIND_DISSIPATION:
    logger.info('Error in kind_lusgs_dissipation of control file:  %s', kind_dissipation)
    logger.info('Available:  %s', ', '.join(AVAIL_KIND_DISSIPATION))
    logger.info('Program stopped')
    exit()

  return kind_dissipation


def check_timestep_positive(name, value, n_cell):
  """
  Stop with a readable message when a time step used as a divisor is not positive

  volume/var_dt をそのまま計算すると 0 除算で対角が inf になり、delta Q が 0 になる。
  例外も警告も出ないまま「反復しても解が動かない」計算になってしまうので、
  ここで止めて原因を示す。
  """

  if value > 0.0 :
    return

  logger.error('Error: %s must be positive but is %s (cell %s)', name, value, n_cell)
  logger.error('--check courant_number / timestep_constant / timestep_outer of the control file')
  logger.error('Program stopped')
  # 呼び出し側のスクリプトから失敗を検知できるよう 0 以外で終了する
  sys.exit(1)


def get_time_term(config, n_cell, volume, var_dt, var_conserv, var_conserv_prev,
                  num_conserv_prev_level=NUM_PREV_LEVEL_REQUIRED_BDF2):
  """
  Time-derivative contribution to the LU-SGS system for one cell

  戻り値:
    --diag_time: 対角に加える時間項（スカラー）
    --dq_unst:   右辺に加える非定常項（保存変数と同じ長さのベクトル。定常計算では 0）
    --face_scale: 面の寄与にかける係数（定常計算では過剰緩和係数 lusgs_beta が入る）

  --num_conserv_prev_level: var_conserv_prev に入っている「本物の」過去の解の段数。
      2 段目 Q^(n-1) が揃っていない最初のステップで 2 次後退差分を使うと、
      Q^(n-1)=Q^(n) となって式が 1.5*(Q^(n+1)-Q^(n))/dt に退化し、
      刻み幅 dt/1.5 の 1 次後退差分を解くことになる。この O(dt) の誤差は
      1 ステップだけでも最後まで残り、全体の時間精度を 1 次に落とす。
      そのため段数が足りないステップは 1 次後退差分で立ち上げる。
  """

  kind_steady_mode         = config['time_integration']['kind_steady_mode']
  lusgs_beta               = config['time_integration']['lusgs_beta']
  kind_backward_difference = config['time_integration']['kind_backward_difference']
  var_dt_const             = config['time_integration']['timestep_outer']

  # 過去の解が足りないステップは 1 次後退差分で立ち上げる
  if kind_backward_difference == '2nd_backward_diff' and \
     num_conserv_prev_level < NUM_PREV_LEVEL_REQUIRED_BDF2 :
    kind_backward_difference = '1st_backward_diff'

  # 0 除算を黙って inf にしないよう、割る前に確かめる
  check_timestep_positive('pseudo time step (var_dt)', var_dt[n_cell], n_cell)

  if kind_steady_mode == 'steady' :
    # Steady flow
    return volume[n_cell]/var_dt[n_cell], 0.0, 0.50*lusgs_beta

  elif kind_steady_mode == 'unsteady' :
    # Unsteady flow with Pseudo-time stepping
    check_timestep_positive('timestep_outer', var_dt_const, n_cell)

    if kind_backward_difference == '2nd_backward_diff' :
      # - 2nd order accuracy backward difference
      # - (Volume*(3/2dt+1/d_tau) + 0.5*eigenvalue)
      diag_time = 1.50*volume[n_cell]/var_dt_const + volume[n_cell]/var_dt[n_cell]
      dq_unst   = ( 1.50*var_conserv[:,n_cell] - 2.0*var_conserv_prev[0,:,n_cell] +0.50*var_conserv_prev[1,:,n_cell] )*volume[n_cell]/var_dt_const
      return diag_time, dq_unst, 0.50

    elif kind_backward_difference == '1st_backward_diff' :
      # - 1st order accuracy backward difference
      # - (Volume*(1/dt+1/d_tau) + 0.5*eigenvalue)
      diag_time = volume[n_cell]/var_dt_const + volume[n_cell]/var_dt[n_cell]
      dq_unst   = ( var_conserv[:,n_cell] - var_conserv_prev[0,:,n_cell] )*volume[n_cell]/var_dt_const
      return diag_time, dq_unst, 0.50

    else :
      logger.info('Error in kind_backward_difference of control file:  %s', kind_backward_difference)
      logger.info('Program stopped')
      exit()

  else :
    logger.info('Error in kind_steady_mode of control file:  %s', kind_steady_mode)
    logger.info('Program stopped')
    exit()
