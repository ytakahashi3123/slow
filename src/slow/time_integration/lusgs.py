#!/usr/bin/env python3

# Program to provide settings shared by the LU-SGS routines

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

# Kind of dissipation used for the LU-SGS implicit operator
# --'scalar': 0.5*(A-lambda*I)*S。lambda=max(lambda_a,lambda_b) のスカラー散逸
# --'matrix': 0.5*(A-|A_Roe|)*S。Roe 平均による行列散逸。
#             対角もブロック行列になる
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
    print('Error in kind_lusgs_dissipation of control file: ', kind_dissipation)
    print('Available: ', ', '.join(AVAIL_KIND_DISSIPATION))
    print('Program stopped')
    exit()

  return kind_dissipation


def get_time_term(config, n_cell, volume, var_dt, var_conserv, var_conserv_prev):
  """
  Time-derivative contribution to the LU-SGS system for one cell

  戻り値:
    --diag_time: 対角に加える時間項（スカラー）
    --dq_unst:   右辺に加える非定常項（保存変数と同じ長さのベクトル。定常計算では 0）
    --face_scale: 面の寄与にかける係数（定常計算では過剰緩和係数 lusgs_beta が入る）

  """

  kind_steady_mode         = config['time_integration']['kind_steady_mode']
  lusgs_beta               = config['time_integration']['lusgs_beta']
  kind_backward_difference = config['time_integration']['kind_backward_difference']
  var_dt_const             = config['time_integration']['timestep_outer']

  if kind_steady_mode == 'steady' :
    # Steady flow
    return volume[n_cell]/var_dt[n_cell], 0.0, 0.50*lusgs_beta

  elif kind_steady_mode == 'unsteady' :
    # Unsteady flow with Pseudo-time stepping
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
      print('Error in kind_backward_difference of control file: ', kind_backward_difference )
      print('Program stopped')
      exit()

  else :
    print('Error in kind_steady_mode of control file: ', kind_steady_mode )
    print('Program stopped')
    exit()
