#!/usr/bin/env python3

# Program to verify the LU-SGS settings shared by lusgs_diagonal and lusgs_sweep
#
# 散逸の種類の選択と、時間項の組み立て（定常／後退差分）を固定する。
# get_time_term は lusgs_diagonal の scalar 版と matrix 版が共有しているので、
# ここが崩れると両方の陰解演算子が同時に狂う。

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

import numpy as np
import pytest

from slow.time_integration import lusgs


NUM_CONSERV = 5
NUM_CELL    = 3


def make_config(kind_steady_mode='unsteady', kind_backward_difference='2nd_backward_diff',
                lusgs_beta=1.01, timestep_outer=1.0e-4, kind_lusgs_dissipation=None):
  config = { 'time_integration': { 'kind_steady_mode': kind_steady_mode,           \
                                   'kind_backward_difference': kind_backward_difference, \
                                   'lusgs_beta': lusgs_beta,                       \
                                   'timestep_outer': timestep_outer } }
  if kind_lusgs_dissipation is not None:
    config['time_integration']['kind_lusgs_dissipation'] = kind_lusgs_dissipation

  return config


def make_state():
  volume          = np.array([2.0, 3.0, 4.0])
  var_dt          = np.array([1.0e-6, 2.0e-6, 4.0e-6])
  var_conserv     = np.arange(1.0, 1.0 + NUM_CONSERV*NUM_CELL).reshape(NUM_CONSERV, NUM_CELL)
  var_conserv_prev = np.zeros((2, NUM_CONSERV, NUM_CELL))
  var_conserv_prev[0,:,:] = var_conserv - 1.0
  var_conserv_prev[1,:,:] = var_conserv - 3.0

  return volume, var_dt, var_conserv, var_conserv_prev


def test_dissipation_defaults_to_scalar_when_the_key_is_absent():
  # 既存の config.yml をそのまま使えること（従来の挙動が既定）

  assert lusgs.get_kind_dissipation(make_config()) == lusgs.KIND_DISSIPATION_SCALAR


@pytest.mark.parametrize('kind', [lusgs.KIND_DISSIPATION_SCALAR, lusgs.KIND_DISSIPATION_MATRIX])
def test_dissipation_is_taken_from_the_control_file(kind):
  assert lusgs.get_kind_dissipation(make_config(kind_lusgs_dissipation=kind)) == kind


def test_unknown_dissipation_stops_the_program():
  # 綴り間違いを黙って既定値で流さないこと

  with pytest.raises(SystemExit):
    lusgs.get_kind_dissipation(make_config(kind_lusgs_dissipation='roe'))


def test_time_term_for_steady_flow():
  # 定常: 対角は V/d_tau、非定常項は無し、面の寄与に過剰緩和係数 lusgs_beta が掛かる

  volume, var_dt, var_conserv, var_conserv_prev = make_state()
  config = make_config(kind_steady_mode='steady', lusgs_beta=1.01)

  for n_cell in range(0, NUM_CELL):
    diag_time, dq_unst, face_scale = lusgs.get_time_term(config, n_cell, volume, var_dt, \
                                                         var_conserv, var_conserv_prev)
    assert diag_time == pytest.approx(volume[n_cell]/var_dt[n_cell], rel=1.0e-14)
    assert dq_unst == 0.0
    assert face_scale == pytest.approx(0.50*1.01, rel=1.0e-14)


def test_time_term_for_first_order_backward_difference():
  # 1 次後退差分: 対角に V/dt を追加、非定常項は (Q-Q^n)*V/dt

  volume, var_dt, var_conserv, var_conserv_prev = make_state()
  timestep_outer = 1.0e-4
  config = make_config(kind_backward_difference='1st_backward_diff', timestep_outer=timestep_outer)

  for n_cell in range(0, NUM_CELL):
    diag_time, dq_unst, face_scale = lusgs.get_time_term(config, n_cell, volume, var_dt, \
                                                         var_conserv, var_conserv_prev)
    assert diag_time == pytest.approx(volume[n_cell]/timestep_outer + volume[n_cell]/var_dt[n_cell], rel=1.0e-14)
    expected = ( var_conserv[:,n_cell] - var_conserv_prev[0,:,n_cell] )*volume[n_cell]/timestep_outer
    np.testing.assert_allclose(dq_unst, expected, rtol=1.0e-14)
    # 非定常では lusgs_beta を掛けない（意図的な仕様）
    assert face_scale == 0.50


def test_time_term_for_second_order_backward_difference():
  # 2 次後退差分 (BDF2): 対角に 1.5*V/dt (=3V/(2dt))、
  # 非定常項は (1.5Q-2Q^n+0.5Q^n-1)*V/dt

  volume, var_dt, var_conserv, var_conserv_prev = make_state()
  timestep_outer = 1.0e-4
  config = make_config(kind_backward_difference='2nd_backward_diff', timestep_outer=timestep_outer)

  for n_cell in range(0, NUM_CELL):
    diag_time, dq_unst, face_scale = lusgs.get_time_term(config, n_cell, volume, var_dt, \
                                                         var_conserv, var_conserv_prev)
    assert diag_time == pytest.approx(1.50*volume[n_cell]/timestep_outer + volume[n_cell]/var_dt[n_cell], rel=1.0e-14)
    expected = ( 1.50*var_conserv[:,n_cell] - 2.0*var_conserv_prev[0,:,n_cell] \
                 + 0.50*var_conserv_prev[1,:,n_cell] )*volume[n_cell]/timestep_outer
    np.testing.assert_allclose(dq_unst, expected, rtol=1.0e-14)
    assert face_scale == 0.50


@pytest.mark.parametrize('kind_steady_mode, kind_backward_difference',
                         [('transient', '2nd_backward_diff'),
                          ('unsteady', '3rd_backward_diff')])
def test_unknown_time_scheme_stops_the_program(kind_steady_mode, kind_backward_difference):
  volume, var_dt, var_conserv, var_conserv_prev = make_state()
  config = make_config(kind_steady_mode=kind_steady_mode, \
                       kind_backward_difference=kind_backward_difference)

  with pytest.raises(SystemExit):
    lusgs.get_time_term(config, 0, volume, var_dt, var_conserv, var_conserv_prev)

