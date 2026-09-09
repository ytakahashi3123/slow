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


def test_second_order_backward_difference_falls_back_while_the_history_is_short():
  """
  過去の解が 1 段しかないステップでは BDF2 を使わず 1 次後退差分で立ち上げること。

  Q^(n-1) が無いまま BDF2 の式を使うと Q^(n-1)=Q^(n) となり、
  (1.5Q-2Q^n+0.5Q^n)/dt = 1.5*(Q-Q^n)/dt すなわち刻み幅 dt/1.5 の
  1 次後退差分になる。この O(dt) の誤差は 1 ステップだけでも最後まで残り、
  全体の時間精度を 2 次から 1 次に落とす（実測: 25 セルのケースで
  観測次数 2.09 -> 1.21、最細の dt で誤差 9.8 倍）。
  """

  volume, var_dt, var_conserv, var_conserv_prev = make_state()
  timestep_outer = 1.0e-4
  config = make_config(kind_backward_difference='2nd_backward_diff', timestep_outer=timestep_outer)

  for n_cell in range(0, NUM_CELL):
    diag_time, dq_unst, face_scale = lusgs.get_time_term(config, n_cell, volume, var_dt, \
                                                         var_conserv, var_conserv_prev, \
                                                         num_conserv_prev_level=1)
    # 1 次後退差分と完全に一致すること
    diag_ref, dq_ref, scale_ref = lusgs.get_time_term(
        make_config(kind_backward_difference='1st_backward_diff', timestep_outer=timestep_outer),
        n_cell, volume, var_dt, var_conserv, var_conserv_prev)

    assert diag_time == pytest.approx(diag_ref, rel=1.0e-14)
    np.testing.assert_allclose(dq_unst, dq_ref, rtol=1.0e-14)
    assert face_scale == scale_ref


def test_second_order_backward_difference_is_used_once_the_history_is_complete():
  # 2 段揃えば BDF2 に戻ること。既定値も 2 段揃った状態とする

  volume, var_dt, var_conserv, var_conserv_prev = make_state()
  config = make_config(kind_backward_difference='2nd_backward_diff')

  for n_cell in range(0, NUM_CELL):
    with_level = lusgs.get_time_term(config, n_cell, volume, var_dt, var_conserv, \
                                     var_conserv_prev, num_conserv_prev_level=2)
    default    = lusgs.get_time_term(config, n_cell, volume, var_dt, var_conserv, var_conserv_prev)

    assert with_level[0] == pytest.approx(default[0], rel=1.0e-14)
    np.testing.assert_allclose(with_level[1], default[1], rtol=1.0e-14)
    assert with_level[2] == default[2]
    # BDF2 の対角は 1 次後退差分より大きい（係数 1.5）
    first_order = lusgs.get_time_term(
        make_config(kind_backward_difference='1st_backward_diff'),
        n_cell, volume, var_dt, var_conserv, var_conserv_prev)
    assert with_level[0] > first_order[0]


def test_first_order_backward_difference_is_unaffected_by_the_history_length():
  # 1 次後退差分は Q^(n) だけで足りるので段数に依らない

  volume, var_dt, var_conserv, var_conserv_prev = make_state()
  config = make_config(kind_backward_difference='1st_backward_diff')

  for n_cell in range(0, NUM_CELL):
    short = lusgs.get_time_term(config, n_cell, volume, var_dt, var_conserv, \
                                var_conserv_prev, num_conserv_prev_level=1)
    full  = lusgs.get_time_term(config, n_cell, volume, var_dt, var_conserv, \
                                var_conserv_prev, num_conserv_prev_level=2)

    assert short[0] == full[0]
    np.testing.assert_allclose(short[1], full[1], rtol=1.0e-14)
    assert short[2] == full[2]


def test_previous_conservative_shift_counts_up_the_available_history():
  # set_conservative_previous が段数を数え上げ、BDF2 に必要な 2 段で飽和すること

  from slow.time_integration.time_integration import time_integration

  volume, var_dt, var_conserv, var_conserv_prev = make_state()
  config = make_config()
  ti     = time_integration()

  var_conserv_prev, level = ti.set_conservative_previous(config, var_conserv, var_conserv_prev, 1)
  assert level == 2
  # Q^(n) には直前の解が入る
  np.testing.assert_allclose(var_conserv_prev[0,:,:], var_conserv, rtol=1.0e-14)

  var_conserv_prev, level = ti.set_conservative_previous(config, var_conserv, var_conserv_prev, level)
  assert level == 2


@pytest.mark.parametrize('kind_steady_mode', ['steady', 'unsteady'])
@pytest.mark.parametrize('var_dt_bad', [0.0, -1.0e-7])
def test_non_positive_pseudo_timestep_stops_the_program(kind_steady_mode, var_dt_bad):
  """
  volume/var_dt の 0 除算で黙って進まないこと。

  ガードが無いと対角が inf になり、delta Q が 0 になる。例外も警告も出ないため
  「反復しても解が動かない」計算になってしまう。設定間違い（courant_number: 0 など）
  で起こりうるので、原因を示して止める。
  """

  volume, var_dt, var_conserv, var_conserv_prev = make_state()
  var_dt = var_dt.copy()
  var_dt[1] = var_dt_bad
  config = make_config(kind_steady_mode=kind_steady_mode)

  # 正常なセルはそのまま計算できる
  lusgs.get_time_term(config, 0, volume, var_dt, var_conserv, var_conserv_prev)

  with pytest.raises(SystemExit):
    lusgs.get_time_term(config, 1, volume, var_dt, var_conserv, var_conserv_prev)


@pytest.mark.parametrize('timestep_outer_bad', [0.0, -1.0e-6])
def test_non_positive_timestep_outer_stops_the_program(timestep_outer_bad):
  # 非定常では volume/timestep_outer でも割るので、こちらも確かめる

  volume, var_dt, var_conserv, var_conserv_prev = make_state()
  config = make_config(kind_steady_mode='unsteady', timestep_outer=timestep_outer_bad)

  with pytest.raises(SystemExit):
    lusgs.get_time_term(config, 0, volume, var_dt, var_conserv, var_conserv_prev)


def test_positive_timesteps_are_left_alone():
  # ガードが正常な計算を邪魔しないこと（対角は volume/var_dt のまま）

  volume, var_dt, var_conserv, var_conserv_prev = make_state()
  config = make_config(kind_steady_mode='steady')

  for n_cell in range(0, NUM_CELL):
    diag_time, _, _ = lusgs.get_time_term(config, n_cell, volume, var_dt, \
                                          var_conserv, var_conserv_prev)
    assert diag_time == pytest.approx(volume[n_cell]/var_dt[n_cell], rel=1.0e-14)


# ---------------------------------------------------------------- 時間項の 2 つの実装が一致すること

TIME_SETTING_CASES = [
  pytest.param('steady',   '2nd_backward_diff', 1.01, id='steady'),
  pytest.param('unsteady', '1st_backward_diff', 1.01, id='unsteady-bdf1'),
  pytest.param('unsteady', '2nd_backward_diff', 1.01, id='unsteady-bdf2'),
]


@pytest.mark.parametrize('kind_steady_mode, kind_backward_difference, lusgs_beta', TIME_SETTING_CASES)
@pytest.mark.parametrize('num_conserv_prev_level', [1, 2])
def test_time_term_kernel_matches_the_reference(kind_steady_mode, kind_backward_difference,
                                                lusgs_beta, num_conserv_prev_level):
  """
  セルループのカーネルが get_time_term と同じ値を出すこと。

  get_time_term はセルごとに config の辞書を 4 回引くので、7,000 セル規模では
  それが対角の計算時間の大半を占めていた。カーネル側は設定を整数に解決してから
  受け取るが、式は同じでなければならない。ここが崩れると陰解演算子が静かに狂う。
  """

  from slow.time_integration import lusgs_diagonal_kernel

  volume, var_dt, var_conserv, var_conserv_prev = make_state()
  timestep_outer = 1.0e-4
  config = make_config(kind_steady_mode=kind_steady_mode,
                       kind_backward_difference=kind_backward_difference,
                       lusgs_beta=lusgs_beta, timestep_outer=timestep_outer)

  var_rhs = np.arange(1.0, 1.0 + NUM_CONSERV*NUM_CELL).reshape(NUM_CONSERV, NUM_CELL)*1.0e-3

  # 面の寄与に相当する値を入れておく（対角がゼロから始まると face_scale が効かない）
  diagonal_face = np.array([3.0e-2, 7.0e-2, 1.1e-1])

  # --参照: get_time_term をセルごとに呼ぶ従来の形
  diag_ref = diagonal_face.copy()
  dq_ref   = np.zeros((NUM_CONSERV, NUM_CELL))
  for n_cell in range(0, NUM_CELL):
    diag_time, dq_unst, face_scale = lusgs.get_time_term(
        config, n_cell, volume, var_dt, var_conserv, var_conserv_prev, num_conserv_prev_level)
    diag_ref[n_cell] = diag_time + face_scale*diag_ref[n_cell]
    dq_ref[:,n_cell] = ( -var_rhs[:,n_cell]-dq_unst )/diag_ref[n_cell]

  # --カーネル
  kind_time, timestep_outer_resolved, face_scale = lusgs.get_time_setting(config, num_conserv_prev_level)
  diag_new = diagonal_face.copy()
  dq_new   = np.zeros((NUM_CONSERV, NUM_CELL))
  lusgs_diagonal_kernel.apply_time_term_scalar(
      kind_time, lusgs.KIND_TIME_STEADY, lusgs.KIND_TIME_BDF2,
      NUM_CELL, NUM_CONSERV, timestep_outer_resolved, face_scale,
      volume, var_dt, var_conserv, var_conserv_prev, var_rhs, diag_new, dq_new)

  np.testing.assert_allclose(diag_new, diag_ref, rtol=0.0, atol=0.0)
  np.testing.assert_allclose(dq_new, dq_ref, rtol=0.0, atol=0.0)


def test_time_setting_resolves_the_startup_fallback():
  # 段数が足りないときに BDF2 が BDF1 に落ちることを、解決した識別子の側でも固定する

  config = make_config(kind_backward_difference='2nd_backward_diff')

  assert lusgs.get_time_setting(config, 2)[0] == lusgs.KIND_TIME_BDF2
  assert lusgs.get_time_setting(config, 1)[0] == lusgs.KIND_TIME_BDF1


@pytest.mark.parametrize('kind_steady_mode, kind_backward_difference',
                         [('transient', '2nd_backward_diff'),
                          ('unsteady', '3rd_backward_diff')])
def test_unknown_time_scheme_stops_get_time_setting(kind_steady_mode, kind_backward_difference):
  config = make_config(kind_steady_mode=kind_steady_mode,
                       kind_backward_difference=kind_backward_difference)

  with pytest.raises(SystemExit):
    lusgs.get_time_setting(config)
