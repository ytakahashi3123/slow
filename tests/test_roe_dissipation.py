#!/usr/bin/env python3

# Program to verify the matrix dissipation |A_Roe| = P|Lambda|P^-1
#
# 固有ベクトル行列 P は要素を直接書き下しているので、
# 取り違えが入っていないことを P*Lambda*P^-1 == A で確かめる。
# ここでの A は tests/test_flux_jacobian.py で厳密性を確認済みの
# time_integration/flux_jacobian.py が作る解析ヤコビアンである。

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

import math

import numpy as np
import pytest

from slow.time_integration.flux_jacobian import set_flux_jacobian
from slow.time_integration.roe_dissipation import (get_eigenvalues, get_roe_average,
                                                   set_absolute_jacobian, set_pmatrix)


SPECIFIC_HEAT_RATIO = 1.40
GAS_CONSTANT        = 8.3144598/28.8e-3

NUM_CONSERV = 5


# (密度, 速度, 温度)
STATES = [
  pytest.param(1.0e-4, (600.0, 120.0, 0.0),   300.0,  id='supersonic'),
  pytest.param(1.225,  (30.0, -10.0, 0.0),    288.15, id='subsonic'),
  pytest.param(1.0e-4, (-450.0, 80.0, 0.0),   250.0,  id='reversed'),
  pytest.param(1.0e-4, (600.0, 120.0, 250.0), 300.0,  id='three-dimensional'),
]

# SLOW の格子は 2 次元なので法線の z 成分は常に 0。3 次元法線も念のため入れる
NORMALS = [
  pytest.param((1.0, 0.0, 0.0), id='axis-aligned'),
  pytest.param((0.6, 0.8, 0.0), id='oblique-2d'),
  pytest.param((-0.6, -0.8, 0.0), id='oblique-2d-flipped'),
  pytest.param((0.4, 0.5, math.sqrt(1.0-0.4**2-0.5**2)), id='oblique-3d'),
]


def set_state(density, velocity, temperature):
  # 原始変数から (密度, 速度, 音速, 全エンタルピー) を組み立てる

  specific_heat_volum = GAS_CONSTANT/(SPECIFIC_HEAT_RATIO-1.0)
  velocity = np.asarray(velocity, dtype=float)

  pres = density*GAS_CONSTANT*temperature
  sos  = math.sqrt( SPECIFIC_HEAT_RATIO*pres/density )
  enth = specific_heat_volum*temperature + 0.50*np.dot(velocity, velocity) + pres/density

  return density, velocity, sos, enth


def get_analytic_jacobian(velocity, enth, area_vec):
  # flux_jacobian.py の A(n)（固有値項なし）

  jacobian = np.zeros((NUM_CONSERV, NUM_CONSERV))
  cvel     = float(np.dot(velocity, area_vec))
  set_flux_jacobian(jacobian, SPECIFIC_HEAT_RATIO, 0.0, cvel, \
                    velocity[0], velocity[1], velocity[2], enth, *area_vec)

  return jacobian


@pytest.mark.parametrize('area_vec', NORMALS)
@pytest.mark.parametrize('density, velocity, temperature', STATES)
def test_pmatrix_diagonalizes_the_flux_jacobian(density, velocity, temperature, area_vec):
  # P*diag(Lambda)*P^-1 == A(n) であること。
  # P の移植ミスと固有値の並び順の食い違いを同時に検出する

  dens, vel, sos, enth = set_state(density, velocity, temperature)

  pmatrix = np.zeros((NUM_CONSERV, NUM_CONSERV))
  set_pmatrix(pmatrix, SPECIFIC_HEAT_RATIO, dens, vel, sos, *area_vec)
  eigenvalues = get_eigenvalues(vel, sos, *area_vec)

  reconstructed = pmatrix @ ( eigenvalues[:,None] * np.linalg.inv(pmatrix) )
  expected      = get_analytic_jacobian(vel, enth, area_vec)

  np.testing.assert_allclose(reconstructed, expected, rtol=1.0e-9, \
                             atol=1.0e-9*np.abs(expected).max())


@pytest.mark.parametrize('area_vec', NORMALS)
@pytest.mark.parametrize('density, velocity, temperature', STATES)
def test_absolute_jacobian_has_the_absolute_eigenvalues(density, velocity, temperature, area_vec):
  # |A| の固有値が |Lambda| になっていること

  dens, vel, sos, _ = set_state(density, velocity, temperature)

  absjacobian = np.zeros((NUM_CONSERV, NUM_CONSERV))
  set_absolute_jacobian(absjacobian, SPECIFIC_HEAT_RATIO, dens, vel, sos, *area_vec)

  expected = np.sort( np.abs( get_eigenvalues(vel, sos, *area_vec) ) )
  actual   = np.sort( np.abs( np.linalg.eigvals(absjacobian) ) )

  np.testing.assert_allclose(actual, expected, rtol=1.0e-8, atol=1.0e-8*expected.max())


@pytest.mark.parametrize('area_vec', NORMALS)
@pytest.mark.parametrize('density, velocity, temperature', STATES)
def test_absolute_jacobian_is_independent_of_the_normal_direction(density, velocity, temperature, area_vec):
  # |A(-n)| == |A(n)| であること。
  # 前進・後退スイープが同じ散逸行列を共有できる根拠になる

  dens, vel, sos, _ = set_state(density, velocity, temperature)
  area_vec_flipped  = tuple(-component for component in area_vec)

  forward  = np.zeros((NUM_CONSERV, NUM_CONSERV))
  backward = np.zeros((NUM_CONSERV, NUM_CONSERV))
  set_absolute_jacobian(forward,  SPECIFIC_HEAT_RATIO, dens, vel, sos, *area_vec)
  set_absolute_jacobian(backward, SPECIFIC_HEAT_RATIO, dens, vel, sos, *area_vec_flipped)

  np.testing.assert_allclose(backward, forward, rtol=1.0e-8, atol=1.0e-8*np.abs(forward).max())


def test_absolute_jacobian_equals_the_jacobian_for_supersonic_normal_flow():
  # 法線方向が超音速なら固有値が全て同符号になるので |A| == A になること。
  # P|Lambda|P^-1 が正しく組めていれば自動的に成り立つ強い条件

  dens, vel, sos, enth = set_state(1.0e-4, (900.0, 0.0, 0.0), 300.0)
  area_vec = (1.0, 0.0, 0.0)
  assert float(np.dot(vel, area_vec)) > sos, 'テスト条件が法線方向超音速になっていない'

  absjacobian = np.zeros((NUM_CONSERV, NUM_CONSERV))
  set_absolute_jacobian(absjacobian, SPECIFIC_HEAT_RATIO, dens, vel, sos, *area_vec)
  expected = get_analytic_jacobian(vel, enth, area_vec)

  np.testing.assert_allclose(absjacobian, expected, rtol=1.0e-8, \
                             atol=1.0e-8*np.abs(expected).max())


@pytest.mark.parametrize('density, velocity, temperature', STATES)
def test_roe_average_reduces_to_the_common_state(density, velocity, temperature):
  # 左右が同じ状態なら Roe 平均もその状態に一致すること

  dens, vel, sos, enth = set_state(density, velocity, temperature)

  dens_roe, vel_roe, sos_roe = get_roe_average(SPECIFIC_HEAT_RATIO, dens, vel, enth, dens, vel, enth)

  assert dens_roe == pytest.approx(dens, rel=1.0e-13)
  np.testing.assert_allclose(vel_roe, vel, rtol=1.0e-13, atol=1.0e-13*np.abs(vel).max())
  assert sos_roe == pytest.approx(sos, rel=1.0e-12)


def test_roe_average_lies_between_the_two_states():
  # Roe 平均が左右の状態の間に入ること

  dens_a, vel_a, _, enth_a = set_state(1.0e-4, (600.0, 120.0, 0.0), 300.0)
  dens_b, vel_b, _, enth_b = set_state(4.0e-4, (200.0, -50.0, 0.0), 900.0)

  dens_roe, vel_roe, sos_roe = get_roe_average(SPECIFIC_HEAT_RATIO, dens_a, vel_a, enth_a, \
                                               dens_b, vel_b, enth_b)

  assert dens_a < dens_roe < dens_b
  for n in range(0, 3):
    assert min(vel_a[n], vel_b[n]) <= vel_roe[n] <= max(vel_a[n], vel_b[n])
  assert sos_roe > 0.0


def test_upwind_offdiagonal_vanishes_for_supersonic_normal_flow():
  # 法線方向が超音速なら非対角ブロック A^- = 0.5*(A-|A|) が厳密に 0 になること。
  # 上流側から下流側への影響しか残らない（上流化が効いている）ことの確認であり、
  # lusgs_sweep の行列散逸版が正しく組めているかの要になる

  dens, vel, sos, enth = set_state(1.0e-4, (900.0, 0.0, 0.0), 300.0)
  area_vec = (1.0, 0.0, 0.0)

  absjacobian = np.zeros((NUM_CONSERV, NUM_CONSERV))
  set_absolute_jacobian(absjacobian, SPECIFIC_HEAT_RATIO, dens, vel, sos, *area_vec)
  jacobian = get_analytic_jacobian(vel, enth, area_vec)

  offdiagonal = 0.50*( jacobian - absjacobian )

  np.testing.assert_allclose(offdiagonal, np.zeros((NUM_CONSERV, NUM_CONSERV)), \
                             rtol=0.0, atol=1.0e-8*np.abs(jacobian).max())
