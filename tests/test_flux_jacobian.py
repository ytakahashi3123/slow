#!/usr/bin/env python3

# Program to verify the convective flux Jacobian used in the LU-SGS sweep
#
# 流束ヤコビアンを法線流束の厳密な微分と突き合わせる。メッシュもソルバ実行も不要なので
# 単体で高速に走る。LU-SGS の陰解演算子に手を入れる際はまずこのテストを通すこと。

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

import math

import numpy as np
import pytest

from slow.time_integration.flux_jacobian import set_flux_jacobian


# Gas properties (src/config.yml の空気に合わせる)
SPECIFIC_HEAT_RATIO = 1.40
GAS_CONSTANT        = 8.3144598/28.8e-3

NUM_CONSERV = 5


def normal_flux(var_conserv, area_vec):
  # 保存変数 Q=(rho, rho*u, rho*v, rho*w, E) に対する法線流束 F(Q).n

  dens, momx, momy, momz, energy = var_conserv
  uvel, vvel, wvel = momx/dens, momy/dens, momz/dens

  cvel = uvel*area_vec[0] + vvel*area_vec[1] + wvel*area_vec[2]
  pres = (SPECIFIC_HEAT_RATIO-1.0)*( energy - 0.50*dens*(uvel**2 + vvel**2 + wvel**2) )

  return np.array([ dens*cvel,                       \
                    momx*cvel + area_vec[0]*pres,    \
                    momy*cvel + area_vec[1]*pres,    \
                    momz*cvel + area_vec[2]*pres,    \
                    (energy + pres)*cvel ])


def jacobian_by_complex_step(var_conserv, area_vec):
  # 複素ステップ微分による厳密ヤコビアン。比較の基準とする
  # 法線流束は四則演算のみで構成され正則なので、刻み幅に依らず機械精度の微分が得られる
  # (中心差分では保存変数の成分が rho=1e-4 から E=1e1 まで数桁にわたるため刻み幅の
  #  選定が難しく、打ち切り誤差が符号誤りと同程度の桁に紛れてしまう)

  step     = 1.0e-20
  jacobian = np.zeros((NUM_CONSERV, NUM_CONSERV))
  for n in range(0, NUM_CONSERV):
    var_perturbed     = np.array(var_conserv, dtype=complex)
    var_perturbed[n]  = var_perturbed[n] + step*1.0j
    jacobian[:,n]     = np.imag( normal_flux(var_perturbed, area_vec) )/step

  return jacobian


def set_state(density, velocity, temperature):
  # 与えられた原始変数から、保存変数と全エンタルピーを組み立てる

  specific_heat_volum = GAS_CONSTANT/(SPECIFIC_HEAT_RATIO-1.0)
  velocity = np.asarray(velocity, dtype=float)

  energy = density*specific_heat_volum*temperature + 0.50*density*np.dot(velocity, velocity)
  pres   = density*GAS_CONSTANT*temperature

  var_conserv = np.array([density, *(density*velocity), energy])
  # 全エンタルピー H=Cv*T+q+p/rho=(E+p)/rho。orbital.get_enthalpy の戻り値と同値
  enth        = (energy + pres)/density

  return var_conserv, enth


# 検査する状態と面法線の組み合わせ
STATES = [
  pytest.param((1.0e-4, (600.0, 120.0, 0.0), 300.0), (0.6, 0.8, 0.0),
               id='supersonic-oblique'),
  pytest.param((1.0e-4, (600.0, 120.0, 0.0), 300.0), (-0.6, -0.8, 0.0),
               id='supersonic-oblique-flipped'),
  pytest.param((1.225, (30.0, -10.0, 0.0), 288.15), (0.0, 1.0, 0.0),
               id='subsonic-axis-aligned'),
  pytest.param((1.0e-4, (-450.0, 80.0, 0.0), 250.0), (1.0, 0.0, 0.0),
               id='supersonic-reversed-flow'),
  pytest.param((1.0e-4, (600.0, 120.0, 250.0), 300.0),
               (0.4, 0.5, math.sqrt(1.0-0.4**2-0.5**2)),
               id='supersonic-three-dimensional'),
]


def get_jacobian_from_state(state, area_vec, eigenvalue=0.0):
  # set_flux_jacobian を原始変数の状態から呼び出す

  var_conserv, enth = set_state(*state)
  uvel, vvel, wvel  = var_conserv[1:4]/var_conserv[0]
  cvel              = uvel*area_vec[0] + vvel*area_vec[1] + wvel*area_vec[2]

  jacobian = np.zeros((NUM_CONSERV, NUM_CONSERV))
  set_flux_jacobian(jacobian, SPECIFIC_HEAT_RATIO, eigenvalue, \
                    cvel, uvel, vvel, wvel, enth, *area_vec)

  return jacobian, var_conserv


@pytest.mark.parametrize('state, area_vec', STATES)
def test_jacobian_matches_exact_derivative(state, area_vec):
  # A(n)=d(F.n)/dQ が厳密な微分と一致すること（固有値項を除いた素のヤコビアン）

  jacobian, var_conserv = get_jacobian_from_state(state, area_vec)
  jacobian_exact        = jacobian_by_complex_step(var_conserv, area_vec)

  np.testing.assert_allclose(jacobian, jacobian_exact, rtol=1.0e-10, \
                             atol=1.0e-10*np.abs(jacobian_exact).max())


@pytest.mark.parametrize('state, area_vec', STATES)
def test_eigenvalue_is_subtracted_from_every_diagonal(state, area_vec):
  # A(n)-lambda*I になっていること。対角 5 成分すべてから引かれている必要がある

  eigenvalue = 1234.5

  jacobian_without, _ = get_jacobian_from_state(state, area_vec, eigenvalue=0.0)
  jacobian_with, _    = get_jacobian_from_state(state, area_vec, eigenvalue=eigenvalue)

  np.testing.assert_allclose(jacobian_without - jacobian_with, \
                             eigenvalue*np.identity(NUM_CONSERV), rtol=1.0e-12, atol=1.0e-9)


@pytest.mark.parametrize('state, area_vec', STATES)
def test_jacobian_is_odd_in_the_face_normal(state, area_vec):
  # A(-n)=-A(n) であること。
  # 後退スイープ (lusgs_sweep.py) は法線を反転させることで自セルの外向き法線に対する
  # A^-=0.5*(A-lambda*I) を得ている。この奇関数性が崩れると後退スイープが壊れる

  jacobian, _         = get_jacobian_from_state(state, area_vec)
  area_vec_flipped    = tuple(-component for component in area_vec)
  jacobian_flipped, _ = get_jacobian_from_state(state, area_vec_flipped)

  np.testing.assert_allclose(jacobian_flipped, -jacobian, rtol=1.0e-12, \
                             atol=1.0e-12*np.abs(jacobian).max())
