#!/usr/bin/env python3

# Program to verify the maximum eigenvalue used by the LU-SGS operator
#
# lambda=|u.n|+c+2*mu/(rho*d) の式は lusgs_diagonal（対角項）と lusgs_sweep（非対角項）
# の双方が general/thermodynamics.py を共有して使う。式を書き換えると
# LU-SGS の陰解演算子が変わるので、定義そのものをここで固定する。

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

import math

import pytest

from slow.general.thermodynamics import get_max_eigenvalue


SPECIFIC_HEAT_RATIO = 1.40
GAS_CONSTANT        = 8.3144598/28.8e-3


# (密度, 速度, 温度, 粘性係数, セル中心間距離)
STATES = [
  pytest.param(1.0e-4, (600.0, 120.0, 0.0), 300.0, 1.85e-5, 1.0e-3, id='supersonic'),
  pytest.param(1.225,  (30.0, -10.0, 0.0),  288.15, 1.79e-5, 1.0e-2, id='subsonic'),
  pytest.param(1.0e-4, (-450.0, 80.0, 250.0), 250.0, 1.60e-5, 5.0e-4, id='reversed-3d'),
]

NORMALS = [
  pytest.param((1.0, 0.0, 0.0), id='axis-aligned'),
  pytest.param((0.6, 0.8, 0.0), id='oblique-2d'),
  pytest.param((0.4, 0.5, math.sqrt(1.0-0.4**2-0.5**2)), id='oblique-3d'),
]


def call_eigenvalue(dens, velocity, temperature, viscosity, length, area_vec):
  pres = dens*GAS_CONSTANT*temperature
  return get_max_eigenvalue(SPECIFIC_HEAT_RATIO, dens, velocity[0], velocity[1], velocity[2], \
                            pres, viscosity, length, *area_vec)


@pytest.mark.parametrize('area_vec', NORMALS)
@pytest.mark.parametrize('dens, velocity, temperature, viscosity, length', STATES)
def test_eigenvalue_matches_definition(dens, velocity, temperature, viscosity, length, area_vec):
  # lambda=|u.n|+c+2*mu/(rho*d) であること

  pres = dens*GAS_CONSTANT*temperature
  cvel = velocity[0]*area_vec[0] + velocity[1]*area_vec[1] + velocity[2]*area_vec[2]
  sos  = math.sqrt( SPECIFIC_HEAT_RATIO*pres/dens )

  expected = abs(cvel) + sos + 2.0*viscosity/(dens*length)
  actual   = call_eigenvalue(dens, velocity, temperature, viscosity, length, area_vec)

  assert actual == pytest.approx(expected, rel=1.0e-14)


@pytest.mark.parametrize('area_vec', NORMALS)
@pytest.mark.parametrize('dens, velocity, temperature, viscosity, length', STATES)
def test_eigenvalue_is_independent_of_the_normal_direction(dens, velocity, temperature, viscosity, length, area_vec):
  # 法線を反転しても同じ値になること。
  # 後退スイープは法線を反転して呼ぶため、これが成り立たないと前進・後退で lambda が食い違う

  area_vec_flipped = tuple(-component for component in area_vec)

  forward  = call_eigenvalue(dens, velocity, temperature, viscosity, length, area_vec)
  backward = call_eigenvalue(dens, velocity, temperature, viscosity, length, area_vec_flipped)

  assert forward == pytest.approx(backward, rel=1.0e-14)


@pytest.mark.parametrize('area_vec', NORMALS)
@pytest.mark.parametrize('dens, velocity, temperature, viscosity, length', STATES)
def test_viscous_contribution_is_positive(dens, velocity, temperature, viscosity, length, area_vec):
  # 粘性の寄与は必ず lambda を増やす（対角優位性を強める向きに働く）こと

  inviscid = call_eigenvalue(dens, velocity, temperature, 0.0, length, area_vec)
  viscous  = call_eigenvalue(dens, velocity, temperature, viscosity, length, area_vec)

  assert viscous > inviscid


@pytest.mark.parametrize('area_vec', NORMALS)
@pytest.mark.parametrize('dens, velocity, temperature, viscosity, length', STATES)
def test_eigenvalue_bounds_the_convective_speed(dens, velocity, temperature, viscosity, length, area_vec):
  # lambda が反変速度の絶対値を必ず上回ること。
  # 非対角項 0.5*(A-lambda*I) が上流化として機能するための条件

  cvel   = velocity[0]*area_vec[0] + velocity[1]*area_vec[1] + velocity[2]*area_vec[2]
  actual = call_eigenvalue(dens, velocity, temperature, viscosity, length, area_vec)

  assert actual > abs(cvel)
