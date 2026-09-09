#!/usr/bin/env python3

# Program to verify the kernel thermodynamics against the orbital methods
#
# 面ループの内側から呼ぶために general/thermodynamics.py に同じ式を置いている。
# 二つが 1 ULP でも食い違うと、numba を掛けた版と掛けない版で結果が変わる。
# ここで同一であることを固定する。

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/10

import numpy as np
import pytest

from slow.general import thermodynamics
from slow.orbital.orbital import orbital


# 検査する状態の範囲。密度と圧力は数桁にわたらせる
NUM_SAMPLE = 20000


def make_states(seed=0):
  rng = np.random.default_rng(seed)

  dens = 10.0**rng.uniform(-5.0, 1.0, NUM_SAMPLE)
  pres = 10.0**rng.uniform(-2.0, 6.0, NUM_SAMPLE)
  temp = rng.uniform(100.0, 4000.0, NUM_SAMPLE)
  vel  = rng.normal(0.0, 500.0, (3, NUM_SAMPLE))
  visc = 10.0**rng.uniform(-6.0, -4.0, NUM_SAMPLE)
  leng = 10.0**rng.uniform(-4.0, -1.0, NUM_SAMPLE)
  vec  = rng.normal(0.0, 1.0, (3, NUM_SAMPLE))
  vec  = vec/np.linalg.norm(vec, axis=0)

  return dens, pres, temp, vel, visc, leng, vec


SPECIFIC_HEAT_RATIO = 1.40
SPECIFIC_HEAT_VOLUM = 8.3144598/28.8e-3/(SPECIFIC_HEAT_RATIO-1.0)


def test_speedofsound_matches_the_orbital_method():
  dens, pres, _, _, _, _, _ = make_states()

  for n in range(0, NUM_SAMPLE):
    expected = orbital.get_speedofsound("self", SPECIFIC_HEAT_RATIO, dens[n], pres[n])
    actual   = thermodynamics.get_speedofsound(SPECIFIC_HEAT_RATIO, dens[n], pres[n])
    assert actual == expected, f'sample {n}: {actual!r} != {expected!r}'


def test_enthalpy_matches_the_orbital_method():
  dens, pres, temp, vel, _, _, _ = make_states()

  for n in range(0, NUM_SAMPLE):
    expected = orbital.get_enthalpy("self", SPECIFIC_HEAT_VOLUM, dens[n], temp[n],
                                    [vel[0,n], vel[1,n], vel[2,n]], pres[n])
    actual   = thermodynamics.get_enthalpy(SPECIFIC_HEAT_VOLUM, dens[n], temp[n],
                                           vel[0,n], vel[1,n], vel[2,n], pres[n])
    assert actual == expected, f'sample {n}: {actual!r} != {expected!r}'


@pytest.mark.parametrize('vec_sign', [1.0, -1.0])
def test_max_eigenvalue_does_not_depend_on_the_normal_direction(vec_sign):
  # lambda は |u.n| を使うので法線の向きに依らない（後退スイープがこれに依存する）

  dens, pres, _, vel, visc, leng, vec = make_states(seed=1)

  for n in range(0, 200):
    base = thermodynamics.get_max_eigenvalue(SPECIFIC_HEAT_RATIO, dens[n],
                                             vel[0,n], vel[1,n], vel[2,n], pres[n],
                                             visc[n], leng[n], vec[0,n], vec[1,n], vec[2,n])
    flip = thermodynamics.get_max_eigenvalue(SPECIFIC_HEAT_RATIO, dens[n],
                                             vel[0,n], vel[1,n], vel[2,n], pres[n],
                                             visc[n], leng[n],
                                             vec_sign*vec[0,n], vec_sign*vec[1,n], vec_sign*vec[2,n])
    assert flip == base
