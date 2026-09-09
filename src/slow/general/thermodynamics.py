#!/usr/bin/env python3

# Program to provide the thermodynamic relations used inside the numerical kernels

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/10

import numpy as np

from slow.general.jit import kernel


# orbital 側の同名メソッドと同じ式を、カーネルから呼べる素の関数として置く。
# orbital のものはスカラー 1 個ずつを扱う汎用の入口で、こちらは面ループの内側から
# 呼ばれる（numba はクラスのメソッドもリストも受け取れないため、速度の配列ではなく
# 成分を 3 つ受け取る形にしてある）。
#
# 二乗は x*x で書く。べき乗 x**2 は 1 ULP 違うことがあり、コンパイルした版で
# 再現できないため（numba の x**2 / math.pow / np.power はいずれも乗算になる）。
# 両者が一致することは tests/test_thermodynamics.py で固定している。


@kernel
def get_speedofsound(specfic_heat_ratio, density, pressure):
  # Speed of sound for ideal gas: c = sqrt( gamma*p/rho )

  return np.sqrt( specfic_heat_ratio*pressure/density )


@kernel
def get_enthalpy(specific_heat_volum, density, temperature, uvel, vvel, wvel, pressure):
  # Specific total enthalpy: H = Cv*T + 0.5*U^2 + p/rho

  return specific_heat_volum*temperature \
       + 0.50*( uvel*uvel + vvel*vvel + wvel*wvel ) \
       + pressure/density


@kernel
def get_max_eigenvalue(specfic_heat_ratio, dens, uvel, vvel, wvel, pres,
                       viscosity, length, vecx, vecy, vecz):
  """
  Maximum eigenvalue of the flux Jacobian on a cell interface: lambda=|u.n|+c+2*mu/(rho*d)

  非粘性では固有値が |u.n|+c, |u.n|-c, |u.n| なので最大値は |u.n|+c となる。
  これに粘性の寄与 2*mu/(rho*d) を加える。lambda は |u.n| を用いるため法線の向きには依らない。

  time_integration/eigenvalue.py と同じ式で、こちらはカーネルから呼ぶための版。
  """

  cvel = uvel*vecx + vvel*vecy + wvel*vecz
  sos  = np.sqrt( specfic_heat_ratio*pres/dens )

  return abs(cvel) + sos + 2.0*viscosity/(dens*length)
