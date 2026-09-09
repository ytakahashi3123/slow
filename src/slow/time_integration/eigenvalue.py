#!/usr/bin/env python3

# Program to calculate the maximum eigenvalue of the flux Jacobian on a cell interface

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

from slow.orbital.orbital import orbital


def get_max_eigenvalue(specfic_heat_ratio, dens, uvel, vvel, wvel, pres, viscosity, length, vecx, vecy, vecz):
  """
  Maximum eigenvalue of the flux Jacobian on a cell interface: lambda=|u.n|+c+2*mu/(rho*d)

  非粘性では固有値が |u.n|+c, |u.n|-c, |u.n| なので最大値は |u.n|+c となる。
  これに粘性の寄与 2*mu/(rho*d) を加える。

  lambda の式そのものを 1 か所に集約するための関数。lusgs_diagonal（対角項）と
  lusgs_sweep（非対角項）の双方から呼ぶ。

  内部面では両者が同じ値を使う。すなわち面の左右セルで評価した max(lambda_a, lambda_b) であり、
  D=V/dt+0.5*sum(lambda*S) と 0.5*(A-lambda*I)*S が整合した D+L+U の分解になる。
  なお面上の「平均」状態で lambda を評価する案は、ノズル・球の両ケース、CFL 2.5/50 の
  いずれでも内部反復の減衰が 0--3% 悪化したため採っていない（衝撃波近傍で散逸が不足する）。
  境界面 (lusgs_diagonal のみ) は仮想セルの有無に応じた界面の状態で評価する。

  --dens, uvel, vvel, wvel, pres: lambda を評価する状態
  --viscosity: 同じ位置での粘性係数
  --length: セル中心間距離
  --vecx, vecy, vecz: 面の単位法線ベクトル。lambda は |u.n| を用いるため向きには依らない
  """

  # Contravariant velocity
  cvel = uvel*vecx + vvel*vecy + wvel*vecz

  # Thermodynamic properties: speed of sound
  sos  = orbital.get_speedofsound("self", specfic_heat_ratio, dens, pres)

  return abs(cvel) + sos + 2.0*viscosity/(dens*length)
