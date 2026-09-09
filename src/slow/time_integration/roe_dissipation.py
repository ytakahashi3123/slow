#!/usr/bin/env python3

# Program to calculate the matrix dissipation |A_Roe| = P|Lambda|P^-1 on a cell interface
#
# LU-SGS の非対角項 0.5*(A-|A_Roe|)*S に使う行列型の散逸。
# スカラー散逸 lambda*I（time_integration/eigenvalue.py）が最大固有値で全成分を
# 一律に潰すのに対し、こちらは成分ごとに固有値の絶対値ぶんだけ散逸させる。
# Roe 型の行列散逸である。
#
# 固有ベクトル行列 P は圧縮性 Euler 方程式の右固有ベクトル行列である。
# SLOW は 2 次元でも保存変数を 5 個持つので、4x4 の 2 次元形式ではなく 5x5 の
# 3 次元形式を使う。P の逆行列は numpy で数値的に求める（面あたり 5x5 なので、
# 閉形式を書き下すより取り違えの危険が小さい）。
#
# config の kind_lusgs_dissipation: matrix を指定したときだけ使われる。
# 面ごとに 5x5 の逆行列と行列積が入るが、純 Python の面ループ自体が重いため
# 実測では内部反復あたり 1.6 倍のコストに収まる（ノズル格子 7734 要素）。
# 一方で内部反復の収束は速くなるので、総合では有利になりうる。

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

import numpy as np

from slow.orbital.orbital import orbital


def get_roe_average(specfic_heat_ratio, dens_a, vel_a, enth_a, dens_b, vel_b, enth_b):
  """
  Roe average between the two cells sharing a face

  重みは Roe 平均の標準形 R=sqrt(rho_b/rho_a) を用いる。

  --vel_a, vel_b: 速度ベクトル（長さ 3）
  --enth_a, enth_b: 全エンタルピー H
  戻り値: (密度, 速度ベクトル, 音速)
  """

  fact = np.sqrt( abs(dens_b/dens_a) )

  dens_roe = fact*dens_a
  vel_roe  = ( fact*np.asarray(vel_b, dtype=float) + np.asarray(vel_a, dtype=float) )/(fact + 1.0)
  enth_roe = ( fact*enth_b + enth_a )/(fact + 1.0)

  # 音速は Roe 平均の全エンタルピーと運動エネルギーから求める: c^2=(gamma-1)*(H-q)
  sqvel_roe = np.dot(vel_roe, vel_roe)
  sos_roe   = np.sqrt( abs( (specfic_heat_ratio-1.0)*(enth_roe - 0.50*sqvel_roe) ) )

  return dens_roe, vel_roe, sos_roe


def set_pmatrix(pmatrix, specfic_heat_ratio, dens, vel, sos, vecx, vecy, vecz):
  """
  Set the matrix of right eigenvectors P of A(n)

  圧縮性 Euler 方程式の右固有ベクトル行列（3 次元形式）。列は固有値
  (u.n, u.n, u.n, u.n+c, u.n-c) の順に対応する。

  --pmatrix: (5,5) の配列。戻り値ではなくこれを書き換える
  --vel: 速度ベクトル（長さ 3）
  --vecx, vecy, vecz: 面の単位法線ベクトル
  """

  uvel, vvel, wvel = vel[0], vel[1], vel[2]
  sqvel  = uvel*uvel + vvel*vvel + wvel*wvel

  rhooc  = dens/sos
  rhoxc  = dens*sos
  gm1    = specfic_heat_ratio - 1.0

  pmatrix[0,0] = vecx
  pmatrix[0,1] = vecy
  pmatrix[0,2] = vecz
  pmatrix[0,3] = 0.50*rhooc
  pmatrix[0,4] = 0.50*rhooc

  pmatrix[1,0] = uvel*vecx
  pmatrix[1,1] = uvel*vecy - dens*vecz
  pmatrix[1,2] = uvel*vecz + dens*vecy
  pmatrix[1,3] = 0.50*( uvel*rhooc + dens*vecx )
  pmatrix[1,4] = 0.50*( uvel*rhooc - dens*vecx )

  pmatrix[2,0] = vvel*vecx + dens*vecz
  pmatrix[2,1] = vvel*vecy
  pmatrix[2,2] = vvel*vecz - dens*vecx
  pmatrix[2,3] = 0.50*( vvel*rhooc + dens*vecy )
  pmatrix[2,4] = 0.50*( vvel*rhooc - dens*vecy )

  pmatrix[3,0] = wvel*vecx - dens*vecy
  pmatrix[3,1] = wvel*vecy + dens*vecx
  pmatrix[3,2] = wvel*vecz
  pmatrix[3,3] = 0.50*( wvel*rhooc + dens*vecz )
  pmatrix[3,4] = 0.50*( wvel*rhooc - dens*vecz )

  cvel = uvel*vecx + vvel*vecy + wvel*vecz

  pmatrix[4,0] = 0.50*sqvel*vecx + dens*( vvel*vecz - wvel*vecy )
  pmatrix[4,1] = 0.50*sqvel*vecy + dens*( wvel*vecx - uvel*vecz )
  pmatrix[4,2] = 0.50*sqvel*vecz + dens*( uvel*vecy - vvel*vecx )
  pmatrix[4,3] = 0.50*( 0.50*sqvel*rhooc + dens*cvel + rhoxc/gm1 )
  pmatrix[4,4] = 0.50*( 0.50*sqvel*rhooc - dens*cvel + rhoxc/gm1 )

  return


def get_eigenvalues(vel, sos, vecx, vecy, vecz):
  """Eigenvalues of A(n) in the column order used by set_pmatrix"""

  cvel = vel[0]*vecx + vel[1]*vecy + vel[2]*vecz

  return np.array([cvel, cvel, cvel, cvel + sos, cvel - sos])


def set_absolute_jacobian(absjacobian, specfic_heat_ratio, dens, vel, sos, vecx, vecy, vecz,
                          pmatrix=None):
  """
  Set |A(n)| = P|Lambda|P^-1

  法線の向きには依らない（|A(-n)|=|-A(n)|=|A(n)|）。

  --absjacobian: (5,5) の配列。戻り値ではなくこれを書き換える
  --pmatrix: 作業用の (5,5) 配列。面ループから呼ぶときに確保を避けるため渡せるようにしている
  """

  if pmatrix is None:
    pmatrix = np.zeros((5,5))

  set_pmatrix(pmatrix, specfic_heat_ratio, dens, vel, sos, vecx, vecy, vecz)
  eigenvalues = np.abs( get_eigenvalues(vel, sos, vecx, vecy, vecz) )

  # |A| = P diag(|Lambda|) P^-1
  absjacobian[:,:] = pmatrix @ ( eigenvalues[:,None] * np.linalg.inv(pmatrix) )

  return


def get_face_dissipation(specfic_heat_ratio, specific_heat_volum, prim_a, prim_b,
                         vecx, vecy, vecz, absjacobian, pmatrix=None):
  """
  Set |A_Roe| on a face from the primitive variables of the two adjacent cells

  --prim_a, prim_b: 原始変数 (rho, u, v, w, T, p)
  """

  enth_a = orbital.get_enthalpy("self", specific_heat_volum, prim_a[0], prim_a[4], prim_a[1:4], prim_a[5])
  enth_b = orbital.get_enthalpy("self", specific_heat_volum, prim_b[0], prim_b[4], prim_b[1:4], prim_b[5])

  dens_roe, vel_roe, sos_roe = get_roe_average(specfic_heat_ratio, \
                                               prim_a[0], prim_a[1:4], enth_a, \
                                               prim_b[0], prim_b[1:4], enth_b)

  set_absolute_jacobian(absjacobian, specfic_heat_ratio, dens_roe, vel_roe, sos_roe, \
                        vecx, vecy, vecz, pmatrix)

  return
