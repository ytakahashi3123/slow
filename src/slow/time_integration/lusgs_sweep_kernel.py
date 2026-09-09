#!/usr/bin/env python3

# Program to sweep the LU-SGS implicit operator with the scalar dissipation

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/10

import numpy as np

from slow.general import thermodynamics
from slow.general.jit import kernel
from slow.time_integration.flux_jacobian import set_flux_jacobian


@kernel
def sweep_scalar(num_face, num_conserv, face2cell, area_vec, length,
                 specfic_heat_ratio, specific_heat_volum,
                 var_primitiv, viscosity, var_diagonal, var_dq):
  """
  Forward and backward Gauss-Seidel sweeps with the scalar dissipation

  前進で (D+L)x* = b を、後退で (D+U)x = D x* を解く。面ループの本体で、
  numba は dict を受け取れないので値の取り出しは呼び出し側で行う。

  非対角項 A^- が掛かるのは隣接セルの dq なので、ヤコビアンは隣接セル側の状態で
  評価する。面の法線 area_vec[1:4] はセル b の外向き（セル a の内向き）なので、
  後退スイープでは法線を反転してセル a の外向きに対する A^- を得る。

  前進が厳密な下三角 Gauss-Seidel になるのは、面が第 0 セル番号の昇順に並んで
  いるおかげである（face2cell[0,n] < face2cell[1,n] かつ第 0 行が単調非減少）。
  """

  jacobian = np.zeros((num_conserv, num_conserv))

  # Forward sweep
  for n_face in range(0,num_face):

    area = area_vec[0,n_face]
    vecx = area_vec[1,n_face]
    vecy = area_vec[2,n_face]
    vecz = area_vec[3,n_face]

    n_cell_a = face2cell[0,n_face]
    n_cell_b = face2cell[1,n_face]

    # --隣接セル (a) 側の状態でヤコビアンを評価する
    dens = var_primitiv[0,n_cell_a]
    uvel = var_primitiv[1,n_cell_a]
    vvel = var_primitiv[2,n_cell_a]
    wvel = var_primitiv[3,n_cell_a]
    temp = var_primitiv[4,n_cell_a]
    pres = var_primitiv[5,n_cell_a]

    cvel = uvel*vecx + vvel*vecy + wvel*vecz
    enth = thermodynamics.get_enthalpy(specific_heat_volum, dens, temp, uvel, vvel, wvel, pres)

    # --対角項 (lusgs_diagonal) と同一の値でなければ D+L+U が整合した分解にならないので、
    #   面の左右セルで評価して大きい方を採る
    lenght_tmp = length[0,n_face] + length[1,n_face]
    eigenvalue = max( thermodynamics.get_max_eigenvalue(specfic_heat_ratio,                  \
                        var_primitiv[0,n_cell_a], var_primitiv[1,n_cell_a],                  \
                        var_primitiv[2,n_cell_a], var_primitiv[3,n_cell_a],                  \
                        var_primitiv[5,n_cell_a], viscosity[n_cell_a],                       \
                        lenght_tmp, vecx, vecy, vecz),                                       \
                      thermodynamics.get_max_eigenvalue(specfic_heat_ratio,                  \
                        var_primitiv[0,n_cell_b], var_primitiv[1,n_cell_b],                  \
                        var_primitiv[2,n_cell_b], var_primitiv[3,n_cell_b],                  \
                        var_primitiv[5,n_cell_b], viscosity[n_cell_b],                       \
                        lenght_tmp, vecx, vecy, vecz) )

    set_flux_jacobian(jacobian, specfic_heat_ratio, eigenvalue, \
                      cvel, uvel, vvel, wvel, enth, vecx, vecy, vecz)

    var_inv = 0.50*area/var_diagonal[n_cell_b]
    for m in range(0,num_conserv):
      dq_tmp = 0.0
      for k in range(0,num_conserv):
        dq_tmp = dq_tmp + jacobian[m,k]*var_dq[k,n_cell_a]
      var_dq[m,n_cell_b] = var_dq[m,n_cell_b] - dq_tmp*var_inv

  # Backward sweep
  for n_face_b in range(0,num_face):
    n_face = num_face-n_face_b-1

    area = area_vec[0,n_face]
    vecx =-area_vec[1,n_face]
    vecy =-area_vec[2,n_face]
    vecz =-area_vec[3,n_face]

    n_cell_a = face2cell[0,n_face]
    n_cell_b = face2cell[1,n_face]

    # --隣接セル (b) 側の状態でヤコビアンを評価する
    dens = var_primitiv[0,n_cell_b]
    uvel = var_primitiv[1,n_cell_b]
    vvel = var_primitiv[2,n_cell_b]
    wvel = var_primitiv[3,n_cell_b]
    temp = var_primitiv[4,n_cell_b]
    pres = var_primitiv[5,n_cell_b]

    cvel = uvel*vecx + vvel*vecy + wvel*vecz
    enth = thermodynamics.get_enthalpy(specific_heat_volum, dens, temp, uvel, vvel, wvel, pres)

    lenght_tmp = length[0,n_face] + length[1,n_face]
    eigenvalue = max( thermodynamics.get_max_eigenvalue(specfic_heat_ratio,                  \
                        var_primitiv[0,n_cell_a], var_primitiv[1,n_cell_a],                  \
                        var_primitiv[2,n_cell_a], var_primitiv[3,n_cell_a],                  \
                        var_primitiv[5,n_cell_a], viscosity[n_cell_a],                       \
                        lenght_tmp, vecx, vecy, vecz),                                       \
                      thermodynamics.get_max_eigenvalue(specfic_heat_ratio,                  \
                        var_primitiv[0,n_cell_b], var_primitiv[1,n_cell_b],                  \
                        var_primitiv[2,n_cell_b], var_primitiv[3,n_cell_b],                  \
                        var_primitiv[5,n_cell_b], viscosity[n_cell_b],                       \
                        lenght_tmp, vecx, vecy, vecz) )

    set_flux_jacobian(jacobian, specfic_heat_ratio, eigenvalue, \
                      cvel, uvel, vvel, wvel, enth, vecx, vecy, vecz)

    var_inv = 0.50*area/var_diagonal[n_cell_a]
    for m in range(0,num_conserv):
      dq_tmp = 0.0
      for k in range(0,num_conserv):
        dq_tmp = dq_tmp + jacobian[m,k]*var_dq[k,n_cell_b]
      var_dq[m,n_cell_a] = var_dq[m,n_cell_a] - dq_tmp*var_inv

  return var_dq
