#!/usr/bin/env python3

# Program to accumulate the LU-SGS diagonal over faces (scalar dissipation)

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/10

from slow.general import thermodynamics
from slow.general.jit import kernel


@kernel
def accumulate_diagonal_scalar(num_face, num_face_bd, num_primitiv,
                               face2cell, face2cell_bd, virtualcell_bd,
                               area_vec, area_vec_bd, length, length_bd,
                               specfic_heat_ratio, var_primitiv, var_primitiv_bd,
                               viscosity, viscosity_bd, var_diagonal):
  """
  Face contribution to the LU-SGS diagonal with the scalar dissipation: sum of lambda*S

  内部面では面の左右セルで評価した max(lambda_a, lambda_b) を採る。lusgs_sweep も
  同じ値を使うので D+L+U が整合した分解になる。面上の平均状態を使うと衝撃波近傍で
  散逸が不足し、内部反復の収束が悪化する（ノズル CFL200 で 1.6 倍悪化。
  max なら逆に 2.4 倍改善する）。

  境界面は仮想セルの有無に応じた界面の状態で評価する。

  時間項（V/dt など）は設定で分岐するので、呼び出し側のセルループで足す。
  """

  # Inner faces
  for n_face in range(0,num_face):
    area = area_vec[0,n_face]
    vecx = area_vec[1,n_face]
    vecy = area_vec[2,n_face]
    vecz = area_vec[3,n_face]

    leng_a = length[0,n_face]
    leng_b = length[1,n_face]

    n_cell_a = face2cell[0,n_face]
    n_cell_b = face2cell[1,n_face]

    lenght_tmp = leng_a + leng_b
    eigenvalue = max( thermodynamics.get_max_eigenvalue(specfic_heat_ratio,        \
                        var_primitiv[0,n_cell_a], var_primitiv[1,n_cell_a],        \
                        var_primitiv[2,n_cell_a], var_primitiv[3,n_cell_a],        \
                        var_primitiv[5,n_cell_a], viscosity[n_cell_a],             \
                        lenght_tmp, vecx, vecy, vecz),                             \
                      thermodynamics.get_max_eigenvalue(specfic_heat_ratio,        \
                        var_primitiv[0,n_cell_b], var_primitiv[1,n_cell_b],        \
                        var_primitiv[2,n_cell_b], var_primitiv[3,n_cell_b],        \
                        var_primitiv[5,n_cell_b], viscosity[n_cell_b],             \
                        lenght_tmp, vecx, vecy, vecz) )

    var_diagonal[n_cell_a] = var_diagonal[n_cell_a] + eigenvalue*area
    var_diagonal[n_cell_b] = var_diagonal[n_cell_b] + eigenvalue*area

  # Boundary faces
  for n_face in range(0,num_face_bd):
    area = area_vec_bd[0,n_face]
    vecx = area_vec_bd[1,n_face]
    vecy = area_vec_bd[2,n_face]
    vecz = area_vec_bd[3,n_face]

    vcell_bd = virtualcell_bd[n_face]

    leng_a = length_bd[n_face]
    leng_b = float(vcell_bd)*leng_a

    n_cell_a = face2cell_bd[0,n_face]

    # Values on cell interface
    dens = float(vcell_bd)*0.50*( var_primitiv[0,n_cell_a] + var_primitiv_bd[0,n_face] ) \
         + float(1-vcell_bd)*var_primitiv_bd[0,n_face]
    uvel = float(vcell_bd)*0.50*( var_primitiv[1,n_cell_a] + var_primitiv_bd[1,n_face] ) \
         + float(1-vcell_bd)*var_primitiv_bd[1,n_face]
    vvel = float(vcell_bd)*0.50*( var_primitiv[2,n_cell_a] + var_primitiv_bd[2,n_face] ) \
         + float(1-vcell_bd)*var_primitiv_bd[2,n_face]
    wvel = float(vcell_bd)*0.50*( var_primitiv[3,n_cell_a] + var_primitiv_bd[3,n_face] ) \
         + float(1-vcell_bd)*var_primitiv_bd[3,n_face]
    pres = float(vcell_bd)*0.50*( var_primitiv[5,n_cell_a] + var_primitiv_bd[5,n_face] ) \
         + float(1-vcell_bd)*var_primitiv_bd[5,n_face]

    lenght_tmp = leng_a + leng_b
    visc_tmp   = float(vcell_bd)*0.50*( viscosity[n_cell_a] + viscosity_bd[n_face] ) \
               + float(1-vcell_bd)*viscosity_bd[n_face]
    eigenvalue = thermodynamics.get_max_eigenvalue(specfic_heat_ratio, dens, uvel, vvel, wvel, \
                                                   pres, visc_tmp, lenght_tmp, vecx, vecy, vecz)

    var_diagonal[n_cell_a] = var_diagonal[n_cell_a] + eigenvalue*area

  return var_diagonal
