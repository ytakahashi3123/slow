#!/usr/bin/env python3

# Program to sweep Jacobian matrix for LU-SGS routine

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/31

import numpy as np
from slow.orbital.orbital import orbital
from slow.time_integration import eigenvalue as eigenvalue_mod
from slow.time_integration import flux_jacobian
from slow.time_integration import lusgs
from slow.time_integration import roe_dissipation

@orbital.time_measurement_decorated
def get_diagonal(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_primitiv_bd, var_conserv, var_conserv_prev, var_rhs, var_dt, var_diagonal, var_dq, num_conserv_prev_level=lusgs.NUM_PREV_LEVEL_REQUIRED_BDF2):

  # Main routine

  # 行列散逸を選んだときは対角がブロック行列になるので別ルーチンで扱う
  if lusgs.get_kind_dissipation(config) == lusgs.KIND_DISSIPATION_MATRIX :
    return get_diagonal_matrix(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, \
                               transport_coefficient_dict, var_primitiv, var_primitiv_bd, \
                               var_conserv, var_conserv_prev, var_rhs, var_dt, var_diagonal, var_dq, \
                               num_conserv_prev_level)

  # Input parameters
  num_conserv    = dimension_dict['num_conservative']

  num_face       = geom_dict['num_face_inner']
  num_face_bd    = geom_dict['num_face_boundary']
  num_cell       = geom_dict['num_cell']
  face2cell      = geom_dict['face2cell_inner']
  face2cell_bd   = geom_dict['face2cell_boundary']
  virtualcell_bd = geom_dict['virtualcell_boundary']

  area_vec    = metrics_dict['area_vec_inner']
  area_vec_bd = metrics_dict['area_vec_boundary']
  length      = metrics_dict['length_inner']
  length_bd   = metrics_dict['length_boundary']
  volume      = metrics_dict['volume_cell']

  specfic_heat_ratio  = gas_property_dict['specfic_heat_ratio']

  viscosity       = transport_coefficient_dict['viscosity']
  viscosity_bd    = transport_coefficient_dict['viscosity_boundary']


  # Initialize
  #var_dq[:,:] = 0.0
  #var_diagonal[:] = 0.0


  # Inner faces
  for n_face in range(0,num_face):

    # Face area vector
    area = area_vec[0,n_face]
    vecx = area_vec[1,n_face]
    vecy = area_vec[2,n_face]
    vecz = area_vec[3,n_face]
    # Length
    leng_a = length[0,n_face]
    leng_b = length[1,n_face]


    # Cell ID
    # --from self cell side
    n_cell_a = face2cell[0,n_face]
    # --from neigboring cell
    n_cell_b = face2cell[1,n_face]


    # Primitive variables
    prim_a = var_primitiv[:,n_cell_a]
    prim_b = var_primitiv[:,n_cell_b]

    # Maximum eigenvalue of Jacobian matrix
    # --面の左右セルで評価して大きい方を採る。
    #   lusgs_sweep も同じ値を使うので D+L+U が整合した分解になる。
    #   面上の平均状態を使うと衝撃波近傍で散逸が不足し、内部反復の収束が悪化する
    #   （ノズル CFL200 で 1.6 倍悪化。max なら逆に 2.4 倍改善する）
    lenght_tmp = leng_a + leng_b
    eigenvalue = max( eigenvalue_mod.get_max_eigenvalue(specfic_heat_ratio, \
                                                        prim_a[0], prim_a[1], prim_a[2], prim_a[3], prim_a[5], \
                                                        viscosity[n_cell_a], lenght_tmp, vecx, vecy, vecz), \
                      eigenvalue_mod.get_max_eigenvalue(specfic_heat_ratio, \
                                                        prim_b[0], prim_b[1], prim_b[2], prim_b[3], prim_b[5], \
                                                        viscosity[n_cell_b], lenght_tmp, vecx, vecy, vecz) )

    # Set diagonal values
    var_diagonal[n_cell_a] = var_diagonal[n_cell_a] + eigenvalue*area
    var_diagonal[n_cell_b] = var_diagonal[n_cell_b] + eigenvalue*area


  # Boundary faces
  for n_face in range(0,num_face_bd):

    # Face area vector
    area = area_vec_bd[0,n_face]
    vecx = area_vec_bd[1,n_face]
    vecy = area_vec_bd[2,n_face]
    vecz = area_vec_bd[3,n_face]

    # Virtual cell identificaton on boudary
    vcell_bd = virtualcell_bd[n_face]

    # Length
    leng_a = length_bd[n_face]
    leng_b = float(vcell_bd)*leng_a

    # Cell ID
    # --from self cell side
    n_cell_a = face2cell_bd[0,n_face]
    # Primitive variables
    prim_a = var_primitiv[:,n_cell_a]
    prim_b = var_primitiv_bd[:,n_face]

    # Values on cell interface
    prim   = float(vcell_bd)*0.50*( prim_a + prim_b ) + float(1-vcell_bd)*prim_b

    # Maximum eigenvalue of Jacobian matrix
    lenght_tmp = leng_a + leng_b
    visc_tmp   = float(vcell_bd)*0.50*( viscosity[n_cell_a] + viscosity_bd[n_face] ) + float(1-vcell_bd)*viscosity_bd[n_face]
    eigenvalue = eigenvalue_mod.get_max_eigenvalue(specfic_heat_ratio, \
                                                   prim[0], prim[1], prim[2], prim[3], prim[5], \
                                                   visc_tmp, lenght_tmp, vecx, vecy, vecz)

    # Set diagonal values
    var_diagonal[n_cell_a] = var_diagonal[n_cell_a] + eigenvalue*area


  #var_dq = var_rhs


  # Diagonal element for factoriization matrix and delta Q
  for n_cell in range(0,num_cell):
    diag_time, dq_unst, face_scale = lusgs.get_time_term(config, n_cell, volume, var_dt, \
                                                         var_conserv, var_conserv_prev, \
                                                         num_conserv_prev_level)
    var_diagonal[n_cell] = diag_time + face_scale*var_diagonal[n_cell]
    var_dq[:,n_cell]     = ( -var_rhs[:,n_cell]-dq_unst )/var_diagonal[n_cell]


  return var_diagonal, var_dq


@orbital.time_measurement_decorated
def get_diagonal_matrix(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_primitiv_bd, var_conserv, var_conserv_prev, var_rhs, var_dt, var_diagonal, var_dq, num_conserv_prev_level=lusgs.NUM_PREV_LEVEL_REQUIRED_BDF2):

  # Block-diagonal version used with the matrix dissipation |A_Roe|
  #
  # 面 (a,b) について、area_vec の法線 n はセル b の外向き（セル a の内向き）である。
  #   D_a += 0.5*( A(u_a,-n) + |A_Roe| )*S
  #   D_b += 0.5*( A(u_b,+n) + |A_Roe| )*S
  # これは面の流束線形化から定まる A^+ である。
  #
  # var_diagonal は (num_conserv, num_conserv, num_cell) 形状で、
  # 最後に D の「逆行列」を格納する（スイープ側を行列ベクトル積だけで済ませるため）。

  # Input parameters
  num_conserv    = dimension_dict['num_conservative']

  num_face       = geom_dict['num_face_inner']
  num_face_bd    = geom_dict['num_face_boundary']
  num_cell       = geom_dict['num_cell']
  face2cell      = geom_dict['face2cell_inner']
  face2cell_bd   = geom_dict['face2cell_boundary']

  area_vec    = metrics_dict['area_vec_inner']
  area_vec_bd = metrics_dict['area_vec_boundary']
  volume      = metrics_dict['volume_cell']

  specfic_heat_ratio  = gas_property_dict['specfic_heat_ratio']
  specific_heat_volum = gas_property_dict['specific_heat_volume']


  # Initialize (面ループ内で確保しないよう作業配列を用意しておく)
  jacobian    = np.zeros((num_conserv, num_conserv))
  absjacobian = np.zeros((num_conserv, num_conserv))
  pmatrix     = np.zeros((num_conserv, num_conserv))
  identity    = np.identity(num_conserv)


  # Inner faces
  for n_face in range(0,num_face):

    area = area_vec[0,n_face]
    vecx = area_vec[1,n_face]
    vecy = area_vec[2,n_face]
    vecz = area_vec[3,n_face]

    n_cell_a = face2cell[0,n_face]
    n_cell_b = face2cell[1,n_face]

    prim_a = var_primitiv[:,n_cell_a]
    prim_b = var_primitiv[:,n_cell_b]

    # Matrix dissipation on the face (法線の向きには依らない)
    roe_dissipation.get_face_dissipation(specfic_heat_ratio, specific_heat_volum, prim_a, prim_b, \
                                         vecx, vecy, vecz, absjacobian, pmatrix)

    # A^+ for cell a: 外向き法線は -n
    enth_a = orbital.get_enthalpy("self", specific_heat_volum, prim_a[0], prim_a[4], prim_a[1:4], prim_a[5])
    cvel_a = -( prim_a[1]*vecx + prim_a[2]*vecy + prim_a[3]*vecz )
    flux_jacobian.set_flux_jacobian(jacobian, specfic_heat_ratio, 0.0, cvel_a, \
                                    prim_a[1], prim_a[2], prim_a[3], enth_a, -vecx, -vecy, -vecz)
    var_diagonal[:,:,n_cell_a] = var_diagonal[:,:,n_cell_a] + ( jacobian + absjacobian )*area

    # A^+ for cell b: 外向き法線は +n
    enth_b = orbital.get_enthalpy("self", specific_heat_volum, prim_b[0], prim_b[4], prim_b[1:4], prim_b[5])
    cvel_b = prim_b[1]*vecx + prim_b[2]*vecy + prim_b[3]*vecz
    flux_jacobian.set_flux_jacobian(jacobian, specfic_heat_ratio, 0.0, cvel_b, \
                                    prim_b[1], prim_b[2], prim_b[3], enth_b, vecx, vecy, vecz)
    var_diagonal[:,:,n_cell_b] = var_diagonal[:,:,n_cell_b] + ( jacobian + absjacobian )*area


  # Boundary faces
  # --境界では仮想セルの有無で状態の重みを変えず、セルと境界値の Roe 平均をそのまま使う
  for n_face in range(0,num_face_bd):

    area = area_vec_bd[0,n_face]
    vecx = area_vec_bd[1,n_face]
    vecy = area_vec_bd[2,n_face]
    vecz = area_vec_bd[3,n_face]

    n_cell_a = face2cell_bd[0,n_face]

    prim_a = var_primitiv[:,n_cell_a]
    prim_b = var_primitiv_bd[:,n_face]

    roe_dissipation.get_face_dissipation(specfic_heat_ratio, specific_heat_volum, prim_a, prim_b, \
                                         vecx, vecy, vecz, absjacobian, pmatrix)

    enth_a = orbital.get_enthalpy("self", specific_heat_volum, prim_a[0], prim_a[4], prim_a[1:4], prim_a[5])
    cvel_a = -( prim_a[1]*vecx + prim_a[2]*vecy + prim_a[3]*vecz )
    flux_jacobian.set_flux_jacobian(jacobian, specfic_heat_ratio, 0.0, cvel_a, \
                                    prim_a[1], prim_a[2], prim_a[3], enth_a, -vecx, -vecy, -vecz)
    var_diagonal[:,:,n_cell_a] = var_diagonal[:,:,n_cell_a] + ( jacobian + absjacobian )*area


  # Diagonal block for factorization and delta Q
  for n_cell in range(0,num_cell):
    diag_time, dq_unst, face_scale = lusgs.get_time_term(config, n_cell, volume, var_dt, \
                                                         var_conserv, var_conserv_prev, \
                                                         num_conserv_prev_level)
    diagonal_block = diag_time*identity + face_scale*var_diagonal[:,:,n_cell]
    # 以降のスイープを行列ベクトル積だけで済ませるため、ここで逆行列にしておく
    var_diagonal[:,:,n_cell] = np.linalg.inv(diagonal_block)
    var_dq[:,n_cell]         = var_diagonal[:,:,n_cell] @ ( -var_rhs[:,n_cell]-dq_unst )


  return var_diagonal, var_dq
