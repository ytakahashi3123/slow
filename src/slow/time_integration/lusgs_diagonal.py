#!/usr/bin/env python3

# Program to sweep Jacobian matrix for LU-SGS routine

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/31

import numpy as np
from slow.orbital.orbital import orbital
from slow.time_integration import flux_jacobian
from slow.time_integration import lusgs
from slow.time_integration import lusgs_diagonal_kernel
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


  # 面ループの本体は lusgs_diagonal_kernel に置いてある
  # （numba を掛けるため、配列とスカラーだけを引数に取る素の関数にしてある）
  var_diagonal = lusgs_diagonal_kernel.accumulate_diagonal_scalar(
                   num_face, num_face_bd, var_primitiv.shape[0],
                   face2cell, face2cell_bd, virtualcell_bd,
                   area_vec, area_vec_bd, length, length_bd,
                   specfic_heat_ratio, var_primitiv, var_primitiv_bd,
                   viscosity, viscosity_bd, var_diagonal)

  #var_dq = var_rhs


  # Diagonal element for factoriization matrix and delta Q
  # 設定は 1 回だけ解決し、時間刻みの正値性もまとめて確かめる
  kind_time, timestep_outer, face_scale = lusgs.get_time_setting(config, num_conserv_prev_level)
  lusgs.check_timestep_array(kind_time, var_dt, timestep_outer)

  var_diagonal, var_dq = lusgs_diagonal_kernel.apply_time_term_scalar(
                           kind_time, lusgs.KIND_TIME_STEADY, lusgs.KIND_TIME_BDF2,
                           num_cell, num_conserv, timestep_outer, face_scale,
                           volume, var_dt, var_conserv, var_conserv_prev, var_rhs,
                           var_diagonal, var_dq)


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
