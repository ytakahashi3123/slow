#!/usr/bin/env python3

# Program to sweep Jacobian matrix for LU-SGS routine

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/31

import numpy as np
from slow.orbital.orbital import orbital
from slow.time_integration import flux_jacobian
from slow.time_integration import lusgs
from slow.time_integration import lusgs_sweep_kernel
from slow.time_integration import roe_dissipation

@orbital.time_measurement_decorated
def sweep_jacobian(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_conserv, var_diagonal, var_dq):

  # Main routine

  # 行列散逸を選んだときは対角がブロック行列になるので別ルーチンで扱う
  if lusgs.get_kind_dissipation(config) == lusgs.KIND_DISSIPATION_MATRIX :
    return sweep_jacobian_matrix(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, \
                                 transport_coefficient_dict, var_primitiv, var_conserv, \
                                 var_diagonal, var_dq)

  # Input parameters。面ループの本体は lusgs_sweep_kernel に置いてある
  # （numba を掛けるため、配列とスカラーだけを引数に取る素の関数にしてある）
  num_conserv = dimension_dict['num_conservative']

  num_face    = geom_dict['num_face_inner']
  face2cell   = geom_dict['face2cell_inner']

  area_vec    = metrics_dict['area_vec_inner']
  length      = metrics_dict['length_inner']

  specfic_heat_ratio  = gas_property_dict['specfic_heat_ratio']
  specific_heat_volum = gas_property_dict['specific_heat_volume']

  viscosity           = transport_coefficient_dict['viscosity']

  var_dq = lusgs_sweep_kernel.sweep_scalar(num_face, num_conserv, face2cell, area_vec, length,
                                           specfic_heat_ratio, specific_heat_volum,
                                           var_primitiv, viscosity, var_diagonal, var_dq)

  return var_dq


@orbital.time_measurement_decorated
def sweep_jacobian_matrix(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_conserv, var_diagonal, var_dq):

  # Sweep with the matrix dissipation |A_Roe|
  #
  # 非対角ブロックは自セルの外向き法線に対する A^- である。
  #   前進 (b を a から更新): 0.5*( A(u_a,+n) - |A_Roe| )*S   （n は b の外向き）
  #   後退 (a を b から更新): 0.5*( A(u_b,-n) - |A_Roe| )*S   （-n は a の外向き）
  # var_diagonal には lusgs_diagonal.get_diagonal_matrix が D の逆行列を入れてある。

  # Input parameters
  num_conserv  = dimension_dict['num_conservative']

  num_face    = geom_dict['num_face_inner']
  face2cell   = geom_dict['face2cell_inner']

  area_vec    = metrics_dict['area_vec_inner']
  volume      = metrics_dict['volume_cell']

  specfic_heat_ratio  = gas_property_dict['specfic_heat_ratio']
  specific_heat_volum = gas_property_dict['specific_heat_volume']


  # Initialize
  jacobian    = np.zeros((num_conserv, num_conserv))
  absjacobian = np.zeros((num_conserv, num_conserv))
  pmatrix     = np.zeros((num_conserv, num_conserv))


  # Forward sweep
  for n_face in range(0,num_face):

    area = area_vec[0,n_face]
    vecx = area_vec[1,n_face]
    vecy = area_vec[2,n_face]
    vecz = area_vec[3,n_face]

    n_cell_a = face2cell[0,n_face]
    n_cell_b = face2cell[1,n_face]

    prim_a = var_primitiv[:,n_cell_a]
    prim_b = var_primitiv[:,n_cell_b]

    roe_dissipation.get_face_dissipation(specfic_heat_ratio, specific_heat_volum, prim_a, prim_b, \
                                         vecx, vecy, vecz, absjacobian, pmatrix)

    # A^- for cell b: 隣接セル a の状態、b の外向き法線 +n
    enth_a = orbital.get_enthalpy("self", specific_heat_volum, prim_a[0], prim_a[4], prim_a[1:4], prim_a[5])
    cvel_a = prim_a[1]*vecx + prim_a[2]*vecy + prim_a[3]*vecz
    flux_jacobian.set_flux_jacobian(jacobian, specfic_heat_ratio, 0.0, cvel_a, \
                                    prim_a[1], prim_a[2], prim_a[3], enth_a, vecx, vecy, vecz)

    offdiagonal = 0.50*( jacobian - absjacobian )*area
    var_dq[:,n_cell_b] = var_dq[:,n_cell_b] \
                       - var_diagonal[:,:,n_cell_b] @ ( offdiagonal @ var_dq[:,n_cell_a] )


  # Backward sweep
  for n_face_b in range(0,num_face):
    n_face = num_face-n_face_b-1

    area = area_vec[0,n_face]
    vecx = area_vec[1,n_face]
    vecy = area_vec[2,n_face]
    vecz = area_vec[3,n_face]

    n_cell_a = face2cell[0,n_face]
    n_cell_b = face2cell[1,n_face]

    prim_a = var_primitiv[:,n_cell_a]
    prim_b = var_primitiv[:,n_cell_b]

    roe_dissipation.get_face_dissipation(specfic_heat_ratio, specific_heat_volum, prim_a, prim_b, \
                                         vecx, vecy, vecz, absjacobian, pmatrix)

    # A^- for cell a: 隣接セル b の状態、a の外向き法線 -n
    enth_b = orbital.get_enthalpy("self", specific_heat_volum, prim_b[0], prim_b[4], prim_b[1:4], prim_b[5])
    cvel_b = -( prim_b[1]*vecx + prim_b[2]*vecy + prim_b[3]*vecz )
    flux_jacobian.set_flux_jacobian(jacobian, specfic_heat_ratio, 0.0, cvel_b, \
                                    prim_b[1], prim_b[2], prim_b[3], enth_b, -vecx, -vecy, -vecz)

    offdiagonal = 0.50*( jacobian - absjacobian )*area
    var_dq[:,n_cell_a] = var_dq[:,n_cell_a] \
                       - var_diagonal[:,:,n_cell_a] @ ( offdiagonal @ var_dq[:,n_cell_b] )


  return var_dq
