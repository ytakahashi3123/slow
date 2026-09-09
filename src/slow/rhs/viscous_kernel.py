#!/usr/bin/env python3

# Program to calculate the viscous flux on a face and accumulate it into the residual

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/10

import numpy as np

from slow.general.jit import kernel


@kernel
def get_stress_tensor(tau, work_stress, grad_face, viscosty_face,
                      uvel_face, vvel_face, wvel_face):
  """
  Stress tensor tau_ij = mu*( du_i/dx_j + du_j/dx_i - 2/3*div(u)*delta_ij ) and its work

  grad_face[i,m] は m 番目の原始変数の i 方向微分。速度は m=1,2,3。
  work_stress は tau と面上の速度の積で、エネルギー式の粘性項に入る。
  """

  tau_tmp  = 1.0/3.0*(grad_face[0,1] + grad_face[1,2] + grad_face[2,3])
  tau[0,0] = 2.0*viscosty_face*( grad_face[0,1] - tau_tmp )
  tau[0,1] =     viscosty_face*( grad_face[0,2] + grad_face[1,1] )
  tau[0,2] =     viscosty_face*( grad_face[0,3] + grad_face[2,1] )
  tau[1,0] = tau[0,1]
  tau[1,1] = 2.0*viscosty_face*( grad_face[1,2] - tau_tmp )
  tau[1,2] =     viscosty_face*( grad_face[1,3] + grad_face[2,2] )
  tau[2,0] = tau[0,2]
  tau[2,1] = tau[1,2]
  tau[2,2] = 2.0*viscosty_face*( grad_face[2,3] - tau_tmp )

  # Viscous stress work
  work_stress[0] = tau[0,0]*uvel_face + tau[0,1]*vvel_face + tau[0,2]*wvel_face
  work_stress[1] = tau[1,0]*uvel_face + tau[1,1]*vvel_face + tau[1,2]*wvel_face
  work_stress[2] = tau[2,0]*uvel_face + tau[2,1]*vvel_face + tau[2,2]*wvel_face

  return


@kernel
def set_viscous_flux(flux_rhs, tau, work_stress, heat_flux, vecx, vecy, vecz):
  # Viscous flux projected on the face normal

  flux_rhs[0] = 0.0
  flux_rhs[1] = - ( tau[0,0]*vecx + tau[0,1]*vecy + tau[0,2]*vecz )
  flux_rhs[2] = - ( tau[1,0]*vecx + tau[1,1]*vecy + tau[1,2]*vecz )
  flux_rhs[3] = - ( tau[2,0]*vecx + tau[2,1]*vecy + tau[2,2]*vecz )
  flux_rhs[4] = - ( (work_stress[0]+heat_flux[0])*vecx    \
                  + (work_stress[1]+heat_flux[1])*vecy    \
                  + (work_stress[2]+heat_flux[2])*vecz )

  return


@kernel
def accumulate_viscous(num_face, num_face_bd, num_conserv, num_primitiv,
                       face2cell, face2cell_bd, area_vec, area_vec_bd,
                       var_primitiv, var_primitiv_bd, var_gradient,
                       viscosity, thermal_cond, viscosity_bd, thermal_cond_bd, var_rhs):
  """
  Accumulate the viscous flux into the residual over all faces

  面上の勾配は内部面では左右セルの単純平均、境界面では隣接セルの値をそのまま使う。
  輸送係数も同様に相加平均で評価する。

  var_rhs の符号規約は移流項と同じ（法線はセル a に対して内向きなので a に引き b に足す）。
  """

  num_spatial = 3

  grad_face   = np.zeros((num_spatial, num_primitiv))
  tau         = np.zeros((num_spatial, num_spatial))
  work_stress = np.zeros(num_spatial)
  heat_flux   = np.zeros(num_spatial)
  flux_rhs    = np.zeros(num_conserv)

  # --Inner loop
  for n_face in range(0,num_face):

    area = area_vec[0,n_face]
    vecx = area_vec[1,n_face]
    vecy = area_vec[2,n_face]
    vecz = area_vec[3,n_face]

    n_cell_a = face2cell[0,n_face]
    n_cell_b = face2cell[1,n_face]

    # Gradient variables on face（単純な平均で評価している）
    for i in range(0,num_spatial):
      for m in range(0,num_primitiv):
        grad_face[i,m] = 0.50*(var_gradient[i,m,n_cell_a] + var_gradient[i,m,n_cell_b])

    # Transport coefficients on face
    # 熱伝導率は相加平均で評価したが、相乗平均でも良いかも
    viscosty_face     = 0.50*( viscosity[n_cell_a] + viscosity[n_cell_b] )
    thermal_cond_face = 0.50*( thermal_cond[n_cell_a] + thermal_cond[n_cell_b] )

    uvel_face = 0.50*(var_primitiv[1,n_cell_a] + var_primitiv[1,n_cell_b])
    vvel_face = 0.50*(var_primitiv[2,n_cell_a] + var_primitiv[2,n_cell_b])
    wvel_face = 0.50*(var_primitiv[3,n_cell_a] + var_primitiv[3,n_cell_b])

    get_stress_tensor(tau, work_stress, grad_face, viscosty_face,
                      uvel_face, vvel_face, wvel_face)

    # Heat flux: lambda*dT/dx
    heat_flux[0] = thermal_cond_face*grad_face[0,4]
    heat_flux[1] = thermal_cond_face*grad_face[1,4]
    heat_flux[2] = thermal_cond_face*grad_face[2,4]

    set_viscous_flux(flux_rhs, tau, work_stress, heat_flux, vecx, vecy, vecz)

    for m in range(0,num_conserv):
      var_rhs[m,n_cell_a] = var_rhs[m,n_cell_a] - flux_rhs[m]*area
      var_rhs[m,n_cell_b] = var_rhs[m,n_cell_b] + flux_rhs[m]*area

  # --Boundary loop
  for n_face in range(0,num_face_bd):

    area = area_vec_bd[0,n_face]
    vecx = area_vec_bd[1,n_face]
    vecy = area_vec_bd[2,n_face]
    vecz = area_vec_bd[3,n_face]

    n_cell_a = face2cell_bd[0,n_face]

    # Gradient variables on face（隣接セルの勾配をそのまま使う）
    for i in range(0,num_spatial):
      for m in range(0,num_primitiv):
        grad_face[i,m] = var_gradient[i,m,n_cell_a]

    viscosty_face     = 0.50*( viscosity[n_cell_a] + viscosity_bd[n_face] )
    thermal_cond_face = 0.50*( thermal_cond[n_cell_a] + thermal_cond_bd[n_face] )

    uvel_face = 0.50*(var_primitiv[1,n_cell_a] + var_primitiv_bd[1,n_face])
    vvel_face = 0.50*(var_primitiv[2,n_cell_a] + var_primitiv_bd[2,n_face])
    wvel_face = 0.50*(var_primitiv[3,n_cell_a] + var_primitiv_bd[3,n_face])

    get_stress_tensor(tau, work_stress, grad_face, viscosty_face,
                      uvel_face, vvel_face, wvel_face)

    heat_flux[0] = thermal_cond_face*grad_face[0,4]
    heat_flux[1] = thermal_cond_face*grad_face[1,4]
    heat_flux[2] = thermal_cond_face*grad_face[2,4]

    set_viscous_flux(flux_rhs, tau, work_stress, heat_flux, vecx, vecy, vecz)

    for m in range(0,num_conserv):
      var_rhs[m,n_cell_a] = var_rhs[m,n_cell_a] - flux_rhs[m]*area

  return var_rhs
