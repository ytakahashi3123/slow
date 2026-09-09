#!/usr/bin/env python3

# Program to calculate the advection flux on a face and accumulate it into the residual

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/10

import numpy as np

from slow.general import thermodynamics
from slow.general.jit import kernel


# 移流スキームの選択。numba のカーネル内では文字列で分岐させないので整数で渡す
KIND_SCHEME_SLAU2  = 0
KIND_SCHEME_HAENEL = 1

# 二乗は x*x で書く（べき乗は 1 ULP 違うことがあり、コンパイルした版で再現できない）


@kernel
def get_flux_slau2(flux_rhs, specfic_heat_ratio, specific_heat_volum, var_a, var_b,
                   vecx, vecy, vecz, vect1x, vect1y, vect1z):
  """
  Advection flux on a face by the SLAU2 scheme

  var_b が面の左側（隣接セル）、var_a が右側（自セル）の原始変数。
  法線方向を u、接線方向を v として流束を組み、最後に直交座標系へ戻す。
  """

  # Primitive variables on face from left (neigboring) cell
  dens_l = var_b[0]
  uvel_l = var_b[1]*vecx   + var_b[2]*vecy   + var_b[3]*vecz
  vvel_l = var_b[1]*vect1x + var_b[2]*vect1y + var_b[3]*vect1z
  wvel_l = 0.0
  temp_l = var_b[4]
  pres_l = var_b[5]

  q2_l   = uvel_l*uvel_l + vvel_l*vvel_l + wvel_l*wvel_l
  sos_l  = thermodynamics.get_speedofsound(specfic_heat_ratio, dens_l, pres_l)
  enth_l = thermodynamics.get_enthalpy(specific_heat_volum, dens_l, temp_l, uvel_l, vvel_l, wvel_l, pres_l)

  # Primitive variables on face from right (self) cell
  dens_r = var_a[0]
  uvel_r = var_a[1]*vecx   + var_a[2]*vecy   + var_a[3]*vecz
  vvel_r = var_a[1]*vect1x + var_a[2]*vect1y + var_a[3]*vect1z
  wvel_r = 0.0
  temp_r = var_a[4]
  pres_r = var_a[5]

  q2_r   = uvel_r*uvel_r + vvel_r*vvel_r + wvel_r*wvel_r
  sos_r  = thermodynamics.get_speedofsound(specfic_heat_ratio, dens_r, pres_r)
  enth_r = thermodynamics.get_enthalpy(specific_heat_volum, dens_r, temp_r, uvel_r, vvel_r, wvel_r, pres_r)

  # Slau scheme and model parameters
  sos_m   = 0.50*(sos_l+sos_r)
  sos_inv = 1.0/sos_m
  mach_l  = uvel_l*sos_inv
  mach_r  = uvel_r*sos_inv
  # --SLAU, SLAU2
  mach_bar = min(1.0, np.sqrt(0.5*(q2_l + q2_r))*sos_inv)

  # Calculate u+-
  absu    =  ( dens_l*abs(uvel_l) + dens_r*abs(uvel_r) )/( dens_l + dens_r )
  gfact   = -max( min( mach_l, 0.0 ), -1.0 )*min( max( mach_r, 0.0 ), 1.0)
  chi_tmp =  1.0 - mach_bar
  chi     =  chi_tmp*chi_tmp
  u_p     =  uvel_l + (1.0-gfact)*absu + gfact*abs(uvel_l)
  u_m     =  uvel_r - (1.0-gfact)*absu - gfact*abs(uvel_r)

  # Mass flux
  ru_av = 0.50 * ( dens_l*u_p + dens_r*u_m - chi*(pres_r-pres_l)*sos_inv  )
  ru_l  = 0.50 * ( ru_av + abs(ru_av) )
  ru_r  = 0.50 * ( ru_av - abs(ru_av) )

  # Pressure flux
  # machf = 0 when abs(mach) < 1, machf = 1 when abs(mach) >= 1
  if abs(mach_l) < 1.0 :
    machf_l = 0.0
  else :
    machf_l = 1.0
  if abs(mach_r) < 1.0 :
    machf_r = 0.0
  else :
    machf_r = 1.0

  pmav     = 0.50*(pres_l + pres_r)
  alpha1   = 0.0
  machp_l  = mach_l + 1.0
  machm_r  = mach_r - 1.0
  machsq_l = mach_l*mach_l - 1.0
  machsq_r = mach_r*mach_r - 1.0
  presf_p  = (1.0 - machf_l)*( 0.25*(2.0-mach_l)*machp_l*machp_l + alpha1*mach_l*machsq_l*machsq_l ) \
           + machf_l*( 0.50*(1.0+np.sign( mach_l )) )
  presf_m  = (1.0 - machf_r)*( 0.25*(2.0+mach_r)*machm_r*machm_r - alpha1*mach_r*machsq_r*machsq_r ) \
           + machf_r*( 0.50*(1.0-np.sign( mach_r )) )

  # SLAU2（SLAU では最後の項が (1.0-chi)*(presf_p+presf_m-1.0)*pmav になる）
  p_av = pmav + 0.50*(presf_p-presf_m)*(pres_l-pres_r) \
       + np.sqrt(0.50*(q2_l + q2_r))*(presf_p+presf_m-1.0)*0.50*(dens_l+dens_r)*sos_m

  flux_tmp1 = ru_l        + ru_r
  flux_tmp2 = ru_l*uvel_l + ru_r*uvel_r + p_av
  flux_tmp3 = ru_l*vvel_l + ru_r*vvel_r
  # flux_tmp4 は 3 次元化したときの第 2 接線方向の成分。現状は使わない
  flux_tmp5 = ru_l*enth_l + ru_r*enth_r

  flux_rhs[0] = flux_tmp1
  flux_rhs[1] = flux_tmp2*vecx + flux_tmp3*vect1x
  flux_rhs[2] = flux_tmp2*vecy + flux_tmp3*vect1y
  flux_rhs[3] = flux_tmp2*vecz + flux_tmp3*vect1z
  flux_rhs[4] = flux_tmp5

  return flux_rhs


@kernel
def get_flux_haenel(flux_rhs, specfic_heat_ratio, specific_heat_volum, var_a, var_b,
                    vecx, vecy, vecz, vect1x, vect1y, vect1z):
  """
  Advection flux on a face by the Haenel scheme

  引数の意味は get_flux_slau2 と同じ。
  """

  # Primitive variables on face from left (neigboring) cell
  dens_l = var_b[0]
  uvel_l = var_b[1]*vecx   + var_b[2]*vecy   + var_b[3]*vecz
  vvel_l = var_b[1]*vect1x + var_b[2]*vect1y + var_b[3]*vect1z
  wvel_l = 0.0
  temp_l = var_b[4]
  pres_l = var_b[5]

  sos_l  = thermodynamics.get_speedofsound(specfic_heat_ratio, dens_l, pres_l)
  enth_l = thermodynamics.get_enthalpy(specific_heat_volum, dens_l, temp_l, uvel_l, vvel_l, wvel_l, pres_l)

  # Primitive variables on face from right (self) cell
  dens_r = var_a[0]
  uvel_r = var_a[1]*vecx   + var_a[2]*vecy   + var_a[3]*vecz
  vvel_r = var_a[1]*vect1x + var_a[2]*vect1y + var_a[3]*vect1z
  wvel_r = 0.0
  temp_r = var_a[4]
  pres_r = var_a[5]

  sos_r  = thermodynamics.get_speedofsound(specfic_heat_ratio, dens_r, pres_r)
  enth_r = thermodynamics.get_enthalpy(specific_heat_volum, dens_r, temp_r, uvel_r, vvel_r, wvel_r, pres_r)

  mach_l = uvel_l/sos_l
  mach_r = uvel_r/sos_r

  if abs(mach_l) <= 1.0 :
    u_tmp = uvel_l + sos_l
    m_tmp = mach_l + 1.0
    u_p = 0.25*( u_tmp*u_tmp )/sos_l
    p_p = 0.25*pres_l*( m_tmp*m_tmp )*(2.0-mach_l)
  else :
    u_p = 0.50*( uvel_l+abs(uvel_l) )
    p_p = 0.50*pres_l*(uvel_l+abs(uvel_l))/uvel_l

  if abs(mach_r) <= 1.0 :
    u_tmp = uvel_r - sos_r
    m_tmp = mach_r - 1.0
    u_m =-0.25*( u_tmp*u_tmp )/sos_r
    p_m = 0.25*pres_r*( m_tmp*m_tmp )*(2.0+mach_r)
  else :
    u_m = 0.50*( uvel_r-abs(uvel_r) )
    p_m = 0.50*pres_r*(uvel_r-abs(uvel_r))/uvel_r

  ru_av = dens_l*u_p + dens_r*u_m
  ru_l  = 0.50*( ru_av + abs(ru_av) )
  ru_r  = 0.50*( ru_av - abs(ru_av) )
  p_av  = p_p + p_m

  flux_tmp1 = ru_l        + ru_r
  flux_tmp2 = ru_l*uvel_l + ru_r*uvel_r + p_av
  flux_tmp3 = ru_l*vvel_l + ru_r*vvel_r
  flux_tmp5 = ru_l*enth_l + ru_r*enth_r

  flux_rhs[0] = flux_tmp1
  flux_rhs[1] = flux_tmp2*vecx + flux_tmp3*vect1x
  flux_rhs[2] = flux_tmp2*vecy + flux_tmp3*vect1y
  flux_rhs[3] = flux_tmp2*vecz + flux_tmp3*vect1z
  flux_rhs[4] = flux_tmp5

  return flux_rhs


@kernel
def accumulate_advection(num_face, num_face_bd, num_conserv, num_primitiv, kind_scheme,
                         face2cell, face2cell_bd, virtualcell_bd,
                         area_vec, area_vec_bd, length,
                         specfic_heat_ratio, specific_heat_volum, eps_muscl,
                         var_primitiv, var_primitiv_bd, var_gradient, var_limiter, var_rhs):
  """
  Accumulate the advection flux into the residual over all faces

  var_rhs は V*dQ/dt の符号を反転した量で、面の法線 area_vec[1:4] はセル a に対して
  内向きである。したがってセル a には引き、セル b には足す。

  面の左右の値は MUSCL でセル中心から外挿する（eps_muscl=0 なら 1 次精度）。
  """

  flux_rhs = np.zeros(num_conserv)
  var_a    = np.zeros(num_primitiv)
  var_b    = np.zeros(num_primitiv)

  # --Inner loop
  for n_face in range(0,num_face):
    area   = area_vec[0,n_face]
    vecx   = area_vec[1,n_face]
    vecy   = area_vec[2,n_face]
    vecz   = area_vec[3,n_face]
    vect1x = area_vec[4,n_face]
    vect1y = area_vec[5,n_face]
    vect1z = area_vec[6,n_face]

    dl_a   = length[0,n_face]
    dl_b   = length[1,n_face]

    n_cell_a = face2cell[0,n_face]
    n_cell_b = face2cell[1,n_face]

    for m in range(0,num_primitiv):
      grad_a = var_gradient[0,m,n_cell_a]*vecx + var_gradient[1,m,n_cell_a]*vecy + var_gradient[2,m,n_cell_a]*vecz
      grad_b = var_gradient[0,m,n_cell_b]*vecx + var_gradient[1,m,n_cell_b]*vecy + var_gradient[2,m,n_cell_b]*vecz
      var_a[m] = var_primitiv[m,n_cell_a] - eps_muscl*var_limiter[m,n_cell_a]*dl_a*grad_a
      var_b[m] = var_primitiv[m,n_cell_b] + eps_muscl*var_limiter[m,n_cell_b]*dl_b*grad_b

    if kind_scheme == KIND_SCHEME_HAENEL :
      get_flux_haenel(flux_rhs, specfic_heat_ratio, specific_heat_volum, var_a, var_b,
                      vecx, vecy, vecz, vect1x, vect1y, vect1z)
    else :
      get_flux_slau2(flux_rhs, specfic_heat_ratio, specific_heat_volum, var_a, var_b,
                     vecx, vecy, vecz, vect1x, vect1y, vect1z)

    for m in range(0,num_conserv):
      var_rhs[m,n_cell_a] = var_rhs[m,n_cell_a] - flux_rhs[m]*area
      var_rhs[m,n_cell_b] = var_rhs[m,n_cell_b] + flux_rhs[m]*area

  # --Boundary loop
  for n_face in range(0,num_face_bd):
    area   = area_vec_bd[0,n_face]
    vecx   = area_vec_bd[1,n_face]
    vecy   = area_vec_bd[2,n_face]
    vecz   = area_vec_bd[3,n_face]
    vect1x = area_vec_bd[4,n_face]
    vect1y = area_vec_bd[5,n_face]
    vect1z = area_vec_bd[6,n_face]

    n_cell_a = face2cell_bd[0,n_face]
    vcell    = virtualcell_bd[n_face]

    # 仮想セルがあるときはセルの値、無いときは境界値をそのまま使う
    for m in range(0,num_primitiv):
      var_a[m] = float(vcell)*var_primitiv[m,n_cell_a] + float(1-vcell)*var_primitiv_bd[m,n_face]
      var_b[m] = var_primitiv_bd[m,n_face]

    if kind_scheme == KIND_SCHEME_HAENEL :
      get_flux_haenel(flux_rhs, specfic_heat_ratio, specific_heat_volum, var_a, var_b,
                      vecx, vecy, vecz, vect1x, vect1y, vect1z)
    else :
      get_flux_slau2(flux_rhs, specfic_heat_ratio, specific_heat_volum, var_a, var_b,
                     vecx, vecy, vecz, vect1x, vect1y, vect1z)

    for m in range(0,num_conserv):
      var_rhs[m,n_cell_a] = var_rhs[m,n_cell_a] - flux_rhs[m]*area

  return var_rhs
