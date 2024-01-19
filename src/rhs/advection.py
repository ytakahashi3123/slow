#!/usr/bin/env python3

# Program to calculate advection flux on interface

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/25

import numpy as np
from orbital.orbital import orbital

@orbital.time_measurement_decorated
#@orbital.parallel_execution_decorated(max_workers=4)
def flux_advection(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, var_primitiv, var_primitiv_bd, var_gradient, var_limiter, var_rhs):


  # Advection scheme: SLAU/SLAU2
  def slau2():

    # Primitive variables on face from left (neigboring) cell
    dens_l = var_b[0]
    uvel_l = var_b[1]*vecx   + var_b[2]*vecy   + var_b[3]*vecz
    vvel_l = var_b[1]*vect1x + var_b[2]*vect1y + var_b[3]*vect1z
    wvel_l = 0.0
    temp_l = var_b[4]
    pres_l = var_b[5]

#    q2_l   = var_b[1]**2 + var_b[2]**2 + var_b[3]**2
    q2_l   = uvel_l**2 + vvel_l**2 + wvel_l**2
    sos_l  = orbital.get_speedofsound('self',specfic_heat_ratio, dens_l, pres_l)
    enth_l = orbital.get_enthalpy('self',specific_heat_volum, dens_l, temp_l, [uvel_l,vvel_l,wvel_l], pres_l)

    # Primitive variables on face from right (self) cell
    dens_r = var_a[0]
    uvel_r = var_a[1]*vecx   + var_a[2]*vecy   + var_a[3]*vecz
    vvel_r = var_a[1]*vect1x + var_a[2]*vect1y + var_a[3]*vect1z
    wvel_r = 0.0
    temp_r = var_a[4] 
    pres_r = var_a[5]

#    q2_r   = var_a[1]**2 + var_a[2]**2 + var_a[3]**2
    q2_r   = uvel_r**2 + vvel_r**2 + wvel_r**2
    sos_r  = orbital.get_speedofsound('self',specfic_heat_ratio, dens_r, pres_r)
    enth_r = orbital.get_enthalpy('self',specific_heat_volum, dens_r, temp_r, [uvel_r,vvel_r,wvel_r], pres_r)

    # Slau scheme and model parameters
    sos_m   = 0.50*(sos_l+sos_r)
    sos_inv = 1.0/sos_m
    mach_l = uvel_l*sos_inv
    mach_r = uvel_r*sos_inv
    # --SLAU, SLAU2
    mach_bar = min(1.0, np.sqrt(0.5*(q2_l + q2_r))*sos_inv)

    # Calculate u+-
    absu  =  ( dens_l*abs(uvel_l) + dens_r*abs(uvel_r) )/( dens_l + dens_r )
    gfact = -max( min( mach_l, 0.0 ), -1.0 )*min( max( mach_r, 0.0 ), 1.0)
    chi   =  (1.0-mach_bar)**2
    u_p   =  uvel_l + (1.0-gfact)*absu + gfact*abs(uvel_l)
    u_m   =  uvel_r - (1.0-gfact)*absu - gfact*abs(uvel_r)

    # Mass flux
    ru_av = 0.50 * ( dens_l*u_p + dens_r*u_m - chi*(pres_r-pres_l)*sos_inv  )
    ru_l  = 0.50 * ( ru_av + abs(ru_av) )
    ru_r  = 0.50 * ( ru_av - abs(ru_av) )

    # Pressure flux
    #machf_l  = 0.50 + 0.50*np.sign( abs(mach_l)-1.0 ) # machf_l = 0 when abs(mach_l) < 1, machf_l = 1 when abs(mach_l) >= 1
    #machf_r  = 0.50 + 0.50*np.sign( abs(mach_r)-1.0 ) # machf_r = 0 when abs(mach_r) < 1, machf_r = 1 when abs(mach_r) >= 1
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
    presf_p  = (1.0 - machf_l)*( 0.25*(2.0-mach_l)*(mach_l+1.0)**2 + alpha1*mach_l*(mach_l**2-1.0)**2 ) + machf_l*( 0.50*(1.0+np.sign( mach_l )) )
    presf_m  = (1.0 - machf_r)*( 0.25*(2.0+mach_r)*(mach_r-1.0)**2 - alpha1*mach_r*(mach_r**2-1.0)**2 ) + machf_r*( 0.50*(1.0-np.sign( mach_r )) )

    # SLAU
    #
    #p_av  = pmav + 0.50*(presf_p-presf_m)*(pres_l-pres_r) + (1.0-chi)*(presf_p+presf_m-1.0)*pmav
    # SLAU2
    p_av  = pmav + 0.50*(presf_p-presf_m)*(pres_l-pres_r) + np.sqrt(0.50*(q2_l + q2_r))*(presf_p+presf_m-1.0)*0.50*(dens_l+dens_r)*sos_m

    flux_tmp1 = ru_l        + ru_r
    flux_tmp2 = ru_l*uvel_l + ru_r*uvel_r + p_av
    flux_tmp3 = ru_l*vvel_l + ru_r*vvel_r 
    flux_tmp4 = ru_l*wvel_l + ru_r*wvel_r 
    flux_tmp5 = ru_l*enth_l + ru_r*enth_r

    flux_rhs[0] = flux_tmp1
    flux_rhs[1] = flux_tmp2*vecx + flux_tmp3*vect1x
    flux_rhs[2] = flux_tmp2*vecy + flux_tmp3*vect1y
    flux_rhs[3] = flux_tmp2*vecz + flux_tmp3*vect1z
    flux_rhs[4] = flux_tmp5

    return


  # Advection scheme: Haenel
  def haenel():

    # Primitive variables on face from left (neigboring) cell
    dens_l = var_b[0]
    uvel_l = var_b[1]*vecx   + var_b[2]*vecy   + var_b[3]*vecz
    vvel_l = var_b[1]*vect1x + var_b[2]*vect1y + var_b[3]*vect1z
    wvel_l = 0.0
    temp_l = var_b[4]
    pres_l = var_b[5]

    sos_l  = orbital.get_speedofsound('self',specfic_heat_ratio, dens_l, pres_l)
    enth_l = orbital.get_enthalpy('self',specific_heat_volum, dens_l, temp_l, [uvel_l,vvel_l,wvel_l], pres_l)

    # Primitive variables on face from right (self) cell
    dens_r = var_a[0]
    uvel_r = var_a[1]*vecx   + var_a[2]*vecy   + var_a[3]*vecz
    vvel_r = var_a[1]*vect1x + var_a[2]*vect1y + var_a[3]*vect1z
    wvel_r = 0.0
    temp_r = var_a[4] 
    pres_r = var_a[5]

    sos_r  = orbital.get_speedofsound('self',specfic_heat_ratio, dens_r, pres_r)
    enth_r = orbital.get_enthalpy('self',specific_heat_volum, dens_r, temp_r, [uvel_r,vvel_r,wvel_r], pres_r)

    # Slau scheme and model parameters
    mach_l = uvel_l/sos_l
    mach_r = uvel_r/sos_r

    if abs(mach_l) <= 1.0 :
      u_p = 0.25*( (uvel_l+sos_l)**2 )/sos_l
      p_p = 0.25*pres_l*((mach_l+1.0)**2)*(2.0-mach_l)
    else :
      u_p = 0.50*( uvel_l+abs(uvel_l) )
      p_p = 0.50*pres_l*(uvel_l+abs(uvel_l))/uvel_l

    if abs(mach_r) <= 1.0 :
      u_m =-0.25*( (uvel_r-sos_r)**2 )/sos_r
      p_m = 0.25*pres_r*((mach_r-1.0)**2)*(2.0+mach_r)
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
    flux_tmp4 = ru_l*wvel_l + ru_r*wvel_r 
    flux_tmp5 = ru_l*enth_l + ru_r*enth_r

    flux_rhs[0] = flux_tmp1
    flux_rhs[1] = flux_tmp2*vecx + flux_tmp3*vect1x
    flux_rhs[2] = flux_tmp2*vecy + flux_tmp3*vect1y
    flux_rhs[3] = flux_tmp2*vecz + flux_tmp3*vect1z
    flux_rhs[4] = flux_tmp5

    return


  # Main routine
  eps_muscl   = config['gradient_setting']['eps_muscl']
  kind_scheme = config['numericalflux_setting']['advection_scheme']
   

  # Input parameters
  num_conserv  = dimension_dict['num_conservative']

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
  specific_heat_volum = gas_property_dict['specific_heat_volume']


  # Initialize
  #flux_rhs = np.zeros(num_conserv).reshape(num_conserv)
  flux_rhs = np.linspace(0, 0, num_conserv, dtype=float)


  # Flux calculation
  # --Inner loop
  for n_face in range(0,num_face):
    # Face area vector
    area   = area_vec[0,n_face]
    vecx   = area_vec[1,n_face]
    vecy   = area_vec[2,n_face]
    vecz   = area_vec[3,n_face]
    vect1x = area_vec[4,n_face]
    vect1y = area_vec[5,n_face]
    vect1z = area_vec[6,n_face]

    dl_a   = length[0,n_face]
    dl_b   = length[1,n_face]

    # Calculate values on face, interpolated from left (self) and right sides (neigboring) cell
    # --from self cell side
    n_cell_a = face2cell[0,n_face]
    var_a  = var_primitiv[:,n_cell_a] \
           - eps_muscl*var_limiter[:,n_cell_a]*dl_a*( var_gradient[0,:,n_cell_a]*vecx + var_gradient[1,:,n_cell_a]*vecy + var_gradient[2,:,n_cell_a]*vecz)
    # --from neigboring cell
    n_cell_b = face2cell[1,n_face]
    var_b  = var_primitiv[:,n_cell_b] \
           + eps_muscl*var_limiter[:,n_cell_b]*dl_b*( var_gradient[0,:,n_cell_b]*vecx + var_gradient[1,:,n_cell_b]*vecy + var_gradient[2,:,n_cell_b]*vecz)

    # Call advection scheme routine
    if kind_scheme == 'slau2' :
      slau2()
    elif kind_scheme == 'haenel' :
      haenel()
    else :
      slau2()

    # Calculate fluxes
    var_rhs[:,n_cell_a] = var_rhs[:,n_cell_a] - flux_rhs[:]*area
    var_rhs[:,n_cell_b] = var_rhs[:,n_cell_b] + flux_rhs[:]*area


  # --Boundary loop
  for n_face in range(0,num_face_bd):
    # Face area vector
    area   = area_vec_bd[0,n_face]
    vecx   = area_vec_bd[1,n_face]
    vecy   = area_vec_bd[2,n_face]
    vecz   = area_vec_bd[3,n_face]
    vect1x = area_vec_bd[4,n_face]
    vect1y = area_vec_bd[5,n_face]
    vect1z = area_vec_bd[6,n_face]

    # Calculate values on face, interpolated from left (self) and right sides (neigboring) cell
    # --from self cell side
    n_cell_a = face2cell_bd[0,n_face]
    vcell    = virtualcell_bd[n_face]
    var_a    = float(vcell)*var_primitiv[:,n_cell_a] + float(1-vcell)*var_primitiv_bd[:,n_face]
    # --from neigboring cell
    var_b    = var_primitiv_bd[:,n_face]

    # Call advection scheme routine
    if kind_scheme == 'slau2' :
      slau2()
    elif kind_scheme == 'haenel' :
      haenel()
    else :
      slau2()

    # Calculate fluxes
    var_rhs[:,n_cell_a] = var_rhs[:,n_cell_a] - flux_rhs[:]*area


  return var_rhs

