#!/usr/bin/env python3

# Program to sweep Jacobian matrix for LU-SGS routine

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/31

import numpy as np
from orbital.orbital import orbital

def sweep_jacobian(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_conserv, var_diagonal, var_dq):

  # Jacobian matrix
  def jacobian_routine():

    # Set derivatives
    # -- Energy derivatives
    pett =  gamm - 1.0
    # -- Total density derivatives
    prho = -pett * q


    # Jacobian matrix (Left side)
    jacobian[0,0] = 0.0 - eigenvalue
    jacobian[0,1] = vecx
    jacobian[0,2] = vecy
    jacobian[0,3] = vecz
    jacobian[0,4] = 0.0

    jacobian[1,0] =  vecx*prho            - cvel*uvel
    jacobian[1,1] = -vecx*(pett-1.0)*uvel + cvel - eigenvalue
    jacobian[1,2] = -vecx*pett*vvel       + vecy*uvel
    jacobian[1,3] = -vecx*pett*wvel       + vecz*uvel
    jacobian[1,4] =  vecx*pett

    jacobian[2,0] =  vecy*prho            - cvel*vvel
    jacobian[2,1] = -vecy*pett*uvel       + vecx*vvel
    jacobian[2,2] = -vecy*(pett-1.0)*vvel + cvel - eigenvalue
    jacobian[2,3] = -vecy*pett*wvel       + vecz*vvel
    jacobian[2,4] =  vecy*pett

    jacobian[3,0] =  vecz*prho            - cvel*wvel
    jacobian[3,1] = -vecz*pett*uvel       + vecx*wvel
    jacobian[3,2] = -vecz*pett*vvel       + vecy*wvel
    jacobian[3,3] = -vecz*(pett-1.0)*wvel + cvel - eigenvalue
    jacobian[3,4] =  vecz*pett

    jacobian[4,0] =  cvel*(prho-enth)
    jacobian[4,1] = -cvel*pett*uvel + vecx*enth
    jacobian[4,2] = -cvel*pett*vvel + vecy*enth
    jacobian[4,3] = -cvel*pett*wvel + vecz*enth
    jacobian[4,4] =  cvel*(pett+1.0) - eigenvalue

    return


  # Main routine
  
  # Input parameters
  num_conserv  = dimension_dict['num_conservative']

  num_face    = geom_dict['num_face_inner']
  num_face_bd = geom_dict['num_face_boundary']
  num_cell    = geom_dict['num_cell']
  face2cell   = geom_dict['face2cell_inner']

  area_vec    = metrics_dict['area_vec_inner']
  area_vec_bd = metrics_dict['area_vec_boundary']
  length      = metrics_dict['length_inner']
  length_bd   = metrics_dict['length_boundary']
  volume      = metrics_dict['volume_cell']

  specfic_heat_ratio  = gas_property_dict['specfic_heat_ratio']
  specific_heat_volum = gas_property_dict['specific_heat_volume']

  viscosity           = transport_coefficient_dict['viscosity']


  # Initialize
  jacobian      = np.zeros(num_conserv*num_conserv).reshape(num_conserv,num_conserv)
  #kroneko_delta = np.identity(num_conserv, dtype=float) 
  eigenvalue = 0.0

  # Forward sweep
  for n_face in range(0,num_face):

    # Face area vector
    area = area_vec[0,n_face]
    vecx = area_vec[1,n_face]
    vecy = area_vec[2,n_face]
    vecz = 0.0
    # Length
    leng_a = length[0,n_face]
    leng_b = length[1,n_face]


    # Cell ID
    # --from self cell side
    n_cell_a = face2cell[0,n_face]
    # --from neigboring cell
    n_cell_b = face2cell[1,n_face]

    # Primitive variables
    prim = var_primitiv[:,n_cell_a]
    dens   = prim[0]
    uvel   = prim[1]
    vvel   = prim[2]
    wvel   = prim[3]
    temp   = prim[4]
    pres   = prim[5]

    # Contravariant velocity
    cvel = uvel*vecx + vvel*vecy + wvel*vecz

    # Thermodynamic properties: Enthalpy, specific ratio, speed of sound
    enth  = orbital.get_enthalpy("self", specific_heat_volum, dens, temp, [uvel,vvel,wvel], pres)
    sos   = orbital.get_speedofsound("self", specfic_heat_ratio, dens, pres)
    gamm  = specfic_heat_ratio

    # Dynamic pressure
    q   = 0.50*(uvel*uvel + vvel*vvel + wvel*wvel)


    # Maximum eigenvalue of Jacobian matrix: lambda=|u|+c + (viscous contribution)
    # (viscous contribution)=2*mu/(rho*d)
    lenght_tmp = leng_a + leng_b
    visc_tmp   = viscosity[n_cell_a]
    eigenvalue = abs(cvel) + sos + 2.0*visc_tmp/(dens*lenght_tmp)


    # Jacobian matrix
    jacobian_routine()


    # Sweep
    var_inv = 0.50*area/var_diagonal[n_cell_b]
    for m in range(0,num_conserv):
      dq_tmp = np.dot( jacobian[m,:],var_dq[:,n_cell_a])
      var_dq[m,n_cell_b] = var_dq[m,n_cell_b] - dq_tmp*var_inv


  # Backward sweep
  for n_face_b in range(0,num_face):
    n_face = num_face-n_face_b-1

    # Face area vector
    area = area_vec[0,n_face]
    vecx =-area_vec[1,n_face]
    vecy =-area_vec[2,n_face]
    vecz = 0.0
    # Length
    leng_a = length[0,n_face]
    leng_b = length[1,n_face]


    # Cell ID
    # --from self cell side
    n_cell_a = face2cell[0,n_face]
    # --from neigboring cell
    n_cell_b = face2cell[1,n_face]

    # Primitive variables
    prim = var_primitiv[:,n_cell_b]
    # Primitive variables
    dens   = prim[0]
    uvel   = prim[1]
    vvel   = prim[2]
    wvel   = prim[3]
    temp   = prim[4]
    pres   = prim[5]

    # Contravariant velocity
    cvel = uvel*vecx + vvel*vecy + wvel*vecz

    # Thermodynamic properties: Enthalpy, specific ratio, speed of sound
    enth  = orbital.get_enthalpy("self", specific_heat_volum, dens, temp, [uvel,vvel,wvel], pres)
    sos   = orbital.get_speedofsound("self", specfic_heat_ratio, dens, pres)
    gamm  = specfic_heat_ratio

    # Dynamic pressure
    q   = 0.50*(uvel*uvel + vvel*vvel + wvel*wvel)


    # Maximum eigenvalue of Jacobian matrix: lambda=|u|+c + (viscous contribution)
    # (viscous contribution)=2*mu/(rho*d)
    lenght_tmp = leng_a + leng_b
    visc_tmp   = viscosity[n_cell_b]
    eigenvalue = abs(cvel) + sos + 2.0*visc_tmp/(dens*lenght_tmp)


    # Jacobian matrix
    jacobian_routine()


    # Sweep
    var_inv = 0.50*area/var_diagonal[n_cell_a]
    for m in range(0,num_conserv):
      dq_tmp = np.dot( jacobian[m,:],var_dq[:,n_cell_b])
      var_dq[m,n_cell_a] = var_dq[m,n_cell_a] - dq_tmp*var_inv

  return var_dq

