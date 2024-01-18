#!/usr/bin/env python3

# Program to sweep Jacobian matrix for LU-SGS routine

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/31

import numpy as np
from orbital.orbital import orbital

@orbital.time_measurement_decorated
def get_diagonal(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_primitiv_bd, var_conserv, var_conserv_prev, var_rhs, var_dt, var_diagonal, var_dq):

  # Main routine
  
  # Input parameters
  kind_steady_mode         = config['time_integration']['kind_steady_mode']  # 'steady': Steady flow computation, or 'unsteady': Unsteady flow computation
  lusgs_beta               = config['time_integration']['lusgs_beta']
  kind_backward_difference = config['time_integration']['kind_backward_difference']
  var_dt_const             = config['time_integration']['timestep_outer'] 

  #num_conserv  = dimension_dict['num_conservative']
  #num_primitiv = dimension_dict['num_primitive']

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

    # Values on cell interface
    prim = 0.50*( prim_a + prim_b )
    dens   = prim[0]
    uvel   = prim[1]
    vvel   = prim[2]
    wvel   = prim[3]
    #temp   = prim[4]
    pres   = prim[5]

    # Contravariant velocity
    cvel = uvel*vecx + vvel*vecy + wvel*vecz

    # Thermodynamic properties: speed of sound
    sos   = orbital.get_speedofsound("self", specfic_heat_ratio, dens, pres)


    # Maximum eigenvalue of Jacobian matrix: lambda=|u|+c + (viscous contribution) (because lambda=|u|+c, |u|-c, |u| in inviscid flow)
    # (viscous contribution)=2*mu/(rho*d)
    lenght_tmp = leng_a + leng_b
    visc_tmp   = 0.50*( viscosity[n_cell_a] + viscosity[n_cell_b] )
    eigenvalue = abs(cvel) + sos + 2.0*visc_tmp/(dens*lenght_tmp)

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
    dens   = prim[0]
    uvel   = prim[1]
    vvel   = prim[2]
    wvel   = prim[3]
    #temp   = prim[4]
    pres   = prim[5]

    # Contravariant velocity
    cvel = uvel*vecx + vvel*vecy + wvel*vecz

    # Thermodynamic properties: speed of sound
    sos   = orbital.get_speedofsound("self", specfic_heat_ratio, dens, pres)


    # Maximum eigenvalue of Jacobian matrix: lambda=|u|+c + (viscous contribution) (because lambda=|u|+c, |u|-c, |u| in inviscid flow)
    # (viscous contribution)=2*mu/(rho*d)
    lenght_tmp = leng_a + leng_b
    visc_tmp   = float(vcell_bd)*0.50*( viscosity[n_cell_a] + viscosity_bd[n_face] ) + float(1-vcell_bd)*viscosity_bd[n_face]
    eigenvalue = abs(cvel) + sos + 2.0*visc_tmp/(dens*lenght_tmp)

    # Set diagonal values
    var_diagonal[n_cell_a] = var_diagonal[n_cell_a] + eigenvalue*area


  #var_dq = var_rhs


  # Diagonal element for factoriization matrix and delta Q
  if kind_steady_mode == 'steady' :
    # Steady flow
    for n_cell in range(0,num_cell):
      var_diagonal[n_cell] = volume[n_cell]/var_dt[n_cell] + 0.50*lusgs_beta*var_diagonal[n_cell]
      var_dq[:,n_cell]     = -var_rhs[:,n_cell]/var_diagonal[n_cell]
  elif kind_steady_mode == 'unsteady' :
    # Unsteady flow with Pseudo-time stepping
    if kind_backward_difference == '2nd_backward_diff' :
      # - 2nd order accuracy backward difference
      # - (Volume*(3/2dt+1/d_tau) + 0.5*eigenvalue)
      for n_cell in range(0,num_cell):
        var_diagonal[n_cell] = 1.50*volume[n_cell]/var_dt_const + volume[n_cell]/var_dt[n_cell] + 0.50*var_diagonal[n_cell]
        dq_unst              = ( 1.50*var_conserv[:,n_cell] - 2.0*var_conserv_prev[0,:,n_cell] +0.50*var_conserv_prev[1,:,n_cell] )*volume[n_cell]/var_dt_const
        var_dq[:,n_cell]     = ( -var_rhs[:,n_cell]-dq_unst )/var_diagonal[n_cell]
    elif kind_backward_difference == '1st_backward_diff' :
      # - 1st order accuracy backward difference
      # - (Volume*(1/dt+1/d_tau) + 0.5*eigenvalue)
      for n_cell in range(0,num_cell):
        var_diagonal[n_cell] = volume[n_cell]/var_dt_const + volume[n_cell]/var_dt[n_cell] + 0.50*var_diagonal[n_cell]
        dq_unst              = ( var_conserv[:,n_cell] - var_conserv_prev[0,:,n_cell] )*volume[n_cell]/var_dt_const
        var_dq[:,n_cell]     = ( -var_rhs[:,n_cell]-dq_unst )/var_diagonal[n_cell]
    else :
      print('Error in kind_backward_difference of control file: ', kind_backward_difference )
      print('Program stopped')
      exit()
  else :
    print('Error in kind_steady_mode of control file: ', kind_steady_mode )
    print('Program stopped')
    exit()


  return var_diagonal, var_dq

