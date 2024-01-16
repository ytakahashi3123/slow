#!/usr/bin/env python3

# Program module to perform time integration

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/25

import numpy as np
from orbital.orbital import orbital
from time_integration import lusgs_diagonal
from time_integration import lusgs_sweep
from time_integration import update
from time_integration import explicit_euler


class time_integration(orbital):

  def __init__(self):

    print("Calling class: time_integration")

    return


  def initialize_time_integratioin(self, config, dimension_dict, geom_dict):

    print('Setting initial time integration variables')

    num_conserv = dimension_dict['num_conservative']
    num_cell    = geom_dict['num_cell']

    var_diagonal = np.zeros(num_cell).reshape(num_cell)
    var_dq = np.zeros(num_conserv*num_cell).reshape(num_conserv, num_cell)

    var_dt =  np.zeros(num_cell).reshape(num_cell)
    character_time = np.zeros(num_cell).reshape(num_cell)
    
    return var_diagonal, var_dq, var_dt, character_time


  def reinitialize_time_integratioin(self, config, var_diagonal, var_dq):

    var_diagonal[:] = 0.0
    var_dq[:,:] = 0.0

    return var_diagonal, var_dq


  def time_integration_routine(self, config, iteration_inner, dimension_dict, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_primitiv_bd, var_gradient, var_conserv, var_conserv_prev, var_rhs, var_dt, var_diagonal, var_dq):

    # 時間積分を行う: LU-SGS method

    kind_time_scheme = config['time_integration']['kind_time_scheme']
    flag_converged_inner = False

    if kind_time_scheme == 'implicit_lusgs':
    # Implicit scheme by LUSGS
    
      # Initialize variables
      var_diagonal, var_dq = self.reinitialize_time_integratioin(config, var_diagonal, var_dq)

      # Calculate diagonal
      var_diagonal, var_dq = lusgs_diagonal.get_diagonal(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_primitiv_bd, var_conserv, var_conserv_prev, var_rhs, var_dt, var_diagonal, var_dq)

      # Sweep
      var_dq = lusgs_sweep.sweep_jacobian(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_conserv, var_diagonal, var_dq)

      # Update
      var_conserv = update.update_solution(config, geom_dict, var_conserv, var_dq)

      # Display delta conservative variables
      sum_dq = self.display_deltaq(config, iteration_inner, dimension_dict, var_dq)

      # Check convergence 
      flag_converged_inner = self.check_convergence_inner(config, flag_converged_inner, sum_dq)

    elif kind_time_scheme == 'explicit_euler' :
      var_conserv = explicit_euler.explicit_euler(config, geom_dict, metrics_dict, var_dt, var_rhs, var_conserv)

    else:
      print('Error in kind_time_shceme of control file.', kind_time_shceme)
      print('Program stopped')
      exit()

    # Primitive variables
    var_primitiv = self.update_primitive(config, geom_dict, metrics_dict, gas_property_dict, var_conserv, var_primitiv)

    return var_conserv, var_primitiv, flag_converged_inner


  def set_conservative_previous(self, config, var_conserv, var_conserv_prev):
    # Set privious conservative variables

    var_conserv_prev[1,:,:] = var_conserv_prev[0,:,:]
    var_conserv_prev[0,:,:] = var_conserv[:,:]

    return var_conserv_prev


  def update_primitive(self, config, geom_dict, metrics_dict, gas_property_dict, var_conserv, var_primitiv):
    # Updating primitive variables from conservative variables:
    # --Densiry: m=rho
    # --Momentum: mu,mv,mw = rho*u, rho*v, rho*w
    # --Total energy: E =rho*Cv*T + 0.5*rho*U^2
    
    # Input parameters
    num_cell            = geom_dict['num_cell']
    
    gas_constant        = gas_property_dict['gas_constant']
    specific_heat_volum = gas_property_dict['specific_heat_volume']

    coord_cellcenter = metrics_dict['coord_cellcenter']

    # Update conservative variables
    for n_cell in range(0,num_cell):
      conserv_tmp = var_conserv[:,n_cell]
      var_primitiv[:,n_cell] = self.get_primitive(gas_constant, specific_heat_volum, conserv_tmp)
    
    # Check variables
    for n_cell in range(0,num_cell):
      if var_primitiv[0,n_cell] < 0.0 or var_primitiv[4,n_cell] < 0.0 or var_primitiv[5,n_cell] < 0.0 :
        print( n_cell, coord_cellcenter[0,n_cell], coord_cellcenter[1,n_cell],  var_primitiv[0,n_cell], var_primitiv[4,n_cell], var_primitiv[5,n_cell] )
        exit()
    #flag_fail = any((x < 0 for x in var_primitiv[4,:]))

    return var_primitiv


  def set_timestep(self, config, geom_dict, character_time, var_dt):
    
    # set time step at each cell

    kind_time_determine = config['time_integration']['kind_time_determine']
    kind_time_stepping  = config['time_integration']['kind_time_stepping']
    courant_number      = config['time_integration']['courant_number']
    timestep_constant   = config['time_integration']['timestep_constant'] 

    num_cell = geom_dict['num_cell']

    if kind_time_determine == 'cfl' :
      # Time step is determined by Courant number
      if kind_time_stepping == 'local' :
        # Local time stepping
        for n_cell in range(0,num_cell):
          var_dt[n_cell] = courant_number*( character_time[n_cell] + 1.e-20)
      elif kind_time_stepping == 'global' :
        # Global time stepping
        dt_global_tmp = 1.e+20
        for n_cell in range(0,num_cell):
          dt_local_tmp  = courant_number*( character_time[n_cell] + 1.e-20)
          dt_global_tmp = min(dt_global_tmp, dt_local_tmp)
        for n_cell in range(0,num_cell):
          var_dt[n_cell] = dt_global_tmp
      else:
        print('Error in kind_time_stepping of control file: ', kind_time_stepping )
        print('Program stopped')
        exit()

    elif kind_time_determine == 'dt' :
      # Time step is determined by time step given
        for n_cell in range(0,num_cell):
          var_dt[n_cell] = timestep_constant

    else :
      print('Error in kind_time_fix of control file: ', kind_time_determine )
      print('Program stopped')
      exit()

    # Display
    print('Maximum time step:', np.max(var_dt), 'Minimum time step:', np.min(var_dt),)

    return var_dt


  def get_characteristic_time(self, config, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_primitiv_bd, character_time):

    # Local time steppingnにおけるTime stepを計算するために各セルでのcharacteristic_timeを取得する

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

    viscosity           = transport_coefficient_dict['viscosity']
    viscosity_bd        = transport_coefficient_dict['viscosity_boundary']

    # Initialization
    character_time[:] = 0.0

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
      prim_a = var_primitiv[:,n_cell_a]
      prim_b = var_primitiv[:,n_cell_b]

      # Values on cell interface
      prim   = 0.50*( prim_a + prim_b )
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
      character_time[n_cell_a] = max(character_time[n_cell_a], eigenvalue*area)
      character_time[n_cell_b] = max(character_time[n_cell_b], eigenvalue*area)


    for n_face in range(0,num_face_bd):
      # Face area vector
      area = area_vec_bd[0,n_face]
      vecx = area_vec_bd[1,n_face]
      vecy = area_vec_bd[2,n_face]
      vecz = 0.0

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
      character_time[n_cell_a] = max(character_time[n_cell_a], eigenvalue*area)

    # characteristic time
    for n_cell in range(0,num_cell):
      character_time[n_cell] = volume[n_cell]/character_time[n_cell]


    return character_time

