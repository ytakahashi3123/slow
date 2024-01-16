#!/usr/bin/env python3

# ***

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/07

import numpy as np
from orbital.orbital import orbital

class flowfield(orbital):


  def __init__(self):
    print("Calling class: flowfield")

    # gas_property_list
    self.id_specfic_heat_ratio  = 0
    self.id_gas_constant        = 1
    self.id_specific_heat_volum = 2
    self.id_specific_heat_press = 3
    # gas_property_dict
    self.key_specfic_heat_ratio  = 'specfic_heat_ratio'
    self.key_gas_constant        = 'gas_constant'
    self.key_specific_heat_volum = 'specific_heat_volume'
    self.key_specific_heat_press = 'specific_heat_press'

    # transport_coefficient_list
    self.id_viscosity    = 0
    self.id_thermal_cond = 1
    # transport_coefficient_dict
    self.key_viscosity             = 'viscosity'
    self.key_thermal_cond          = 'thermal_conductivity'
    self.key_viscosity_boundary    = 'viscosity_boundary'
    self.key_thermal_cond_boundary = 'thermal_conductivity_boundary'

    # var_conserv
    self.id_var_conserv_density     = 0
    self.id_var_conserv_momentumx   = 1
    self.id_var_conserv_momentumy   = 2
    self.id_var_conserv_momentumz   = 3
    self.id_var_conserv_totalenergy = 4

    # var_primitiv
    self.id_var_primitiv_density     = 0
    self.id_var_primitiv_velocityx   = 1
    self.id_var_primitiv_velocityy   = 2
    self.id_var_primitiv_velocityz   = 3
    self.id_var_primitiv_temperature = 4
    self.id_var_primitiv_pressure    = 5

    #


  def set_gas_properties(self, config):
    # Set gas properties 
    # Variables:
    # -- ID0:specfic_heat_ratio: Specfic heat ratio [-]
    # -- ID1: gas_constant: Gas constant calculated by R_u/m_w [J/K.kg]
    # ---- R_u: Universal gas constant [J/K.mol]
    # ---- m_w: Molecular weight [kg/mol]
    # -- ID2: specific_heat_volum: Specific heat at constant volume [J/K.kg]
    # -- ID3: specific_heat_volum: Specific heat at constant pressure [J/K.kg]

    print('Setting gas properties: specfic heat ratio and gas constant')

    specfic_heat_ratio = float( config['gas_properties']['gamma'] )

    gas_constant_univ = float( config['gas_properties']['gas_constant_univ'] )
    molecular_weight  = float( config['gas_properties']['molecular_weight'] )
    gas_constant = gas_constant_univ/molecular_weight

    specific_heat_volum = 2.5*gas_constant
    specific_heat_press = 3.5*gas_constant

    gas_property_dict = { 'specfic_heat_ratio':specfic_heat_ratio, \
                          'gas_constant':gas_constant, \
                          'specific_heat_volume':specific_heat_volum, \
                          'specific_heat_press':specific_heat_press }

    return gas_property_dict
    

  def define_variables(self, config, dimension_dict, geom_dict):
    
    print('Defining variables ')

    num_conserv      = dimension_dict['num_conservative']
    num_primitiv     = dimension_dict['num_primitive']

    num_cell         = geom_dict['num_cell']
    num_face_bd      = geom_dict['num_face_boundary']

    var_conserv      = np.zeros(num_conserv*num_cell).reshape(num_conserv,num_cell)
    var_primitiv     = np.zeros(num_primitiv*num_cell).reshape(num_primitiv,num_cell)
    var_primitiv_bd  = np.zeros(num_primitiv*num_face_bd).reshape(num_primitiv,num_face_bd)

    var_conserv_prev = np.zeros(2*num_conserv*num_cell).reshape(2,num_conserv,num_cell)

    viscosity        = np.zeros(num_cell).reshape(num_cell)
    thermal_cond     = np.zeros(num_cell).reshape(num_cell)
    viscosity_bd     = np.zeros(num_face_bd).reshape(num_face_bd)
    thermal_cond_bd  = np.zeros(num_face_bd).reshape(num_face_bd)

    transport_coefficient_dict =  {'viscosity':viscosity,                \
                                   'thermal_conductivity': thermal_cond, \
                                   'viscosity_boundary':viscosity_bd,    \
                                   'thermal_conductivity_boundary': thermal_cond_bd }

    return transport_coefficient_dict, var_primitiv, var_primitiv_bd, var_conserv, var_conserv_prev


  def initialize_flowfield(self, config, dimension_dict, geom_dict, metrics_dict, meshnode_dict, meshelem_dict, gas_property_dict, var_primitiv, var_conserv, var_conserv_prev):
    # Initialize flow field
    # Input: config, num_coord_list, coord_list
    # Output: conserv_list, primitiv_list
    # Input and output: -

    # Variables: 
    # -- Conservative variables
    # ---- ID0: rho: Density [kg/m3]
    # ---- ID1: rho*u: Momentum_x [kg/m2.s]
    # ---- ID2: rho*v: Momentum_y [kg/m2.s]
    # ---- ID3: rho*w: Momentum_z [kg/m2.s]
    # ---- ID4: E: Total energy [J/m3]
    # -- Primitive variables
    # ---- ID0: rho:  Density [kg/m3]
    # ---- ID1: u: Velocity_x [m/s]
    # ---- ID2: v: Velocity_y [m/s]
    # ---- ID3: w: Velocity_z [m/s]
    # ---- ID4: T: Temperature [K]
    # ---- ID5: p: Pressure [Pa]

    print('Setting initial conditions')

    flag_initial     = config['computational_setup']['flag_initial']
    kind_steady_mode = config['time_integration']['kind_steady_mode']

    num_conserv  = dimension_dict['num_conservative']
    num_cell     = geom_dict['num_cell']

    gas_constant        = gas_property_dict['gas_constant']
    specific_heat_volum = gas_property_dict['specific_heat_volume']

    if flag_initial :

      print('--from initial condition set in control file')

      iteration = 0
      #coord_cellcenter  = metrics_dict['coord_cellcenter']
      cell2tag          = geom_dict['cell2tag']
      tag2name_physics  = geom_dict['tag2name_physics']

      for i in range(0,num_cell):
        keyname_cell         = str( tag2name_physics[ cell2tag[i] ] )
        density_initial      = float( config['initial_settings'][keyname_cell]['density_initial'] )
        velocity_initial     =[float(ele) for ele in config['initial_settings'][keyname_cell]['velocity_initial']] 
        temperature_initial  = float( config['initial_settings'][keyname_cell]['temperature_initial'] )

        var_primitiv[0,i] = density_initial
        var_primitiv[1,i] = velocity_initial[0]
        var_primitiv[2,i] = velocity_initial[1]
        var_primitiv[3,i] = velocity_initial[2]
        var_primitiv[4,i] = temperature_initial
        var_primitiv[5,i] = self.get_pressure_eos(var_primitiv[0,i], gas_constant, var_primitiv[4,i])

        var_conserv[0,i]  = var_primitiv[0,i]
        var_conserv[1,i]  = var_primitiv[0,i]*var_primitiv[1,i]
        var_conserv[2,i]  = var_primitiv[0,i]*var_primitiv[2,i]
        var_conserv[3,i]  = var_primitiv[0,i]*var_primitiv[3,i]
        var_conserv[4,i]  = self.get_total_energy(var_primitiv[0,i], specific_heat_volum, var_primitiv[4,i], [var_primitiv[1,i], var_primitiv[2,i], var_primitiv[3,i] ])

      if kind_steady_mode == 'unsteady' :
        for n in range(0,num_conserv):
          var_conserv_prev[0,n,:] = var_conserv[n,:]
          var_conserv_prev[1,n,:] = var_conserv[n,:]

      # Output initial state
      if config['post_process']['flag_output_initial_state'] :
        self.routine_postprocess(config, iteration, meshnode_dict, meshelem_dict, metrics_dict, gas_property_dict, var_primitiv)

    else:

      print('--from restart file')

      # Reading restart data
      iteration, self.sum_rhs_init, var_conserv, var_conserv_prev = self.read_restart(config, dimension_dict, var_conserv, var_conserv_prev)

      # Primitive variables
      for n_cell in range(0,num_cell):
       conserv_tmp = var_conserv[:,n_cell]
       var_primitiv[:,n_cell] = self.get_primitive(gas_constant, specific_heat_volum, conserv_tmp)


    print('--Interation: ',iteration)

    return var_primitiv, var_conserv,  var_conserv_prev, iteration

  
  def set_transport_coefficients(self, config, geom_dict, gas_property_dict, var_primitiv, var_primitiv_bd, transport_coefficient_dict):
    # Transport _coefficients, viscosity and thermal conductivity, are calculated.
    # Variables:
    # -- ID:0. viscosity: Viscosity is evaluated based on Sutherland's law: mu=mu0*(T/T0)^(3/2) * (T0 + C)/(T + C)
    # -- ID:1. thermal_cond: Thermal conductivity is calculated by the viscosity, specific heat at constatn pressure, and Prandtl number: lambda=mu*Cp/Pr

    visc_mu_sutherland = float( config['transport_coefficients']['viscosity_sutherland']['mu0'] )
    visc_t_sutherland  = float( config['transport_coefficients']['viscosity_sutherland']['t0'] )
    visc_c_sutherland  = float( config['transport_coefficients']['viscosity_sutherland']['c'] )
    prandlt_number     = float( config['transport_coefficients']['prandlt_number'] )

    num_cell     = geom_dict['num_cell']
    num_face_bd  = geom_dict['num_face_boundary']

    specific_heat_press = gas_property_dict['specific_heat_press']

    temperature         = var_primitiv[4]
    temperature_bd      = var_primitiv_bd[4]

    viscosity           = transport_coefficient_dict['viscosity']
    thermal_cond        = transport_coefficient_dict['thermal_conductivity']
    viscosity_bd        = transport_coefficient_dict['viscosity_boundary']
    thermal_cond_bd     = transport_coefficient_dict['thermal_conductivity_boundary']

    for i in range(0,num_cell):
      viscosity[i]    = visc_mu_sutherland*(temperature[i]/visc_t_sutherland)**(1.50)*( ( visc_t_sutherland + visc_c_sutherland)/(temperature[i] + visc_c_sutherland) )
      thermal_cond[i] = viscosity[i]*specific_heat_press/prandlt_number

    for i in range(0,num_face_bd):
      viscosity_bd[i]    = visc_mu_sutherland*(temperature_bd[i]/visc_t_sutherland)**(1.50)*( ( visc_t_sutherland + visc_c_sutherland)/(temperature_bd[i] + visc_c_sutherland) )
      thermal_cond_bd[i] = viscosity_bd[i]*specific_heat_press/prandlt_number

    #transport_coefficient_list = [viscosity, thermal_cond]
    #transport_coefficient_list_bd = [viscosity_bd, thermal_cond_bd]

    return transport_coefficient_dict
