# coding:utf-8
#!/usr/bin/env python3

# SLOW: Supersonic flow code with low-learning cost
# Version 0.05

# Author: Yusuke Takahashi, Hokkaido University
# Date: 2024/01/31


import numpy as np
from orbital.orbital import orbital
from meshdata.meshdata import meshdata
from flowfield.flowfield import flowfield
from boundary.boundary import boundary 
from gradient.gradient import gradient
from rhs.rhs import rhs
from time_integration.time_integration import time_integration


def main():


  # 設定ファイルの読み込み
  file_control_default = orbital.file_control_default
  arg                  = orbital.argument(file_control_default)
  file_control         = arg.file
  config               = orbital.read_config_yaml(file_control)


  # Set dimension parameters
  dimension_dict = orbital.set_dimension()


  # Set orbital parameters
  orbital.set_orbital_parameters(dimension_dict)


  # Set geometry variables
  # --Boundary data
  #bd_list = geometry.read_boundary_data(config)
  # --Geometry data
  #grid_list, coord_node_list = geometry.read_geometry_data(config)
  # --Set boundary attribution and index
  #grid_list = geometry.set_boundary_attribute(config, bd_list, grid_list)
  # --Set face and cell variables
  #geom_list = geometry.set_geometry_face(config, grid_list, coord_node_list)
  # --Metrics
  #metrics_list = geometry.set_metrics(config, coord_node_list, geom_list)
  #
  meshnode_dict, meshelem_dict, geom_dict, metrics_dict = meshdata.set_mesh_routine(config)
  cellcenter = metrics_dict['coord_cellcenter']

  # Make directories for output data
  orbital.make_directory_output(config)


  # Set gas properties
  gas_property_dict = flowfield.set_gas_properties(config)


  # Define flow field variables (inner-domain/boundary)
  transport_coefficient_dict, \
  var_primitiv, var_primitiv_bd, \
  var_conserv, \
  var_conserv_prev = flowfield.define_variables(config, dimension_dict, geom_dict)

  # Set initial condtions
  var_primitiv, \
  var_conserv, \
  var_conserv_prev, \
  iteration = flowfield.initialize_flowfield(config,            \
                                             dimension_dict,    \
                                             geom_dict,         \
                                             metrics_dict,      \
                                             meshnode_dict,     \
                                             meshelem_dict,     \
                                             gas_property_dict, \
                                             var_primitiv,      \
                                             var_conserv,       \
                                             var_conserv_prev)


  # Initialize gradient variables
  var_gradient, var_limiter, var_neig_maxmin = gradient.initialize_gradient(config, dimension_dict, geom_dict)

  # Initialize RHS variables
  var_rhs = rhs.initialize_rhs(config, dimension_dict, geom_dict)

  # Initialize LUSGS variables
  var_diagonal, var_dq, var_dt, character_time = time_integration.initialize_time_integratioin(config, dimension_dict, geom_dict)


  # Computational parameters
  iteration_maximum            = config['computational_setup']['iteration_maximum']
  frequency_output_postprocess = config['post_process']['frequency_output']
  frequency_output_restart     = config['restart_process']['frequency_output']

  kind_steady_mode        = config['time_integration']['kind_steady_mode']
  kind_time_scheme        = config['time_integration']['kind_time_scheme']
  iteration_inner_maximum = config['time_integration']['iteration_inner_maximum']
  if kind_steady_mode == 'steady': 
    iteration_inner_maximum = 1
  if kind_time_scheme == 'explicit_euler' :
    iteration_inner_maximum = 1

  # Time marching
  while iteration < iteration_maximum :

    flag_converged_outer = False

    for iteration_inner in range(0,iteration_inner_maximum ):
    
      # Set boundary condition 
      var_primitiv_bd = boundary.boundary_condition(config, geom_dict, metrics_dict, \
                                                    gas_property_dict, var_primitiv, var_conserv, \
                                                    var_primitiv_bd)

      # Set transport coefficients
      transport_coefficient_dict = flowfield.set_transport_coefficients(config, geom_dict, \
                                                                        gas_property_dict, var_primitiv, var_primitiv_bd, \
                                                                        transport_coefficient_dict)

      # Time step
      character_time = time_integration.get_characteristic_time(config, geom_dict, metrics_dict, \
                                                                gas_property_dict, transport_coefficient_dict, \
                                                                var_primitiv, var_primitiv_bd, \
                                                                character_time)
      var_dt = time_integration.set_timestep(config, geom_dict, character_time, var_dt)

      # Spatial gradients and slope limiter
      var_gradient, var_limiter = gradient.gradient_routine(config, dimension_dict, geom_dict,
                                                            metrics_dict, var_primitiv, var_primitiv_bd, 
                                                            var_gradient, var_neig_maxmin, var_limiter)

      # RHS calculation routine
      var_rhs = rhs.rhs_routine(config, dimension_dict, geom_dict, metrics_dict, \
                                gas_property_dict, transport_coefficient_dict, \
                                var_primitiv, var_primitiv_bd, var_gradient, var_limiter, \
                                var_rhs)

      # Time integration and update solution
      var_conserv, \
      var_primitiv, \
      flag_converged_inner = time_integration.time_integration_routine(config, \
                                                                       iteration_inner, \
                                                                       dimension_dict, \
                                                                       geom_dict, \
                                                                       metrics_dict, \
                                                                       gas_property_dict, \
                                                                       transport_coefficient_dict, \
                                                                       var_primitiv, \
                                                                       var_primitiv_bd, \
                                                                       var_gradient, \
                                                                       var_conserv, \
                                                                       var_conserv_prev,  \
                                                                       var_rhs, \
                                                                       var_dt, \
                                                                       var_diagonal, \
                                                                       var_dq)

      # Residuals
      sum_rhs = orbital.display_residual(config, iteration, dimension_dict, var_rhs)


      print('Done inner iteration: ', iteration_inner)

      if flag_converged_inner :
        print('Inner iteration converged ', )
        break


    # Set privous conservative variables
    var_conserv_prev = time_integration.set_conservative_previous(config, var_conserv, var_conserv_prev)

    # Increment
    iteration = iteration + 1

    # Output in restart file
    if iteration%frequency_output_restart == 0:
      orbital.output_restart(config, dimension_dict, geom_dict, iteration, var_conserv, var_conserv_prev)

    # Output in visualization file
    if iteration%frequency_output_postprocess == 0:
      orbital.routine_postprocess(config, iteration, meshnode_dict, meshelem_dict, metrics_dict, gas_property_dict, var_primitiv)

    print('Done iteration: ', iteration)


    if kind_steady_mode == 'steady': 
      flag_converged_outer = orbital.check_convergence_outer(config, flag_converged_outer, sum_rhs)


  # Final results
  orbital.output_restart(config, dimension_dict, geom_dict, iteration, var_conserv, var_conserv_prev)
  #orbital.output_tecplot(config, dimension_list, grid_list, geom_dict, iteration, var_primitiv, var_primitiv_bd, var_gradient, var_limiter)
  orbital.routine_postprocess(config, iteration, meshnode_dict, meshelem_dict, metrics_dict, gas_property_dict, var_primitiv)

  return


if __name__ == '__main__':

  print('Initializing Slow solver')

  # Calling classes
  orbital  = orbital()
  meshdata = meshdata()
  flowfield = flowfield()
  boundary = boundary()
  gradient = gradient()
  rhs = rhs()
  time_integration = time_integration()

  # main 
  main()

  print('Finalizing Slow solver')

  exit()