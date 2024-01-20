#!/usr/bin/env python3

# Program to calculate numerical fluxes on cell interface

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/25

import numpy as np
from orbital.orbital import orbital
from rhs import advection
from rhs import viscous

class rhs(orbital):

  def __init__(self):

    print("Calling class: rhs")

    return


  def initialize_rhs(self, config, dimension_dict, geom_dict):

    print('Setting initial RHS variables')

    num_conserv = dimension_dict['num_conservative']
    num_cell    = geom_dict['num_cell']
    var_rhs     = np.zeros(num_conserv*num_cell).reshape(num_conserv, num_cell)

    self.flux_rhs = np.linspace(0, 0, num_conserv, dtype=float)

    return var_rhs


  def reinitialize_rhs(self, config, var_rhs):

    #np.zeros_like(var_rhs, dtype=float)

    var_rhs[:,:] = 0.0

    return var_rhs


  #@orbital.time_measurement_decorated
  def rhs_routine(self, config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_primitiv_bd, var_gradient, var_limiter, var_rhs):

    # Initialize variables
    var_rhs = self.reinitialize_rhs(config, var_rhs)

    # Advection term
    var_rhs = advection.flux_advection(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, var_primitiv, var_primitiv_bd, var_gradient, var_limiter, var_rhs)

    # Viscous term
    var_rhs = viscous.flux_viscous(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_primitiv_bd, var_gradient, var_rhs)

    # 生成項を考慮するならここが良い
    # var_rhs = source.flux_source(config, var_rhs)

    return var_rhs

