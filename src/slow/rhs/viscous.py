#!/usr/bin/env python3

# Program to calculate viscous flux on interface

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/25

from slow.orbital.orbital import orbital
from slow.rhs import viscous_kernel


@orbital.time_measurement_decorated
def flux_viscous(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_primitiv_bd, var_gradient, var_rhs):

  # 面ループと応力・熱流束の本体は viscous_kernel に置いてある
  # （numba を掛けるため、配列とスカラーだけを引数に取る素の関数にしてある）

  # Input parameters
  num_conserv  = dimension_dict['num_conservative']
  num_primitiv = dimension_dict['num_primitive']

  num_face     = geom_dict['num_face_inner']
  num_face_bd  = geom_dict['num_face_boundary']
  face2cell    = geom_dict['face2cell_inner']
  face2cell_bd = geom_dict['face2cell_boundary']

  area_vec    = metrics_dict['area_vec_inner']
  area_vec_bd = metrics_dict['area_vec_boundary']

  viscosity           = transport_coefficient_dict['viscosity']
  thermal_cond        = transport_coefficient_dict['thermal_conductivity']
  viscosity_bd        = transport_coefficient_dict['viscosity_boundary']
  thermal_cond_bd     = transport_coefficient_dict['thermal_conductivity_boundary']

  var_rhs = viscous_kernel.accumulate_viscous(
              num_face, num_face_bd, num_conserv, num_primitiv,
              face2cell, face2cell_bd, area_vec, area_vec_bd,
              var_primitiv, var_primitiv_bd, var_gradient,
              viscosity, thermal_cond, viscosity_bd, thermal_cond_bd, var_rhs)

  return var_rhs
