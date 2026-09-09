#!/usr/bin/env python3

# Program to calculate advection flux on interface

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/25

from slow.orbital.orbital import orbital
from slow.rhs import advection_kernel


# config の文字列から、カーネルへ渡す整数の識別子へ
KIND_SCHEME = { 'slau2':  advection_kernel.KIND_SCHEME_SLAU2,
                'haenel': advection_kernel.KIND_SCHEME_HAENEL }


@orbital.time_measurement_decorated
def flux_advection(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, var_primitiv, var_primitiv_bd, var_gradient, var_limiter, var_rhs):

  # 面ループと各スキームの本体は advection_kernel に置いてある
  # （numba を掛けるため、配列とスカラーだけを引数に取る素の関数にしてある）

  eps_muscl   = config['gradient_setting']['eps_muscl']
  kind_scheme = config['numericalflux_setting']['advection_scheme']
  # 未知の名前は従来どおり slau2 として扱う
  kind_scheme = KIND_SCHEME.get(kind_scheme, advection_kernel.KIND_SCHEME_SLAU2)

  # Input parameters
  num_conserv    = dimension_dict['num_conservative']
  num_primitiv   = dimension_dict['num_primitive']

  num_face       = geom_dict['num_face_inner']
  num_face_bd    = geom_dict['num_face_boundary']
  face2cell      = geom_dict['face2cell_inner']
  face2cell_bd   = geom_dict['face2cell_boundary']
  virtualcell_bd = geom_dict['virtualcell_boundary']

  area_vec    = metrics_dict['area_vec_inner']
  area_vec_bd = metrics_dict['area_vec_boundary']
  length      = metrics_dict['length_inner']

  specfic_heat_ratio  = gas_property_dict['specfic_heat_ratio']
  specific_heat_volum = gas_property_dict['specific_heat_volume']

  var_rhs = advection_kernel.accumulate_advection(
              num_face, num_face_bd, num_conserv, num_primitiv, kind_scheme,
              face2cell, face2cell_bd, virtualcell_bd,
              area_vec, area_vec_bd, length,
              specfic_heat_ratio, specific_heat_volum, eps_muscl,
              var_primitiv, var_primitiv_bd, var_gradient, var_limiter, var_rhs)

  return var_rhs
