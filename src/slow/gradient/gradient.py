#!/usr/bin/env python3

# ***

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/04/07

import logging
import numpy as np
from slow.gradient import gradient_kernel
from slow.orbital.orbital import orbital

logger = logging.getLogger(__name__)

class gradient(orbital):


  def __init__(self):
    logger.info('Calling class: gradient')

    # Gradient variables

    self.avail_gradient_scheme = ['GG', 'WGG']
    self.avail_slope_limiter = ['minmod', 'none']

  @orbital.time_measurement_decorated
  def initialize_gradient(self, config, dimension_dict, geom_dict):

    logger.info('Setting initial gradient variables')

    num_spatial  = 3
    num_primitiv = dimension_dict['num_primitive']
    num_cell = geom_dict['num_cell']

    var_gradient = np.zeros(num_spatial*num_primitiv*num_cell).reshape(num_spatial ,num_primitiv, num_cell)


    logger.info('--Checking gradient method')
    kind_gradient = str( config['gradient_setting']['kind_gradient'] )
    flag_kind_gradient = False
    for n in range(0, len(self.avail_gradient_scheme) ):
      if kind_gradient == self.avail_gradient_scheme[n]:
        flag_kind_gradient = True
    if not flag_kind_gradient:
      logger.info('Error, gradient scheme is not implemented. Please check the contronl file')
      exit()


    # Slope limiter
    var_limiter      = np.zeros(num_primitiv*num_cell).reshape(num_primitiv, num_cell)
    var_limiter[:,:] = 1.0
    var_neig_maxmin  = np.zeros(2*num_primitiv*num_cell).reshape(2,num_primitiv, num_cell)

    logger.info('--Checking limiter model')
    kind_limiter = str( config['gradient_setting']['kind_limiter'] )
    flag_kind_limiter = False
    for n in range(0, len(self.avail_slope_limiter) ):
      if kind_limiter == self.avail_slope_limiter[n]:
        flag_kind_limiter = True
    if not flag_kind_limiter:
      logger.info('Error, slope limiter is not implemented. Please check the contronl file')
      exit()

    return var_gradient, var_limiter, var_neig_maxmin


  @orbital.time_measurement_decorated
  def gradient_routine(self, config, dimension_dict, geom_dict, metrics_dict, var_primitiv, var_primitiv_bd, var_gradient, var_neig_maxmin, var_limiter):

     # Spatial gradients
    var_gradient = self.get_gradient(config, dimension_dict, geom_dict, \
                                    metrics_dict, var_primitiv, var_primitiv_bd, \
                                     var_gradient)

    # Slope limiter
    flag_muscl = config['gradient_setting']['flag_muscl']
    if flag_muscl :
      var_limiter = self.get_slopelimiter(config, dimension_dict, geom_dict, metrics_dict, \
                                          var_primitiv, var_primitiv_bd, var_gradient, var_neig_maxmin, \
                                          var_limiter)

    return var_gradient, var_limiter


  @orbital.time_measurement_decorated
  def get_gradient(self, config, dimension_dict, geom_dict, metrics_dict, var_primitiv, var_primitiv_bd, var_gradient):

    # 面ループの本体は gradient_kernel に置いてある（numba を掛けるため、
    # 配列とスカラーだけを引数に取る素の関数にしてある）

    num_primitiv   = dimension_dict['num_primitive']

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

    var_gradient = gradient_kernel.accumulate_gradient(
                     num_face, num_face_bd, num_cell, num_primitiv,
                     face2cell, face2cell_bd, virtualcell_bd,
                     area_vec, area_vec_bd, length, length_bd, volume,
                     var_primitiv, var_primitiv_bd, var_gradient)

    return var_gradient


  @orbital.time_measurement_decorated
  def get_slopelimiter(self, config, dimension_dict, geom_dict, metrics_dict, var_primitiv, var_primitiv_bd, var_gradient, var_neig_maxmin, var_limiter):

    #def delta_minmod():
    #  grad_face[:] = grad_face[:] + 1.e-20
    #  grad_face_sign[:] = 0.50*np.sign(grad_face[:])
    #  del_face[:] = (0.50+grad_face_sign[:])*var_max[:] + (0.50-grad_face_sign[:])*var_min[:] - var_tmp[:]
    #  del_face[:] = del_face[:]/grad_face[:]
    #  return

    kind_limiter = str( config['gradient_setting']['kind_limiter'] )

    num_face       = geom_dict['num_face_inner']
    num_face_bd    = geom_dict['num_face_boundary']
    num_cell       = geom_dict['num_cell']
    face2cell      = geom_dict['face2cell_inner']
    face2cell_bd   = geom_dict['face2cell_boundary']
    virtualcell_bd = geom_dict['virtualcell_boundary']

    num_primitiv = dimension_dict['num_primitive']

    area_vec    = metrics_dict['area_vec_inner']
    area_vec_bd = metrics_dict['area_vec_boundary']
    length      = metrics_dict['length_inner']
    length_bd   = metrics_dict['length_boundary']

    # 隣接セルの最大・最小と minmod の面ループは gradient_kernel に置いてある
    var_neig_maxmin = gradient_kernel.find_neighbour_maxmin(
                        num_face, num_face_bd, num_primitiv,
                        face2cell, face2cell_bd,
                        var_primitiv, var_primitiv_bd, var_neig_maxmin)

    if kind_limiter == 'minmod' :
      var_limiter = gradient_kernel.apply_minmod_limiter(
                      num_face, num_face_bd, num_primitiv,
                      face2cell, face2cell_bd, virtualcell_bd,
                      area_vec, area_vec_bd, length, length_bd,
                      var_primitiv, var_gradient, var_neig_maxmin, var_limiter)

    elif kind_limiter == 'none' :
      var_limiter[:,:] = 1.0


    return var_limiter