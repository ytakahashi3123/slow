#!/usr/bin/env python3

# ***

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/04/07

import numpy as np
from orbital.orbital import orbital
import time as time

class gradient(orbital):


  def __init__(self):
    print("Calling class: gradient")

    # Gradient variables

    self.avail_gradient_scheme = ['GG', 'WGG']
    self.avail_slope_limiter = ['minmod', 'none']


  def initialize_gradient(self, config, dimension_dict, geom_dict):

    print('Setting initial gradient variables')

    num_spatial  = 3
    num_primitiv = dimension_dict['num_primitive']
    num_cell = geom_dict['num_cell']

    var_gradient = np.zeros(num_spatial*num_primitiv*num_cell).reshape(num_spatial ,num_primitiv, num_cell)


    print('--Checking gradient method')
    kind_gradient = str( config['gradient_setting']['kind_gradient'] )
    flag_kind_gradient = False
    for n in range(0, len(self.avail_gradient_scheme) ):
      if kind_gradient == self.avail_gradient_scheme[n]:
        flag_kind_gradient = True
    if not flag_kind_gradient:
      print('Error, gradient scheme is not implemented. Please check the contronl file')
      exit()


    # Slope limiter
    var_limiter      = np.zeros(num_primitiv*num_cell).reshape(num_primitiv, num_cell)
    var_limiter[:,:] = 1.0
    var_neig_maxmin  = np.zeros(2*num_primitiv*num_cell).reshape(2,num_primitiv, num_cell)

    print('--Checking limiter model')
    kind_limiter = str( config['gradient_setting']['kind_limiter'] )
    flag_kind_limiter = False
    for n in range(0, len(self.avail_slope_limiter) ):
      if kind_limiter == self.avail_slope_limiter[n]:
        flag_kind_limiter = True
    if not flag_kind_limiter:
      print('Error, slope limiter is not implemented. Please check the contronl file')
      exit()

    return var_gradient, var_limiter, var_neig_maxmin


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


  def get_gradient(self, config, dimension_dict, geom_dict, metrics_dict, var_primitiv, var_primitiv_bd, var_gradient):

    num_face       = geom_dict['num_face_inner']
    num_face_bd    = geom_dict['num_face_boundary']
    num_cell       = geom_dict['num_cell']
    face2cell      = geom_dict['face2cell_inner']
    #face2node_bd   = geom_list[5]
    face2cell_bd   = geom_dict['face2cell_boundary']
    #cell2node      = geom_dict[7]
    virtualcell_bd = geom_dict['virtualcell_boundary']

    area_vec    = metrics_dict['area_vec_inner']
    area_vec_bd = metrics_dict['area_vec_boundary']
    length      = metrics_dict['length_inner']
    length_bd   = metrics_dict['length_boundary']
    volume      = metrics_dict['volume_cell']

    # Initialize
    var_gradient[:,:,:] = 0.0

    # --Inner loop
    fact_m = 0.5
    fact_p = 0.5
    for n_face in range(0,num_face):
      n_cell_self = face2cell[0,n_face]
      n_cell_neig = face2cell[1,n_face]
      area   = area_vec[0,n_face]
      vec_x  = area*area_vec[1,n_face]
      vec_y  = area*area_vec[2,n_face]
      dl_s   = length[0,n_face]
      dl_n   = length[1,n_face]
      fact_m = dl_s/(dl_s+dl_n)
      fact_p = dl_n/(dl_s+dl_n)
      var_face = fact_p*var_primitiv[:,n_cell_self] + fact_m*var_primitiv[:,n_cell_neig]
      var_gradient[0,:,n_cell_self] = var_gradient[0,:,n_cell_self] - var_face*vec_x
      var_gradient[1,:,n_cell_self] = var_gradient[1,:,n_cell_self] - var_face*vec_y
      var_gradient[0,:,n_cell_neig] = var_gradient[0,:,n_cell_neig] + var_face*vec_x
      var_gradient[1,:,n_cell_neig] = var_gradient[1,:,n_cell_neig] + var_face*vec_y

    # --Boundaryr loop
    fact_m = 0.5
    fact_p = 0.5
    for n_face in range(0,num_face_bd):
      n_cell_self = face2cell_bd[0,n_face]
      area        = area_vec_bd[0,n_face]
      vec_x       = area*area_vec_bd[1,n_face]
      vec_y       = area*area_vec_bd[2,n_face]
      dl_s   = length[0,n_face]
      dl_n   = length[0,n_face]*float(virtualcell_bd[n_face])
      fact_m = dl_s/(dl_s+dl_n)
      fact_p = dl_n/(dl_s+dl_n)
      var_face = fact_p*var_primitiv[:,n_cell_self]+fact_m*var_primitiv_bd[:,n_face]
      var_gradient[0,:,n_cell_self] = var_gradient[0,:,n_cell_self] - var_face*vec_x
      var_gradient[1,:,n_cell_self] = var_gradient[1,:,n_cell_self] - var_face*vec_y


    for n_cell in range(0,num_cell):
      inv_volume = 1.0/volume[n_cell]
      var_gradient[0,:,n_cell] = var_gradient[0,:,n_cell]*inv_volume
      var_gradient[1,:,n_cell] = var_gradient[1,:,n_cell]*inv_volume

    return var_gradient


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

    # Initialize
    var_limiter[:,:] = 1.0
    #grad_face = np.zeros(num_primitiv).reshape(num_primitiv)
    #grad_face_sign = np.zeros(num_primitiv).reshape(num_primitiv)
    #del_face = np.zeros(num_primitiv).reshape(num_primitiv)

    # Searching maximum and minimum variables on neigboring cell
    #start_time = time.time()
    var_neig_maxmin[0,:,:] = var_primitiv[:,:]
    var_neig_maxmin[1,:,:] = var_primitiv[:,:]
    for n_face in range(0,num_face):
      n_cell_self = face2cell[0,n_face]
      n_cell_neig = face2cell[1,n_face]
      #var_neig_maxmin[0,:,n_cell_self] = [max(var_neig_maxmin[0,m,n_cell_self], var_primitiv[m,n_cell_neig]) for m in range(num_primitiv)]
      #var_neig_maxmin[0,:,n_cell_neig] = [max(var_neig_maxmin[0,m,n_cell_neig], var_primitiv[m,n_cell_self]) for m in range(num_primitiv)]
      #var_neig_maxmin[1,:,n_cell_self] = [min(var_neig_maxmin[0,m,n_cell_self], var_primitiv[m,n_cell_neig]) for m in range(num_primitiv)]
      #var_neig_maxmin[1,:,n_cell_neig] = [min(var_neig_maxmin[0,m,n_cell_neig], var_primitiv[m,n_cell_self]) for m in range(num_primitiv)]
      for m in range(0,num_primitiv):
        var_neig_maxmin[0,m,n_cell_self] = max(var_neig_maxmin[0,m,n_cell_self], var_primitiv[m,n_cell_neig]) 
        var_neig_maxmin[0,m,n_cell_neig] = max(var_neig_maxmin[0,m,n_cell_neig], var_primitiv[m,n_cell_self]) 
        var_neig_maxmin[1,m,n_cell_self] = min(var_neig_maxmin[1,m,n_cell_self], var_primitiv[m,n_cell_neig]) 
        var_neig_maxmin[1,m,n_cell_neig] = min(var_neig_maxmin[1,m,n_cell_neig], var_primitiv[m,n_cell_self]) 
    for n_face in range(0,num_face_bd):
      n_cell_self = face2cell_bd[0,n_face]
      for m in range(0,num_primitiv):
        var_neig_maxmin[0,m,n_cell_self] = max(var_neig_maxmin[0,m,n_cell_self], var_primitiv_bd[m,n_face]) 
        var_neig_maxmin[1,m,n_cell_self] = min(var_neig_maxmin[1,m,n_cell_self], var_primitiv_bd[m,n_face]) 
    #elapsed_time = time.time() - start_time
    #print(elapsed_time)


    #start_time = time.time()
    if kind_limiter == 'minmod' :
      # Inner faces
      for n_face in range(0,num_face):
        n_cell_self = face2cell[0,n_face]
        n_cell_neig = face2cell[1,n_face]    
        dl_s   = length[0,n_face]
        dl_n   = length[1,n_face]
        vec_x  = area_vec[1,n_face]*(dl_s+dl_n)
        vec_y  = area_vec[2,n_face]*(dl_s+dl_n)
        vec_z  = area_vec[3,n_face]*(dl_s+dl_n)

        # Slope limiter
        for m in range(0,num_primitiv):
          # --selfside cell
          grad_face =-(var_gradient[0,m,n_cell_self]*vec_x+var_gradient[1,m,n_cell_self]*vec_y+var_gradient[2,m,n_cell_self]*vec_z)
          if grad_face >= 0.0:
            grad_face = grad_face + 1.e-20
          else :
            grad_face = grad_face - 1.e-20
          if grad_face >= 0.0:
            del_face  = var_neig_maxmin[0,m,n_cell_self] - var_primitiv[m,n_cell_self]
          else :
            del_face  = var_neig_maxmin[1,m,n_cell_self] - var_primitiv[m,n_cell_self]
          del_face  =  max( 0.0, min(1.0, del_face/grad_face) )
          var_limiter[m,n_cell_self] = min(var_limiter[m,n_cell_self], del_face)

          # --neigboring side cell
          grad_face = (var_gradient[0,m,n_cell_neig]*vec_x+var_gradient[1,m,n_cell_neig]*vec_y+var_gradient[2,m,n_cell_neig]*vec_z)
          if grad_face >= 0.0:
            grad_face = grad_face + 1.e-20
          else :
            grad_face = grad_face - 1.e-20
          if grad_face >= 0.0:
            del_face  = var_neig_maxmin[0,m,n_cell_neig] - var_primitiv[m,n_cell_neig]
          else :
            del_face  = var_neig_maxmin[1,m,n_cell_neig] - var_primitiv[m,n_cell_neig]
          del_face  =  max( 0.0, min(1.0, del_face/grad_face) )
          var_limiter[m,n_cell_neig] = min(var_limiter[m,n_cell_neig], del_face)


      # Boundary faces
      for n_face in range(0,num_face_bd):
        n_cell_self = face2cell_bd[0,n_face] 
        dl_s   = length_bd[n_face]
        dl_n   = length_bd[n_face]*float(virtualcell_bd[n_face])
        vec_x  = area_vec_bd[1,n_face]*(dl_s+dl_n)
        vec_y  = area_vec_bd[2,n_face]*(dl_s+dl_n)
        vec_z  = area_vec_bd[3,n_face]*(dl_s+dl_n)

        # Slope limiter
        for m in range(0,num_primitiv):
          grad_face =-(var_gradient[0,m,n_cell_self]*vec_x+var_gradient[1,m,n_cell_self]*vec_y+var_gradient[2,m,n_cell_self]*vec_z)
          if grad_face >= 0.0:
            grad_face = grad_face + 1.e-20
          else :
            grad_face = grad_face - 1.e-20
          if grad_face >= 0.0:
            del_face  = var_neig_maxmin[0,m,n_cell_self] - var_primitiv[m,n_cell_self]
          else :
            del_face  = var_neig_maxmin[1,m,n_cell_self] - var_primitiv[m,n_cell_self]
          del_face  =  max( 0.0, min(1.0, del_face/grad_face) )
          var_limiter[m,n_cell_self] = min(var_limiter[m,n_cell_self], del_face)

    elif kind_limiter == 'none' :
      var_limiter[:,:] = 1.0

    #elapsed_time = time.time() - start_time
    #print(elapsed_time)

    return var_limiter