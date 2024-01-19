#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2022/03/21

import numpy as np
import time as time
import concurrent.futures
from functools import wraps
from general.general import general

class orbital(general):

  # Constants
  file_control_default = "config.yml"

  # Boundary name
  KIND_BOUNDARY_FREESTREAM         = 'freestream'
  KIND_BOUNDARY_TOTAL_PRESS_TEMP   = 'inlet_total_press_temp'
  KIND_BOUNDARY_MASSFLOWRATE       = 'inlet_massflowrate'
  KIND_BOUNDARY_SYMMETRY           = 'symmetry'
  KIND_BOUNDARY_OUTLET             = 'outlet'
  KIND_BOUNDARY_OUTLET_PRESS_FIXED = 'outlet_press_fixed'
  KIND_BOUNDARY_WALL_FIX           = 'wall'
  KIND_BOUNDARY_WALL_SLIP          = 'wall_slip'
  KIND_BOUNDARY_AMBIENT            = 'ambient'

  # Convergence
  flag_converged_dq = False
  flag_converged_rhs = False

  
  def __init__(self):
    print("Calling class: orbital")

    return


  def set_dimension(self):

    # -- Conservative variables (rho, rho*u, rho*v, rho*w, E)
    # -- Primitive variables (rho, u, v, w, T, p)

    print("Setting dimensions of variables")

    num_conserv  = 5
    num_primitiv = 6

    dimension_dict = { 'num_conservative':num_conserv, \
                       'num_primitive':num_primitiv }

    return dimension_dict


  def set_orbital_parameters(self, dimension_dict):

    print("Setting orbital variables")

    num_conserv = dimension_dict['num_conservative']

    # Variables for convergence check
    self.sum_dq_init  = np.zeros(num_conserv).reshape(num_conserv)
    self.sum_rhs_init = np.zeros(num_conserv).reshape(num_conserv)

    return


  def write_tecplotdata( self, filename, print_message, header, delimiter, comments, output_data ):
    
    print(print_message,':',filename)
    np.savetxt(filename, output_data, header=header, delimiter=delimiter, comments=comments )

    return


  def get_pressure_eos(self, density, gas_constant, temperature):
    # Pressure calculated based on equation of state for ideal gas
    # p=rho*R*T

    pressure = density*gas_constant*temperature

    return pressure


  def get_density_eos(self, pressure, gas_constant, temperature):
    # Density calculated based on equation of state for ideal gas
    # rho=p/(R*T)

    density = pressure/(gas_constant*temperature)

    return density


  def get_speedofsound(self, specfic_heat_ratio, density, pressure):
    # Speed of sound for ideal gas
    # sos = sqrt( gamma*p/rho)

    speedofsound = np.sqrt( specfic_heat_ratio*pressure/density )

    return speedofsound


  def get_total_energy(self, density, specific_heat_volum, temperature, velocity):
    # Total energy: rho*Cv*T + 0.5*rho*U^2

    total_energy = density*specific_heat_volum*temperature + 0.50*density*( velocity[0]**2+velocity[1]**2+velocity[2]**2 )

    return total_energy


  def get_enthalpy(self, specific_heat_volum, density, temperature, velocity, pressure):
    # Specific enthalpy: Cv*T + 0.5*U^2 + p/rho

    enthalpy = specific_heat_volum*temperature + 0.50*( velocity[0]**2+velocity[1]**2+velocity[2]**2 ) + pressure/density

    return enthalpy


  def get_conservative(self, density, velocity, temperature, specific_heat_volum):

    cons_rho = density
    cons_mu = density*velocity[0]
    cons_mv = density*velocity[1]
    cons_mw = density*velocity[2]
    cons_e = self.get_total_energy(density, specific_heat_volum, temperature, velocity)
    conservative = [cons_rho,cons_mu,cons_mv,cons_mw,cons_e]

    return conservative
    

  def get_primitive(self, gas_constant, specific_heat_volum, conservative):
    # Densiry: m=rho
    # Momentum: mu,mv,mw = rho*u, rho*v, rho*w
    # Total energy: E =rho*Cv*T + 0.5*rho*U^2

    prim_rho  =  conservative[0]
    prim_u    =  conservative[1]/conservative[0]
    prim_v    =  conservative[2]/conservative[0]
    prim_w    =  conservative[3]/conservative[0]
    prim_temp =( conservative[4] - 0.50*prim_rho*(prim_u**2 + prim_v**2 + prim_w**2) )/(prim_rho*specific_heat_volum)
    prim_pres = self.get_pressure_eos(prim_rho, gas_constant, prim_temp)

    primtive = [prim_rho, prim_u, prim_v, prim_w, prim_temp, prim_pres]

    return primtive

  
  def make_directory_output(self, config):
    # Make directory

    print('Making directories for output...')

    dir_restart  = config['restart_process']['directory_output']
    self.make_directory(dir_restart)

    dir_result  = config['post_process']['directory_output']
    self.make_directory(dir_result)

    return


  def output_restart(self, config, dimension_dict, geom_dict, iteration, var_conserv, var_conserv_prev):

    dir_restart      = config['restart_process']['directory_output']
    file_restart     = config['restart_process']['file_restart']
    flag_time_series = config['restart_process']['flag_time_series']   # --True: stored individually as time series, False: stored by overwriting
    digid_step       = config['restart_process']['digid_step']

    kind_steady_mode = config['time_integration']['kind_steady_mode']

    num_conserv = dimension_dict['num_conservative']
    num_cell    = geom_dict['num_cell']

    if flag_time_series :
      addfile = '_'+str(iteration).zfill(digid_step)
      filename_tmp = dir_restart + '/' + self.split_file(file_restart,addfile,'.')
    else :
      filename_tmp = dir_restart + '/' + file_restart

    res_init = ' '
    for n in range(0,num_conserv):
      res_init = res_init+' '+str( self.sum_rhs_init[n] )

    print_message = 'Writing restart data'
    header  = '# Restart data \n ' + '# Iteration, Number of cell, Number of conservative variable, Initial residuals \n '+'# '+ str(iteration)+' '+str(num_cell)+' '+ str(num_conserv)+' '+res_init
    output_data = np.c_[ var_conserv[0,:], var_conserv[1,:], var_conserv[2,:], var_conserv[3,:], var_conserv[4,:] ]
    delimiter = '\t'
    comments = ''
    self.write_tecplotdata( filename_tmp, print_message, header, delimiter, comments, output_data )


    if kind_steady_mode == 'unsteady' :
      file_unsteady = config['restart_process']['file_unsteady']
      if flag_time_series :
        addfile = '_'+str(iteration).zfill(digid_step)
        filename_tmp = dir_restart + '/' + self.split_file(file_unsteady,addfile,'.')
      else :
        filename_tmp = dir_restart + '/' + file_unsteady

      print_message = 'Writing unsteady data'
      header  = '# Unsteady data \n ' + '# Iteration, Number of cell, Number of conservative variable, Initial residual \n '+'# '+ str(iteration)+' '+str(num_cell) + ' '+str(num_conserv)+' '+res_init
      output_data = np.c_[ var_conserv_prev[0,0,:], var_conserv_prev[0,1,:], var_conserv_prev[0,2,:], var_conserv_prev[0,3,:], var_conserv_prev[0,4,:],  var_conserv_prev[1,0,:], var_conserv_prev[1,1,:], var_conserv_prev[1,2,:], var_conserv_prev[1,3,:], var_conserv_prev[1,4,:]]
      delimiter = '\t'
      comments = ''
      self.write_tecplotdata( filename_tmp, print_message, header, delimiter, comments, output_data )

    return


  def read_restart(self, config, dimension_dict, var_conserv, var_conserv_prev):

    dir_restart      = config['restart_process']['directory_output']
    file_restart     = config['restart_process']['file_restart']
    flag_time_series = config['restart_process']['flag_time_series']   # --True: stored individually as time series, False: stored by overwriting
    digid_step       = config['restart_process']['digid_step']
    restart_step     = config['restart_process']['restart_step']

    kind_steady_mode = config['time_integration']['kind_steady_mode']

    num_conserv = dimension_dict['num_conservative']

    if flag_time_series :
      addfile = '_s'+str(restart_step).zfill(digid_step)
      filename_tmp = dir_restart + '/' + self.split_file(file_restart,addfile,'.')
    else :
      filename_tmp = dir_restart + '/' + file_restart
    
    # Iteration 
    with open(filename_tmp) as f:
      lines = f.readlines()[2]
    # リストとして取得 
    words = lines.split()
    iteration = int(words[1])
    sum_rhs_init_tmp = np.zeros(num_conserv).reshape(num_conserv)
    for n in range(0,num_conserv):
      sum_rhs_init_tmp[n] = float(words[4+n])

    # Conservative data
    delimiter = None
    comments = '#'
    skiprows = 0 
    data_input = np.loadtxt(filename_tmp, delimiter=delimiter, comments=comments, skiprows=skiprows)
    for n in range(0,num_conserv):
      var_conserv[n,:] = data_input[:,n]


    # Previous conservative for unsteady simulation
    if kind_steady_mode == 'unsteady' :
      file_unsteady = config['restart_process']['file_unsteady']
      if flag_time_series :
        addfile = '_'+str(restart_step).zfill(digid_step)
        filename_tmp = dir_restart + '/' + self.split_file(file_unsteady,addfile,'.')
      else :
        filename_tmp = dir_restart + '/' + file_unsteady

      check_file_exist = self.check_file_exist(filename_tmp)
      if check_file_exist :
        # Unsteady fileが存在する場合
        delimiter = None
        comments = '#'
        skiprows = 0 
        data_input = np.loadtxt(filename_tmp, delimiter=delimiter, comments=comments, skiprows=skiprows)
        for n in range(0,num_conserv):
          var_conserv_prev[0,n,:] = data_input[:,n]
          var_conserv_prev[1,n,:] = data_input[:,n+num_conserv]
      else :
        # Unsteady fileが存在しない場合は最新ステップの保存量で代用
         for n in range(0,num_conserv):
          var_conserv_prev[0,n,:] = var_conserv[n,:]
          var_conserv_prev[1,n,:] = var_conserv[n,:]


    return iteration, sum_rhs_init_tmp, var_conserv, var_conserv_prev


  def output_tecplot(self, config, dimension_dict, grid_list, geom_list, iteration, var_primitiv, var_primitiv_bd, var_gradient, var_limiter):

    coord_grid = grid_list[6]
    var_primitiv_grid, var_gradient_grid, var_limiter_grid = self.convert_var_geom2grid(dimension_dict, grid_list, geom_list, var_primitiv, var_primitiv_bd, var_gradient, var_limiter)
    var_r = var_primitiv_grid[0,:]
    var_u = var_primitiv_grid[1,:]
    var_v = var_primitiv_grid[2,:]
    var_w = var_primitiv_grid[3,:]
    var_t = var_primitiv_grid[4,:]
    var_p = var_primitiv_grid[5,:]
    #var_t_gradx  = var_gradient_grid[0,4,:]
    #var_t_grady  = var_gradient_grid[1,4,:]
    var_t_gradx  = var_limiter_grid[1,:]
    var_t_grady  = var_limiter_grid[4,:]

    dir_result       = config['post_process']['directory_output']
    file_tecplot     = config['post_process']['filename_output_tecplot']
    flag_time_series = config['post_process']['flag_time_series']
    digid_step       = config['post_process']['digid_step']

    if flag_time_series :
      addfile = '_s'+str(iteration).zfill(digid_step)
      filename_tmp = dir_result  + '/' + self.split_file(file_tecplot,addfile,'.')
      header  = 'Variables = x,y,z,dens,u,v,w,temp,pres,t_gradx,t_grady \n zone t=step'+addfile+' i= '+str(grid_list[0])+' j= '+str(grid_list[1])+' k= '+str(grid_list[2])+' f=point'
    else :
      filename_tmp = dir_result  + '/' + file_tecplot
      header  = 'Variables = x,y,z,dens,u,v,w,temp,pres,t_gradx,t_grady \n zone t=step i= '+str(grid_list[0])+' j= '+str(grid_list[1])+' k= '+str(grid_list[2])+' f=point'

    print_message = 'Writing output data'
    output_data = np.c_[ coord_grid[0],
                          coord_grid[1],
                          coord_grid[2],
                          var_r,
                          var_u,
                          var_v,
                          var_w,
                          var_t,
                          var_p,
                          var_t_gradx,
                          var_t_grady
                      ]
    delimiter = '\t'
    comments = ''
    self.write_tecplotdata( filename_tmp, print_message, header, delimiter, comments, output_data )

    return


  def convert_var_geom2grid(self, dimension_dict, grid_list, geom_list, var_primitiv, var_primitiv_bd, var_gradient, var_limiter):

    num_coord_chain = grid_list[3]
    index2chain     = grid_list[4]
    chain2index     = grid_list[5]

    num_face_bd  = geom_list[1]
    num_cell     = geom_list[2]
    face2node_bd = geom_list[5]
    face2cell_bd = geom_list[6]
    cell2node    = geom_list[7]

    num_primitiv = dimension_dict['num_primitive']
    num_spatial  = 3
    var_primitiv_grid = np.zeros(num_primitiv*num_coord_chain).reshape(num_primitiv, num_coord_chain)
    var_gradient_grid = np.zeros(num_spatial*num_primitiv*num_coord_chain).reshape(num_spatial, num_primitiv, num_coord_chain)
    var_limiter_grid  = np.zeros(num_primitiv*num_coord_chain).reshape(num_primitiv, num_coord_chain)

    # Inner loop
    for n_cell in range(0,num_cell):
      i = cell2node[0,n_cell]
      j = cell2node[1,n_cell]
      k = cell2node[2,n_cell]
      n = index2chain[i,j,k]
      var_primitiv_grid[:,n]   = var_primitiv[:,n_cell]
      var_gradient_grid[:,:,n] = var_gradient[:,:,n_cell]
      var_limiter_grid[:,n]  = var_limiter[:,n_cell]
    # Boundary
    for n_face in range(0,num_face_bd):
      n_cell = face2cell_bd[0,n_face]
      i = face2node_bd[0,0,n_face]
      j = face2node_bd[0,1,n_face]
      k = face2node_bd[0,2,n_face]
      n = index2chain[i,j,k]
      #var_primitiv_grid[:,n]   = var_primitiv[:,n_cell]
      var_primitiv_grid[:,n]   = var_primitiv_bd[:,n_face]
      var_gradient_grid[:,:,n] = var_gradient[:,:,n_cell] 
      var_limiter_grid[:,n]  = var_limiter[:,n_cell] 

    return var_primitiv_grid, var_gradient_grid, var_limiter_grid


  def routine_postprocess(self, config, iteration, meshnode_dict, meshelem_dict, metrics_dict, gas_property_dict, var_primitiv):

    # Tecplot (not implemented yet)
    #if ( config['post_process']['flag_output_tecplot'] ):
    #  self.output_tecplot(config, dimension_dict, grid_list, geom_dict, iteration, var_primitiv, var_primitiv_bd, var_gradient, var_limiter)

    # VTK
    if config['post_process']['flag_output_vtk'] :
      # -- Setting variables
      scalar_dict, vector_dict = self.prepare_postprocess(config, gas_property_dict, metrics_dict, var_primitiv)

      # -- Output
      self.write_gmsh2vtk(config, iteration, meshnode_dict, meshelem_dict, scalar_dict=scalar_dict, vector_dict=vector_dict)

    return


  def prepare_postprocess(self, config, gas_property_dict, metrics_dict, var_primitiv):

    specfic_heat_ratio = gas_property_dict['specfic_heat_ratio']

    density      = var_primitiv[0,:]
    temperature  = var_primitiv[4,:]
    pressure     = var_primitiv[5,:]
    velocity     = [ var_primitiv[1,:],var_primitiv[2,:],var_primitiv[3,:] ]
    velocity_mag = np.sqrt( velocity[0]**2 + velocity[1]**2 + velocity[2]**2 )
    speedofsound = self.get_speedofsound(specfic_heat_ratio, density, pressure)
    machnumber   = velocity_mag/speedofsound

    #volume = metrics_dict['volume_cell']
    #cellcenter_list = metrics_dict['coord_cellcenter']
    #cellcenter  = [ cellcenter_list[0],cellcenter_list[1],cellcenter_list[2] ]

    scalar_dict = { 'density': density, 'temperature': temperature, 'pressure': pressure, 'machnumber': machnumber }
    #scalar_dict = { 'density': density, 'temperature': temperature, 'pressure': pressure, 'machnumber': machnumber, 'volume':volume }
    vector_dict = { 'velocity':velocity }
    #vector_dict = { 'velocity':velocity, 'cellcenter':cellcenter }

    return scalar_dict, vector_dict


  def write_gmsh2vtk(self, config, iteration, meshnode_dict, meshelem_dict, **kwargs):

    dir_result       = config['post_process']['directory_output']
    filename_vtk     = config['post_process']['filename_output_vtk']
    flag_time_series = config['post_process']['flag_time_series']
    digid_step       = config['post_process']['digid_step']

    if flag_time_series :
      addfile = '_s'+str(iteration).zfill(digid_step)
      filename_tmp = dir_result  + '/' + self.split_file(filename_vtk,addfile,'.')
    else :
      filename_tmp = dir_result  + '/' + filename_vtk


    num_dim    = meshnode_dict['num_dim']
    num_node   = meshnode_dict['num_node']
    coord_node = meshnode_dict['coord_node']

    num_elem       = meshelem_dict['num_elem']
    type_elem      = meshelem_dict['type_elem']
    elem2node_dict = meshelem_dict['elem2node_dict']
    num_elembytype = meshelem_dict['num_elembytype']
    elem_index_l2g = meshelem_dict['elem_index_l2g' ]
    num_type_elem  = meshelem_dict['num_type_elem']
    num_nodebytype = meshelem_dict['num_nodebytype']


    # 1: points, 2: Triangle, 9: Quadrangle, 10, Tetrahedron, 12: Hexahedron, 13: Triangular prism, 14: Pyramid
    vtk_id_cell_type = {'Point':1, 'Triangle':5, 'Quadrangle':9, 'Tetrahedron':10, 'Hexahedron':12, 'Prism':13, 'Pyramid':14}
    kind_nodebytype  = ['Point', 'Triangle', 'Quadrangle', 'Tetrahedron', 'Hexahedron', 'Prism', 'Pyramid']


    newline_code='\n'
    blank_code=' '
    
    # VTK header
    vtk_header_version='# vtk DataFile Version 2.0'
    vtk_header_title='SLOW results'
    vtk_header_encode='Ascii'
    vtk_header_dataset_type='Dataset Unstructured_grid'

    # Points
    vtk_point_header   = 'Points'+blank_code
    vtk_point_number   = str(num_node)+blank_code
    vtk_point_datatype = 'float'

    # Cells
    vtk_cell_header = 'Cells'+blank_code
    # --count elements
    cell_counts     = 0
    celldata_counts = 0
    flag_cell = [False]*num_type_elem
    for i in range(0,num_type_elem):
      if num_elembytype[i] > 0 and (kind_nodebytype[i] == 'Triangle' or kind_nodebytype[i] == 'Quadrangle'):
        flag_cell[i]    = True
        # Element数
        cell_counts     = cell_counts + num_elembytype[i]
        # Line上に存在する数の総数（三角形の場合：(三角形タイプ)+(要素ID)*３=4こ/1lineになる）
        celldata_counts = celldata_counts + (num_nodebytype[i]+1)*num_elembytype[i]
    vtk_cell_number     = str(cell_counts)+blank_code
    vtk_celldata_number = str(celldata_counts)+blank_code

    # Cell types
    vtk_celltype_header = 'Cell_Types'+blank_code
    vtk_celltype_number = vtk_cell_number

    # Cell data 
    vtk_celldata_scalars_header   = 'Cell_data' + blank_code
    # --scalars
    try:
      flag_scalar_output   = True
      num_flowfield_tmp    = len( kwargs['scalar_dict'] )
      flowfield_data_tmp   = kwargs['scalar_dict']
      flowfield_name_tmp   = list( kwargs['scalar_dict'].keys() )
      vtk_celldata_scalars_name     = 'Scalars'   + blank_code 
      vtk_celldata_scalars_datatype = 'float'
      vtk_celldata_scalars_lookup   = 'Lookup_table default'
    except :
      flag_scalar_output = False   

    # --vectors
    try:
      flag_vector_output = True
      num_vector_tmp            = len( kwargs['vector_dict'] )
      num_vector_direction_tmp  = 3 
      flowfield_vector_tmp      = kwargs['vector_dict']
      flowfield_vector_name_tmp = list( kwargs['vector_dict'].keys() )
      vtk_celldata_vectors_name = 'Vectors' + blank_code 
      vtk_celldata_vectors_datatype = 'float'
    except :
      flag_vector_output = False


    print('Writing VTK file: ', filename_tmp)
    # File open  
    file = open(filename_tmp, "w")

    # --Header
    file.write( vtk_header_version      + newline_code )
    file.write( vtk_header_title        + newline_code )
    file.write( vtk_header_encode       + newline_code )
    file.write( vtk_header_dataset_type + newline_code )

    # --Points
    file.write( vtk_point_header + vtk_point_number + vtk_point_datatype + newline_code )
    for n in range(0,num_node) :
      txt_tmp = str(coord_node[n,0]) + blank_code + str(coord_node[n,1]) + blank_code + str(coord_node[n,2]) + newline_code
      file.write( txt_tmp )

    # --Cells
    file.write( vtk_cell_header + vtk_cell_number + vtk_celldata_number + newline_code )
    for i in range(0,num_type_elem):
      if flag_cell[i] :
        for n in range(0,num_elembytype[i]):
          index_tmp = elem_index_l2g[i][n]
          txt_tmp = str( num_nodebytype[i] ) + blank_code
          for x in elem2node_dict[index_tmp] :
            # VTKでは(Element_ID-1)_gmshで読み込む
            txt_tmp = txt_tmp + str(x-1) + blank_code
          txt_tmp = txt_tmp + newline_code
          file.write( txt_tmp )

    # --Cell types
    file.write( vtk_celltype_header + vtk_celltype_number + newline_code )
    for i in range(0,num_type_elem):
      if flag_cell[i] :
        for n in range(0,num_elembytype[i]):
          txt_tmp = str( vtk_id_cell_type[ kind_nodebytype[i] ] ) + newline_code
          file.write( txt_tmp )

    # --Cell data
    if flag_scalar_output or flag_vector_output :
      file.write( vtk_celldata_scalars_header + vtk_celltype_number + newline_code )
    # ----scalars
    if flag_scalar_output :
      for i in range(0,num_flowfield_tmp):
        str_tmp = flowfield_name_tmp[i]
        file.write( vtk_celldata_scalars_name + flowfield_name_tmp[i] + blank_code + vtk_celldata_scalars_datatype + newline_code )
        file.write( vtk_celldata_scalars_lookup + newline_code )
        #print(str_tmp)
        for n in range(0,cell_counts):
          file.write( str( flowfield_data_tmp[str_tmp][n] ) + newline_code )
    # ----vectors
    if flag_vector_output :
      for i in range(0,num_vector_tmp):
        str_tmp = flowfield_vector_name_tmp[i]
        file.write( vtk_celldata_vectors_name + flowfield_vector_name_tmp[i] + blank_code + vtk_celldata_vectors_datatype + newline_code )
        for n in range(0,cell_counts):
          txt_tmp = ''
          for m in range(0,num_vector_direction_tmp):
            txt_tmp = txt_tmp + str( flowfield_vector_tmp[str_tmp][m][n] ) + blank_code
          txt_tmp = txt_tmp + newline_code
          file.write( txt_tmp )

    # File close
    file.close()

    return


  def check_convergence_inner(self, config, flag_converged_inner, sum_dq):

    criterion_convergence = config['time_integration']['criterion_convergence_innerloop']
    flag_conv_relative    = config['time_integration']['flag_convergence_relative_innerloop']

    if flag_conv_relative :
      if sum_dq[4]/self.sum_dq_init[4] <= criterion_convergence :
        flag_converged_inner = True
    else:
      if sum_dq[4] <= criterion_convergence :
        flag_converged_inner = True

    return flag_converged_inner


  def check_convergence_outer(self, config, flag_converged_outer, sum_rhs):

    criterion_convergence = config['time_integration']['criterion_convergence_outerloop']
    flag_conv_relative    = config['time_integration']['flag_convergence_relative_outerloop']

    if flag_conv_relative :
      if sum_rhs[4] <= criterion_convergence :
        flag_converged_outer = True
    else :
      if sum_rhs[4]/self.sum_rhs_init[4] <= criterion_convergence :
        flag_converged_outer = True

    return flag_converged_outer


  def display_residual(self, config, iteration, dimension_dict, var_rhs):

    num_conserv = dimension_dict['num_conservative']
    sum_rhs = np.sum(var_rhs**2,axis=1)
    print('Residuals: ', [sum_rhs[m] for m in range(0,num_conserv)])
    if iteration == 0:
      self.sum_rhs_init = sum_rhs

    #flag_conv_relative = config['time_integration']['flag_convergence_relative_outerloop']
    #if flag_conv_relative :
    #  sum_rhs = sum_rhs/self.sum_rhs_init

    return sum_rhs


  def display_deltaq(self, config, iteration_inner, dimension_dict, var_dq):

    num_conserv = dimension_dict['num_conservative']
    sum_dq  = np.sum(var_dq**2, axis=1)
    print('Delta Q  : ', [sum_dq[m] for m in range(0,num_conserv)])
    if iteration_inner == 0:
      self.sum_dq_init = sum_dq

    #flag_conv_relative = config['time_integration']['flag_convergence_relative_innerloop']
    #if flag_conv_relative :
    #  sum_dq = sum_dq/self.sum_dq_init

    return sum_dq

  # Decorator for time measurement
  def time_measurement_decorated(func):
    @wraps(func)
    def wrapper(*args, **kargs) :
      flag_time_measurement = True
      if flag_time_measurement :
        start_time = time.time()
        result = func(*args,**kargs)
        elapsed_time = time.time() - start_time
        print('Elapsed time of '+str(func.__name__)+str(':'),elapsed_time,'s')
      else :
        result = func(*args,**kargs)
      return result 
    return wrapper
  
  # Decorator for parallel computation
  #def parallel_execution_decorated(func):
  #  @wraps(func)
  #  def wrapper(*args, **kwargs):
  #    # Create a ThreadPoolExecutor
  #    with concurrent.futures.ThreadPoolExecutor() as executor:
  #      # Submit the function to the executor
  #      result = executor.submit(func, *args, **kwargs).result()
  #    return result
  #  return wrapper

  #def parallel_execution_decorated(max_workers=None):
  #  def decorator(func):
  #    @wraps(func)
  #    def wrapper(*args, **kwargs):
  #      # Create a ThreadPoolExecutor with the specified number of workers
  #      with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
  #        # Submit the function to the executor
  #        result = executor.submit(func, *args, **kwargs).result()
  #        return result
  #    return wrapper
  #  return decorator

  def parallel_execution_decorated(func=None, max_workers=None):
    if func is None:
      return lambda f: parallel_execution_decorated(f, max_workers)
    @wraps(func)
    def wrapper(*args, **kwargs):
      with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        result = executor.submit(func, *args, **kwargs).result()
      return result
    return wrapper