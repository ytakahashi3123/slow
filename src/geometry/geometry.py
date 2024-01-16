#!/usr/bin/env python3

# ***

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/04

import numpy as np
from orbital.orbital import orbital

class geometry(orbital):


  def __init__(self):
    print("Calling class: geometry")

    # bd_list
    self.id_num_boundary   = 0
    self.id_bd_index_min_x = 1
    self.id_bd_index_max_x = 2
    self.id_bd_index_min_y = 3
    self.id_bd_index_max_y = 4
    self.id_bd_id_attri    = 5
    self.id_bd_name_attri  = 6

    # grid_list
    self.id_num_grid_x = 0
    self.id_num_grid_y = 1
    self.id_num_grid_z = 2
    self.id_num_grid_chain = 3
    self.id_index2chain = 4 
    self.id_chain2index = 5
    self.id_coord_grid = 6
    self.id_boundary_attri = 7

    # geom_list
    self.id_num_face     = 0
    self.id_num_face_bd  = 1
    self.id_num_cell     = 2
    self.id_face2node    = 3 
    self.id_face2cell    = 4 
    self.id_face2node_bd = 5
    self.id_face2cell_bd = 6
    self.id_cell2node    = 7
    self.id_node2cell    = 8
    self.virtualcell_bd  = 9

    # metrics_list
    self.id_area_vec    = 0
    self.id_area_vec_bd = 1
    self.id_lenght      = 2
    self.id_length_bd   = 3 
    self.id_volume      = 4

    # Boundary ID
    self.ID_BOUNDARY_FREESTREAM   = 1
    self.ID_BOUNDARY_AXISYMMETRIC = 6
    self.ID_BOUNDARY_OUTLET       = 10
    self.ID_BOUNDARY_WALL_FIX     = 20

  def read_boundary_data(self, config):

    print( 'Setting boundary data based on boundary file...')

    # Boundary condition (*.inp)
    filename_tmp = config['geometry']['filename_bound']
    print( '--Reading computational boundary data (.inp): ', filename_tmp)
    with open(filename_tmp) as f:
      lines = f.readlines()
    #-- リストとして取得
    lines_strip = [line.strip() for line in lines]

    num_block  = int(lines_strip[0])
    num_domain = int(lines_strip[1])

    lines_ele_tmp = lines_strip[2].split()
    num_coord_x   = int(lines_ele_tmp[0])
    num_coord_y   = int(lines_ele_tmp[1])

    num_boundary  = int(lines_strip[4])
    bd_index_min_x=[]
    bd_index_max_x=[]
    bd_index_min_y=[]
    bd_index_max_y=[]
    bd_id_attri=[]
    bd_name_attri=[]
    for n in range(0, num_boundary):
      lines_ele_tmp = lines_strip[5+n].split()
      bd_index_min_x.append( int( lines_ele_tmp[0] ) )
      bd_index_max_x.append( int( lines_ele_tmp[1] ) )
      bd_index_min_y.append( int( lines_ele_tmp[2] ) )
      bd_index_max_y.append( int( lines_ele_tmp[3] ) )
      bd_id_attri.append(    int( lines_ele_tmp[4] ) )
      bd_name_attri.append( lines_ele_tmp[5] )

    #num_coord_list = [num_coord_x, num_coord_y]
    bd_list = [num_boundary, bd_index_min_x, bd_index_max_x, bd_index_min_y, bd_index_max_y, bd_id_attri, bd_name_attri]

    return bd_list


  def read_geometry_data(self, config):

    print( 'Setting geometry data based on grids file...')

    # Geometry data (*.dat)
    filename_tmp = config['geometry']['filename_geom']
    print( '--Reading computational girs data (.dat): ', filename_tmp)
    with open(filename_tmp) as f:
      lines = f.readlines()
    f.close()

    #-- リストとして取得
    lines_strip = [line.strip() for line in lines]

    #-- Indexes settings
    num_var_per_line = int( config['geometry']['num_varpline'] )

    line_count = 0
    lines_ele_tmp   = lines_strip[line_count].split()
    num_grid_x     = int(lines_ele_tmp[0])
    num_grid_y     = int(lines_ele_tmp[1])
    num_grid_z     = int(lines_ele_tmp[2])
    num_grid_chain = num_grid_x*num_grid_y*num_grid_z
    print( '--Numbers of indexes (x,y,z): ', num_grid_x, num_grid_y, num_grid_z )
    print( '--Number of chain (x*y*z): ', num_grid_chain )


    # indentifying i-j-k indexes
    index2chain = np.zeros(num_grid_x*num_grid_y*num_grid_z).reshape(num_grid_x,num_grid_y,num_grid_z).astype(int)
    chain2index = np.zeros(3*num_grid_chain).reshape(3,num_grid_chain).astype(int)
    i_count_tmp = 0
    for k in range(0,num_grid_z):
      for j in range(0,num_grid_y):
        for i in range(0,num_grid_x):
          index2chain[i,j,k] = i_count_tmp
          chain2index[0,i_count_tmp] = i
          chain2index[1,i_count_tmp] = j
          chain2index[2,i_count_tmp] = k
          i_count_tmp = i_count_tmp + 1


    #-- Grids data
    num_dim = 3
    coord_grid = np.zeros(num_dim*num_grid_chain).reshape(num_dim,num_grid_chain)

    #---- X
    index_init_tmp = line_count + 1
    for i in range(0,num_grid_chain):
      index_tmp       = i//num_var_per_line + index_init_tmp
      mod_tmp         = i%num_var_per_line
      ele_tmp         = lines_strip[index_tmp].split() 
      coord_grid[0,i] = float( ele_tmp[mod_tmp] )
      line_count      = index_tmp

    #---- Y
    index_init_tmp = line_count + 1
    for i in range(0,num_grid_chain):
      index_tmp       = i//num_var_per_line + index_init_tmp
      mod_tmp         = i%num_var_per_line
      ele_tmp         = lines_strip[index_tmp].split() 
      coord_grid[1,i] = float( ele_tmp[mod_tmp] )
      line_count      = index_tmp

    #---- Z
    index_init_tmp = line_count + 1
    for i in range(0,num_grid_chain):
      index_tmp       = i//num_var_per_line + index_init_tmp
      mod_tmp         = i%num_var_per_line
      ele_tmp         = lines_strip[index_tmp].split() 
      coord_grid[2,i] = float( ele_tmp[mod_tmp] )
      line_count     = index_tmp


    # Coordinates of computational node
    num_dim = 3
    coord_node = np.zeros(num_dim*num_grid_x*num_grid_y*num_grid_z).reshape(num_dim,num_grid_x,num_grid_y,num_grid_z)
    for n in range(0,num_grid_chain):
      i = chain2index[0,n]
      j = chain2index[1,n]
      k = chain2index[2,n]
      coord_node[0,i,j,k] = coord_grid[0,n]
      coord_node[1,i,j,k] = coord_grid[1,n]
      coord_node[2,i,j,k] = coord_grid[2,n]


    grid_list       = [num_grid_x, num_grid_y, num_grid_z, num_grid_chain, index2chain, chain2index, coord_grid]
    coord_node_list = coord_node

    return grid_list, coord_node_list


  def set_boundary_attribute(self, config, bd_list, grid_list):

    print( 'Setting boundary attribution and indexes...' )

    num_boundary    = bd_list[0]
    bd_index_min_x  = bd_list[1]
    bd_index_max_x  = bd_list[2]
    bd_index_min_y  = bd_list[3]
    bd_index_max_y  = bd_list[4]
    bd_id_attri     = bd_list[5]

    num_grid_x     = grid_list[0]
    num_grid_y     = grid_list[1]
    num_grid_z     = grid_list[2]
    num_grid_chain = grid_list[3]
    index2chain    = grid_list[4]
    chain2index    = grid_list[5]

    boundary_attri = np.zeros(num_grid_chain).reshape(num_grid_chain).astype(int)


    # Set boundary ID
    print( '--'+str(num_boundary)+' boundary(-ies) found')
    for n in range(0, num_boundary):
      if ( bd_index_min_x[n] == bd_index_max_x[n] ):
        i = bd_index_min_x[n] - 1
        k = 0
        for j in range( bd_index_min_y[n]-1, bd_index_max_y[n] ):
          boundary_attri[ index2chain[i,j,k] ] = bd_id_attri[n]
      elif ( bd_index_min_y[n] == bd_index_max_y[n] ): 
        j = bd_index_min_y[n] - 1
        k = 0
        for i in range( bd_index_min_x[n]-1, bd_index_max_x[n] ):
          boundary_attri[ index2chain[i,j,k] ] = bd_id_attri[n]
      else:
        print( 'Error, check boundary index' )
        exit()


    # Gerating boundary index
    boundary_chain = []
    for n in range(0, num_grid_chain):
      if( boundary_attri[n] != 0 ):
        boundary_chain.append( n )
    boundary_chain = np.array( boundary_chain )

    num_boundary_chain = len( boundary_chain )
    print( '--Number of boundary chains: ',num_boundary_chain)


    # Finding Neighboring cell of boundary cells
    #boundary_neighbor = np.zeros(num_boundary_chain).reshape(num_boundary_chain)
    #boundary_neighbor = boundary_neighbor.astype(int)
    #for n_bd in range(0, num_boundary_chain ):
    #  n_in = boundary_chain[n_bd]
    #  i = chain2index[0,n_in]
    #  j = chain2index[1,n_in]
    #  if i == 0 and j != 0 :
    #    boundary_neighbor[n_bd] = index2chain[i+1,j,k]
    #  elif i == num_coord_x-1 and j != 0 :
    #    boundary_neighbor[n_bd] = index2chain[i-1,j,k]
    #  elif i != 0 and j == 0 :
    #    boundary_neighbor[n_bd] = index2chain[i,j+1,k]
    #  elif i != 0 and j == num_coord_y-1 :
    #    boundary_neighbor[n_bd] = index2chain[i,j-1,k]
    #  else:
    #    boundary_neighbor[n_bd] = 0


    grid_list.append( boundary_attri )

    return grid_list


  def set_geometry_face(self, config, grid_list, coord_node_list):

    print( 'Setting geometry variables...' )

    num_grid_x     = grid_list[0]
    num_grid_y     = grid_list[1]
    num_grid_z     = grid_list[2]
    index2chain    = grid_list[4]
    boundary_attri = grid_list[7]

    coord_node = coord_node_list

    num_dim = 3


    # Count number of faces
    # --Inner faces
    face_id = []
    n_count_tmp = 0
    for k in range(0,num_grid_z):
      for j in range(1,num_grid_y):
        for i in range(1,num_grid_x-1):
          face_id.append(n_count_tmp)
          n_count_tmp = n_count_tmp + 1
    for k in range(0,num_grid_z):
      for i in range(1,num_grid_x):
        for j in range(1,num_grid_y-1):
          face_id.append(n_count_tmp)
          n_count_tmp = n_count_tmp + 1
    face_id  = np.array( face_id )
    num_face = len(face_id)

    #--Boundary faces
    face_bd_id = []
    n_count_tmp = 0
    i=0
    for j in range(1,num_grid_y):
      face_bd_id.append(n_count_tmp)
      n_count_tmp = n_count_tmp + 1
    j=0
    for i in range(1,num_grid_x):
      face_bd_id.append(n_count_tmp)
      n_count_tmp = n_count_tmp + 1
    i=num_grid_y-1
    for j in range(1,num_grid_y):
      face_bd_id.append(n_count_tmp)
      n_count_tmp = n_count_tmp + 1
    j=num_grid_x-1
    for i in range(1,num_grid_x):
      face_bd_id.append(n_count_tmp)
      n_count_tmp = n_count_tmp + 1
    face_bd_id  = np.array( face_bd_id )
    num_face_bd = len(face_bd_id)

    print('--Number of faces (inner and boundary): ',num_face, num_face_bd)


    # Set face to index (=node)
    print('--Set face to node...')
    #face2node[l,m,n], face2node_bd[l,m,n]
    # --l: l=0: 自身のノード, l=1: Face(コネクタin2D)を構成する隣のノード
    # --m: i,j,k
    # --n: 連番（Face ID)
    face2node    = np.zeros(2*num_dim*num_face).reshape(2,num_dim,num_face).astype(int)
    face2node_bd = np.zeros(2*num_dim*num_face_bd).reshape(2,num_dim,num_face_bd).astype(int)

    #--Inner
    n_count_tmp = 0
    for k in range(0,num_grid_z):
      for j in range(1,num_grid_y):
        for i in range(1,num_grid_x-1):
          face2node[0,0,n_count_tmp] = i
          face2node[0,1,n_count_tmp] = j
          face2node[0,2,n_count_tmp] = k
          face2node[1,0,n_count_tmp] = i
          face2node[1,1,n_count_tmp] = j-1
          face2node[1,2,n_count_tmp] = k
          n_count_tmp = n_count_tmp + 1
    for k in range(0,num_grid_z):
      for i in range(1,num_grid_x):
        for j in range(1,num_grid_y-1):
          face2node[0,0,n_count_tmp] = i-1
          face2node[0,1,n_count_tmp] = j
          face2node[0,2,n_count_tmp] = k
          face2node[1,0,n_count_tmp] = i
          face2node[1,1,n_count_tmp] = j
          face2node[1,2,n_count_tmp] = k
          n_count_tmp = n_count_tmp + 1

    #--Boundary faces
    n_count_tmp = 0
    k=0
    i=0
    for j in range(1,num_grid_y):
      face2node_bd[0,0,n_count_tmp] = i
      face2node_bd[0,1,n_count_tmp] = j-1
      face2node_bd[0,2,n_count_tmp] = k
      face2node_bd[1,0,n_count_tmp] = i
      face2node_bd[1,1,n_count_tmp] = j
      face2node_bd[1,2,n_count_tmp] = k
      n_count_tmp = n_count_tmp + 1
    j=0
    for i in range(1,num_grid_x):
      face2node_bd[0,0,n_count_tmp] = i
      face2node_bd[0,1,n_count_tmp] = j
      face2node_bd[0,2,n_count_tmp] = k
      face2node_bd[1,0,n_count_tmp] = i-1
      face2node_bd[1,1,n_count_tmp] = j
      face2node_bd[1,2,n_count_tmp] = k
      n_count_tmp = n_count_tmp + 1
    i=num_grid_x-1
    for j in range(1,num_grid_y):
      face2node_bd[0,0,n_count_tmp] = i
      face2node_bd[0,1,n_count_tmp] = j
      face2node_bd[0,2,n_count_tmp] = k
      face2node_bd[1,0,n_count_tmp] = i
      face2node_bd[1,1,n_count_tmp] = j-1
      face2node_bd[1,2,n_count_tmp] = k
      n_count_tmp = n_count_tmp + 1
    j=num_grid_y-1
    for i in range(1,num_grid_x):
      face2node_bd[0,0,n_count_tmp] = i-1
      face2node_bd[0,1,n_count_tmp] = j
      face2node_bd[0,2,n_count_tmp] = k
      face2node_bd[1,0,n_count_tmp] = i
      face2node_bd[1,1,n_count_tmp] = j
      face2node_bd[1,2,n_count_tmp] = k
      n_count_tmp = n_count_tmp + 1


    # Index from face to cell
    print('--Set indexes from face to cell...')

    cell_id = []
    n_count_tmp = 0
    for k in range(0,num_grid_z):
      for j in range(1,num_grid_y):
        for i in range(1,num_grid_x):
          cell_id.append(n_count_tmp)
          n_count_tmp = n_count_tmp + 1
    cell_id  = np.array( cell_id )
    num_cell = len(cell_id)
    print('--Number of cells: ',num_cell )

    # --Cell to indexes of (i,j,j)
    cell2node  = np.zeros(num_dim*num_cell).reshape(num_dim,num_cell).astype(int)
    node2cell  = np.zeros(num_grid_x*num_grid_y*num_grid_z).reshape(num_grid_x,num_grid_y,num_grid_z).astype(int)

    n_count_tmp = 0
    for k in range(0,num_grid_z):
      for j in range(1,num_grid_y):
        for i in range(1,num_grid_x):
          cell2node[0,n_count_tmp] = i
          cell2node[1,n_count_tmp] = j
          cell2node[2,n_count_tmp] = k
          node2cell[i,j,k] = n_count_tmp
          n_count_tmp = n_count_tmp + 1

    # face2cell[0,:]-->自分自身を構成するセル
    # face2cell[1,:]-->隣のセル
    face2cell    = np.zeros(2*num_face).reshape(2,num_face).astype(int)
    # face2cell_bd[0,:]-->自分自身を構成するセル
    # face2cell_bd[1,:]-->Boudary attributionが入っている
    face2cell_bd = np.zeros(2*num_face_bd).reshape(2,num_face_bd).astype(int)

    #--Inner faces
    n_count_tmp = 0
    for k in range(0,num_grid_z):
      for j in range(1,num_grid_y):
        for i in range(1,num_grid_x-1):
          n = node2cell[i,j,k]
          face2cell[0,n_count_tmp] = n
          n = node2cell[i+1,j,k]
          face2cell[1,n_count_tmp] = n
          n_count_tmp = n_count_tmp + 1
    for k in range(0,num_grid_z):
      for i in range(1,num_grid_x):
        for j in range(1,num_grid_y-1):
          n = node2cell[i,j,k]
          face2cell[0,n_count_tmp] = n
          n = node2cell[i,j+1,k]
          face2cell[1,n_count_tmp] = n
          n_count_tmp = n_count_tmp + 1

    #--Boundary faces
    n_count_tmp = 0
    k=0
    i=0
    for j in range(1,num_grid_y):
      n = node2cell[i+1,j,k]
      face2cell_bd[0,n_count_tmp] = n
      face2cell_bd[1,n_count_tmp] = boundary_attri[ index2chain[i,j,k] ]
      n_count_tmp = n_count_tmp + 1
    j=0
    for i in range(1,num_grid_x):
      n = node2cell[i,j+1,k]
      face2cell_bd[0,n_count_tmp] = n
      face2cell_bd[1,n_count_tmp] = boundary_attri[ index2chain[i,j,k] ]
      n_count_tmp = n_count_tmp + 1
    i=num_grid_x-1
    for j in range(1,num_grid_y):
      n = node2cell[i,j,k]
      face2cell_bd[0,n_count_tmp] = n
      face2cell_bd[1,n_count_tmp] = boundary_attri[ index2chain[i,j,k] ]
      n_count_tmp = n_count_tmp + 1
    j=num_grid_y-1
    for i in range(1,num_grid_x):
      n = node2cell[i,j,k]
      face2cell_bd[0,n_count_tmp] = n
      face2cell_bd[1,n_count_tmp] = boundary_attri[ index2chain[i,j,k] ]
      n_count_tmp = n_count_tmp + 1


    # Virtual cells identification on boundary
    # Freestream, symmetric, outlet (gradient-free): virtual (virtualcell_bd=1)
    # Wall: non-virtual (virtualcell_bd=0)
    virtualcell_bd = np.zeros(num_face_bd).reshape(num_face_bd).astype(int)
    for n_face in range(0,num_face_bd):
      bd_attri = face2cell_bd[1,n_face]
      if bd_attri == self.ID_BOUNDARY_FREESTREAM :
        virtualcell_bd[n_face] = 1
      elif bd_attri == self.ID_BOUNDARY_AXISYMMETRIC :
        virtualcell_bd[n_face] = 1
      elif bd_attri == self.ID_BOUNDARY_OUTLET :
        virtualcell_bd[n_face] = 1
      elif bd_attri == self.ID_BOUNDARY_WALL_FIX :
        virtualcell_bd[n_face] = 0
      else:
        print( 'No boundary ID', 'N_Face:',n_face)
        print( 'Check boundary condition. Program stopped')
        exit()

    #for n in range(0,num_face_bd):
    #  print(n,face2cell_bd[0,n],face2cell_bd[1,n], cell2node[0,face2cell_bd[0,n]], cell2node[1,face2cell_bd[0,n]])
    #exit()
    # 境界の４隅で問題があるかも？要確認!! 2022/03/29

    geom_list = [num_face, num_face_bd, num_cell, face2node, face2cell, face2node_bd, face2cell_bd, cell2node, node2cell, virtualcell_bd]

    return geom_list


  def set_metrics(self, config, coord_node_list, geom_list):

    # Metrics
    # Calculate area vectors
    print( 'Setting metrics...' )

    num_dim = 3

    coord_node = coord_node_list

    num_face     = geom_list[0]
    num_face_bd  = geom_list[1]
    num_cell     = geom_list[2]
    face2node    = geom_list[3]
    face2cell    = geom_list[4]
    face2node_bd = geom_list[5]
    face2cell_bd = geom_list[6]
    cell2node    = geom_list[7]
    node2cell    = geom_list[8]


    area_vec    = np.zeros(10*num_face).reshape(10,num_face)
    area_vec_bd = np.zeros(10*num_face_bd).reshape(10,num_face_bd)
    dz = 1.0
    #　内向きを正としている
    # --Inner faces
    for n in range(0,num_face):
      i0 = face2node[0,0,n]
      j0 = face2node[0,1,n]
      k0 = face2node[0,2,n]
      i1 = face2node[1,0,n]
      j1 = face2node[1,1,n]
      k1 = face2node[1,2,n]
      dx = coord_node[0,i0,j0,k0] - coord_node[0,i1,j1,k1]
      dy = coord_node[1,i0,j0,k0] - coord_node[1,i1,j1,k1]
      dl = np.sqrt( dx**2 + dy**2 )
      area_vec[0,n] = dl*dz
      #area_vec[1,n] = dy/dl
      #area_vec[2,n] =-dx/dl
      area_vec[1,n] =-dy/dl
      area_vec[2,n] = dx/dl
      area_vec[3,n] = 0.0
      area_vec[4,n] = dx/dl
      area_vec[5,n] = dy/dl
      area_vec[6,n] = 0.0
      area_vec[7,n] = 0.0
      area_vec[8,n] = 0.0
      area_vec[9,n] = 0.0

    # --Boundary faces
    for n in range(0,num_face_bd):
      i0 = face2node_bd[0,0,n]
      j0 = face2node_bd[0,1,n]
      k0 = face2node_bd[0,2,n]
      i1 = face2node_bd[1,0,n]
      j1 = face2node_bd[1,1,n]
      k1 = face2node_bd[1,2,n]
      dx = coord_node[0,i0,j0,k0] - coord_node[0,i1,j1,k1]
      dy = coord_node[1,i0,j0,k0] - coord_node[1,i1,j1,k1]
      dl = np.sqrt( dx**2 + dy**2 )
      area_vec_bd[0,n] = dl*dz
      #area_vec_bd[1,n] = dy/dl
      #area_vec_bd[2,n] =-dx/dl
      area_vec_bd[1,n] =-dy/dl
      area_vec_bd[2,n] = dx/dl
      area_vec_bd[3,n] = 0.0
      area_vec_bd[4,n] = dx/dl
      area_vec_bd[5,n] = dy/dl
      area_vec_bd[6,n] = 0.0
      area_vec_bd[7,n] = 0.0
      area_vec_bd[8,n] = 0.0
      area_vec_bd[9,n] = 0.0
      #print(area_vec_bd[0,n],area_vec_bd[1,n], area_vec_bd[2,n],i0,j0)

    # Lengths between face center and cell center
    print('--Calculate lengths...')
    # lenght[0,:]-->自分自身を構成するセル側への距離
    # lenght[1,:]-->隣のセル側への距離
    length    = np.zeros(2*num_face).reshape(2,num_face)
    # lenght_bd[:]-->自分自身を構成するセル側への距離
    length_bd = np.zeros(num_face_bd).reshape(num_face_bd)
    # --Inner faces
    for n_face in range(0,num_face):
      # Face center
      i0 = face2node[0,0,n_face]
      j0 = face2node[0,1,n_face]
      k0 = face2node[0,2,n_face]
      i1 = face2node[1,0,n_face]
      j1 = face2node[1,1,n_face]
      k1 = face2node[1,2,n_face]
      x_face = 0.5*( coord_node[0,i0,j0,k0]+coord_node[0,i1,j1,k1] )
      y_face = 0.5*( coord_node[1,i0,j0,k0]+coord_node[1,i1,j1,k1] )
      z_face = 0.5*( coord_node[2,i0,j0,k0]+coord_node[2,i1,j1,k1] )

      # Cell center
      n_cell = face2cell[0,n_face]
      i = cell2node[0,n_cell]
      j = cell2node[1,n_cell]
      k = cell2node[2,n_cell]
      x_cell = 0.25*( coord_node[0,i,j,k]+coord_node[0,i-1,j,k]+coord_node[0,i,j-1,k]+coord_node[0,i-1,j-1,k] )
      y_cell = 0.25*( coord_node[1,i,j,k]+coord_node[1,i-1,j,k]+coord_node[1,i,j-1,k]+coord_node[1,i-1,j-1,k] )
      z_cell = 0.25*( coord_node[2,i,j,k]+coord_node[2,i-1,j,k]+coord_node[2,i,j-1,k]+coord_node[2,i-1,j-1,k] )

      # Lenght between cell center and face center
      length[0,n_face] = np.sqrt( (x_cell-x_face)**2 + (y_cell-y_face)**2 + (z_cell-z_face)**2 )

      # Counter cell center
      n_cell = face2cell[1,n_face]
      i = cell2node[0,n_cell]
      j = cell2node[1,n_cell]
      k = cell2node[2,n_cell]
      x_cell = 0.25*( coord_node[0,i,j,k]+coord_node[0,i-1,j,k]+coord_node[0,i,j-1,k]+coord_node[0,i-1,j-1,k] )
      y_cell = 0.25*( coord_node[1,i,j,k]+coord_node[1,i-1,j,k]+coord_node[1,i,j-1,k]+coord_node[1,i-1,j-1,k] )
      z_cell = 0.25*( coord_node[2,i,j,k]+coord_node[2,i-1,j,k]+coord_node[2,i,j-1,k]+coord_node[2,i-1,j-1,k] )

      # Lenght between counter cell center and face center
      length[1,n_face] = np.sqrt( (x_cell-x_face)**2 + (y_cell-y_face)**2 + (z_cell-z_face)**2 )

      #print(i,j,length[0,n_face],length[1,n_face])

    # --Boundary faces
    for n_face in range(0,num_face_bd):
      # Face center
      i0 = face2node_bd[0,0,n_face]
      j0 = face2node_bd[0,1,n_face]
      k0 = face2node_bd[0,2,n_face]
      i1 = face2node_bd[1,0,n_face]
      j1 = face2node_bd[1,1,n_face]
      k1 = face2node_bd[1,2,n_face]
      x_face = 0.5*( coord_node[0,i0,j0,k0]+coord_node[0,i1,j1,k1] )
      y_face = 0.5*( coord_node[1,i0,j0,k0]+coord_node[1,i1,j1,k1] )
      z_face = 0.5*( coord_node[2,i0,j0,k0]+coord_node[2,i1,j1,k1] )

      # Cell center
      n_cell = face2cell_bd[0,n_face]
      i = cell2node[0,n_cell]
      j = cell2node[1,n_cell]
      k = cell2node[2,n_cell]
      x_cell = 0.25*( coord_node[0,i,j,k]+coord_node[0,i-1,j,k]+coord_node[0,i,j-1,k]+coord_node[0,i-1,j-1,k] )
      y_cell = 0.25*( coord_node[1,i,j,k]+coord_node[1,i-1,j,k]+coord_node[1,i,j-1,k]+coord_node[1,i-1,j-1,k] )
      z_cell = 0.25*( coord_node[2,i,j,k]+coord_node[2,i-1,j,k]+coord_node[2,i,j-1,k]+coord_node[2,i-1,j-1,k] )

      # Lenght between cell center and face center
      length_bd[n_face] = np.sqrt( (x_cell-x_face)**2 + (y_cell-y_face)**2 + (z_cell-z_face)**2 )


    #  Calculate volume
    print('--Calculate volimes...')
    volume = np.zeros(num_face).reshape(num_face)
    dz = 1.0
    for n_cell in range(0,num_cell):
      i = cell2node[0,n_cell]
      j = cell2node[1,n_cell]
      k = cell2node[2,n_cell]

      x4 = coord_node[0,i-1,j  ,k]
      x3 = coord_node[0,i  ,j,  k]
      x2 = coord_node[0,i  ,j-1,k]
      x1 = coord_node[0,i-1,j-1,k]
      y4 = coord_node[1,i-1,j  ,k]
      y3 = coord_node[1,i  ,j,  k]
      y2 = coord_node[1,i  ,j-1,k]
      y1 = coord_node[1,i-1,j-1,k]
      volume[n_cell] = 0.50*( (x3-x1)*(y4-y2) - (y3-y1)*(x4-x2) ) * dz
      
      #print(volume[n_cell],n_cell)

    metrics_list = [area_vec, area_vec_bd, length, length_bd, volume]

    return metrics_list


  def set_walldisctance(self, config, coord_node_list, geom_list):

    print( 'Calculating distance from wall...' )

    return