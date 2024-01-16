#!/usr/bin/env python3

# ***

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/03

import numpy as np
from orbital.orbital import orbital


if __name__ == '__main__':

  # Class orbital
  orbital = orbital()


  # 設定ファイルの読み込み
  # file_control_default = "config_bopt.yml"
  file_control_default = orbital.file_control_default
  arg          = orbital.argument(file_control_default)
  file_control = arg.file
  config       = orbital.read_config_yaml(file_control)


  # Boundary condition (*.inp)
  filename_tmp = config['boundary']['filename_bound']
  print( 'Reading computational boundary (.inp) data: ', filename_tmp)
  with open(filename_tmp) as f:
    lines = f.readlines()
  #-- リストとして取得
  lines_strip = [line.strip() for line in lines]
  
  num_block  = int(lines_strip[0])
  num_domain = int(lines_strip[1])

  lines_count = 1
  num_coord_x_list = []
  num_coord_y_list = []
  id_domain_list   = []
  num_bound_list   = []
  for n in range(0,num_domain):
#    lines_count   = index_init_tmp + 1
    lines_count   = lines_count + 1
    lines_ele_tmp = lines_strip[lines_count].split()
    num_coord_x   = int(lines_ele_tmp[0])
    num_coord_y   = int(lines_ele_tmp[1])
    num_coord_x_list.append( num_coord_x )
    num_coord_y_list.append( num_coord_y )

    lines_count   = lines_count + 1
    lines_ele_tmp = lines_strip[lines_count].split()
    id_domain     = int(lines_ele_tmp[1])
    id_domain_list.append( id_domain )

    lines_count   = lines_count + 1
    lines_ele_tmp = lines_strip[lines_count].split()
    num_bound     = int(lines_ele_tmp[0])
    num_bound_list.append( num_bound )

    flag_loop = False
    i = 0
    while not flag_loop :
      i = i + 1
      if (lines_count+i) >= len(lines_strip): break

      # Reading boundary indexes and their attributions
      lines_ele_tmp = lines_strip[lines_count+i].split()
      if len(lines_ele_tmp) != 5:
        flag_loop = True
        i = i - 1
        break
      index_min_bd_x = lines_ele_tmp[0]
      index_max_bd_x = lines_ele_tmp[1]
      index_min_bd_y = lines_ele_tmp[2]
      index_max_bd_y = lines_ele_tmp[3]
      id_bd_attri    = lines_ele_tmp[4]
    
    lines_count = lines_count + i
  print(id_domain_list, num_bound_list)
  f.close()


  # Reading boundary information (*.fvbnd)
  filename_tmp = config['boundary']['filename_fvbnd']
  print( 'Reading computational bound (.fvbnd) data: ', filename_tmp)
  with open(filename_tmp) as f:
    lines = f.readlines()
  #-- リストとして取得
  lines_strip = [line.strip() for line in lines]

  flag_loop = False
  i = 0
  boundary_name=[]
  bd_attribute=[]
  bd_id=[]
  bd_index_min_x=[]
  bd_index_max_x=[]
  bd_index_min_y=[]
  bd_index_max_y=[]
  bd_index_min_z=[]
  bd_index_max_z=[]
  while not flag_loop :
    i = i + 1
    # Reading boundary names
    lines_ele_tmp = lines_strip[i].split()
    if len(lines_ele_tmp) != 1:
      flag_loop = True
      i = i - 1
      break
    boundary_name.append(lines_ele_tmp)
  
  for n in range(len(boundary_name)+1,len(lines_strip)):
    lines_ele_tmp = lines_strip[n].split()
    bd_attribute.append( lines_ele_tmp[0] )
    bd_id.append( lines_ele_tmp[1] )
    bd_index_min_x.append( lines_ele_tmp[2] )
    bd_index_max_x.append( lines_ele_tmp[3] )
    bd_index_min_y.append( lines_ele_tmp[4] )
    bd_index_max_y.append( lines_ele_tmp[5] )
    bd_index_min_z.append( lines_ele_tmp[6] )
    bd_index_max_z.append( lines_ele_tmp[7] )

  print(bd_attribute,bd_id)
  print(bd_index_min_x,bd_index_max_x)
  print(bd_index_min_y,bd_index_max_y)
  print(bd_index_min_z,bd_index_max_z)
  f.close()
  

  # Read computational grids (*.dat)
  filename_tmp = config['grids']['filename_grids']
  print( 'Reading computational grids data: ', filename_tmp)
  with open(filename_tmp) as f:
    lines = f.readlines()
  #-- リストとして取得
  lines_strip = [line.strip() for line in lines]

  #-- Indexes settings
  num_domain    = int(lines_strip[0])
  print( 'Number of domains: ', num_domain )
  num_var_per_line = int( config['grids']['num_varpline'] )

  line_count = 0
  coord_x_list = []
  coord_y_list = []
  coord_z_list = []
  for n in range(0,num_domain):
    line_count    = line_count + 1
    lines_ele_tmp = lines_strip[line_count].split()
    num_coord_x   = int(lines_ele_tmp[0])
    num_coord_y   = int(lines_ele_tmp[1])
    num_coord_z   = int(lines_ele_tmp[2])
    num_coord_list = num_coord_x*num_coord_y*num_coord_z
    print( 'Numbers of indexes (x,y,z): ', num_coord_x, num_coord_y, num_coord_z )

    #-- Grids data
    coord_x = np.zeros(num_coord_list).reshape(num_coord_list)
    coord_y = np.zeros(num_coord_list).reshape(num_coord_list)
    coord_z = np.zeros(num_coord_list).reshape(num_coord_list)
    #num_line = round( num_coord_x*num_coord_y*num_coord_z/num_var_per_line )

  #---- X
    index_init_tmp = line_count + 1
    for i in range(0,num_coord_list):
      index_tmp  = i//num_var_per_line + index_init_tmp
      mod_tmp    = i%num_var_per_line
      ele_tmp    = lines_strip[index_tmp].split() 
      coord_x[i] = float( ele_tmp[mod_tmp] )
      line_count = index_tmp

    #---- Y
    index_init_tmp = line_count + 1
    for i in range(0,num_coord_list):
      index_tmp  = i//num_var_per_line + index_init_tmp
      mod_tmp    = i%num_var_per_line
      ele_tmp    = lines_strip[index_tmp].split() 
      coord_y[i] = float( ele_tmp[mod_tmp] )
      line_count = index_tmp

    #---- Z
    index_init_tmp = line_count + 1
    for i in range(0,num_coord_list):
      index_tmp  = i//num_var_per_line + index_init_tmp
      mod_tmp    = i%num_var_per_line
      ele_tmp    = lines_strip[index_tmp].split() 
      coord_z[i] = float( ele_tmp[mod_tmp] )
      line_count = index_tmp

    coord_x_list.append(coord_x)
    coord_y_list.append(coord_y)
    coord_z_list.append(coord_z)
    index_init_tmp = line_count + 1

  f.close()


  # Merge computational grids
  print( 'Merging computational grids' )
  #num_coord_x_list
  #num_coord_y_list
  #id_domain_list

  # --それぞれのドメインの要素が同じかどうか確かめる。すべて同じならその方向にはマージしない。
  flag_x_tmp = all(val == num_coord_x_list[0] for val in num_coord_x_list)
  flag_y_tmp = all(val == num_coord_y_list[0] for val in num_coord_y_list)

  if flag_x_tmp: 
    num_coord_x_global = num_coord_x_list[0]
    num_coord_y_global = sum( num_coord_y_list ) - 1*(num_domain-1)
  elif flag_y_tmp:
    num_coord_x_global = sum( num_coord_x_list ) - 1*(num_domain-1)
    num_coord_y_global = num_coord_y_list[0]
  else:
    print( 'Error, indexes are incorrect.' )
    exit()

  print( 'Number of index in x direction (global): ', num_coord_x_global )
  print( 'Number of index in y direction (global): ', num_coord_y_global )

  num_coord_list_global = num_coord_x_global*num_coord_y_global
  coord_x_global = np.zeros(num_coord_list_global).reshape(num_coord_list_global)
  coord_y_global = np.zeros(num_coord_list_global).reshape(num_coord_list_global)
  coord_z_global = np.zeros(num_coord_list_global).reshape(num_coord_list_global)
  #for n in range(0,num_domain):
  coord_x_global = np.append(coord_x_list[0],  coord_x_list[1])
  coord_y_global = np.append(coord_y_list[0],  coord_y_list[1])
  coord_z_global = np.append(coord_z_list[0],  coord_z_list[1])

  # Write grids for vizualization
  if ( config['grids']['flag_output_tecplot'] ):
    filename_tmp = config['grids']['filename_output_tecplot']
    header  = 'Variables = x,y,z,dummy \n zone i= '+str(num_coord_x_global)+' j= '+str(num_coord_y_global)+' k= '+str(num_coord_z)+' f=point'
    cp_data = np.c_[ coord_x_global,
                     coord_y_global,
                     coord_z_global,
                     coord_x_global
                    ]
    np.savetxt(filename_tmp, cp_data, header=header, delimiter='\t', comments='' )

