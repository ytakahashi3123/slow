#!/usr/bin/env python3

# Module to read griddata

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/04/30

#
# "gmsh" package is required
#

import numpy as np
#from platform import python_version
from orbital.orbital import orbital

class meshdata(orbital):

  def __init__(self):

    print("Calling class: meshdata")

    self.id_2d = 2
    self.id_3d = 3

    self.id_type_bar     = 1
    self.id_type_tri     = 2
    self.id_type_quad    = 3
    self.id_type_tet     = 4
    self.id_type_hex     = 5
    self.id_type_prism   = 6
    self.id_type_pyramid = 7
    self.num_type_elem   = 7
    self.num_nodebytype  = [2,3,4,4,8,6,5]
    self.kind_nodebytype  = ['Bar', 'Triangle', 'Quadrangle', 'Tetrahedron', 'Hexahedron', 'Prism', 'Pyramid']

    # Boundary name-->orbital/orbital.py

    #meshnode_dict ID
    self.kind_num_dim    = 'num_dim'
    self.kind_num_node   = 'num_node'
    self.kind_coord_node = 'coord_node'

    # meshelem_dict ID
    self.kind_num_elem       = 'num_elem'
    self.kind_type_elem      = 'type_elem'
    self.kind_dim_elem       = 'dim_elem'
    self.kind_tag_elem       = 'tag_elem'
    self.kind_elem2node_dict = 'elem2node_dict'
    self.kind_num_elembytype = 'num_elembytype'
    self.kind_elem_index_l2g = 'elem_index_l2g' 
    self.kind_num_type_elem  = 'num_type_elem'
    self.kind_num_nodebytype = 'num_nodebytype'

    # geom_dict ID
    self.id_num_face_inner       =  0
    self.id_num_face_boundary    =  1
    self.id_num_cell             =  2
    self.id_face2node_inner      =  3 
    self.id_face2node_boundary   =  4 
    self.id_face2cell_inner      =  5
    self.id_face2cell_boundary   =  6 
    self.id_virtualcell_boundary =  7 
    self.id_face2tag_boundary    =  8
    self.id_tag2name_physics     =  9
    #self.id_cell2dim             =  10
    #self.kind_cell2tag             = 11
    self.kind_num_face_inner       = 'num_face_inner'
    self.kind_num_face_boundary    = 'num_face_boundary'
    self.kind_num_cell             = 'num_cell'
    self.kind_face2node_inner      = 'face2node_inner'
    self.kind_face2node_boundary   = 'face2node_boundary'
    self.kind_face2cell_inner      = 'face2cell_inner'
    self.kind_face2cell_boundary   = 'face2cell_boundary'
    self.kind_cell2tag             = 'cell2tag'
    #self.kind_cell2dim             = 'cell2dim'
    self.kind_virtualcell_boundary = 'virtualcell_boundary'
    self.kind_face2tag_boundary    = 'face2tag_boundary'
    self.kind_tag2name_physics     = 'tag2name_physics'

    # metrics_dict ID
    self.kind_area_vec_inner    = 'area_vec_inner'
    self.kind_area_vec_boundary = 'area_vec_boundary'
    self.kind_length_inner      = 'length_inner'
    self.kind_length_boundary   = 'length_boundary'
    self.kind_volume_cell       = 'volume_cell'
    self.kind_coord_cellcenter  = 'coord_cellcenter'

    return


  def compare_versions(self, version1, version2):
    parts1 = [int(part) for part in version1.split('.')]
    parts2 = [int(part) for part in version2.split('.')]

    # バージョン番号の各部分を比較
    for part1, part2 in zip(parts1, parts2):
        if part1 < part2:
            return -1
        elif part1 > part2:
            return 1

    # 部分的な比較が同じだった場合、長い方が大きい
    if len(parts1) < len(parts2):
        return -1
    elif len(parts1) > len(parts2):
        return 1

    # 完全に同じバージョン
    return 0


  def set_mesh_routine(self, config):

    # Read Gmsh format
    meshnode_dict, meshelem_dict, facecell_list, celltmp_list = self.read_gmsh( config )

    # Set virtual cell
    virtualcell_boundary = self.set_virtualcell_boundary(config, facecell_list)

    # Set geometry 
    geom_dict = self.set_geometry(facecell_list, virtualcell_boundary)

    # Define metrics
    metrics_dict = self.set_metrics(meshnode_dict, meshelem_dict, celltmp_list, geom_dict)

    return meshnode_dict, meshelem_dict, geom_dict, metrics_dict


  def read_gmsh(self, config):
    
    # This function is copied from gmsh-4.9.5-source/tutorials/python/x1.py and modified

    import gmsh
    #import sys

    # File name
    filename = config['meshdata_io']['filename_mesh']

    print('Reading Gmsh file: ', filename)

    # Start GMSH functions
    gmsh.initialize()

    gmsh.open(filename)

    print('Model ' + gmsh.model.getCurrent() + ' (' + str(gmsh.model.getDimension()) + 'D)')

    entities = gmsh.model.getEntities()
    num_node_list   = []
    num_elem_list   = [] 
    for e in entities:
      # Dimension and tag of the entity:
      dim = e[0]
      tag = e[1]

      # Get the mesh nodes for the entity (dim, tag):
      nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(dim, tag)

      # Get the mesh elements for the entity (dim, tag):
      elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)
    
      # Definitions
      num_node_list.append( len( nodeTags) )
      num_elem_list.append( sum(len(i) for i in elemTags) )

    # Dimension
    num_dim = gmsh.model.getDimension()
    print('Dimension:', num_dim)

    # Nodes
    print('Number of nodes: '   , num_node_list)
    print('Number of elements: ', num_elem_list)

    # Node coordinateを取得
    print('Getting node data...')
    num_node = sum(num_node_list)
    coord_node = np.zeros(num_node*3).reshape(num_node,3)
    dim_node = np.zeros(num_node).reshape(num_node).astype(int)
    tag_node = np.zeros(num_node).reshape(num_node).astype(int)
    for n in range(0,num_node):
      nodecoord_tmp, parametriccoord_tmp, dim_tmp, tag_tmp = gmsh.model.mesh.getNode(n+1)
      coord_node[n,0] = nodecoord_tmp[0]
      coord_node[n,1] = nodecoord_tmp[1]
      coord_node[n,2] = nodecoord_tmp[2]
      dim_node[n] =  dim_tmp # 使わない？
      tag_node[n] =  tag_tmp # 使わない？
      #print('Nodes: ', coord_node[n,0], coord_node[n,1], coord_node[n,2], dim_tmp, tag_tmp)

    # Elements情報を取得
    # Types: bar(1) tri(2) quad(3) tet(4) hex(5) prism(6) pyramid(7)
    print('Getting element data...')
    num_elem       = sum(num_elem_list)
    type_elem      = np.zeros(num_elem).reshape(num_elem).astype(int)
    dim_elem       = np.zeros(num_elem).reshape(num_elem).astype(int)
    tag_elem       = np.zeros(num_elem).reshape(num_elem).astype(int)
    num_elembytype = np.zeros(self.num_type_elem).reshape(self.num_type_elem).astype(int)
    nodetag_elem_list = []
    index_l2g_tmp  = [ [],[],[],[],[],[],[] ]
    for n in range(0,num_elem):
      elemtype_tmp, nodetag_elem_tmp, elemdim_tmp, elemtag_tmp = gmsh.model.mesh.getElement(n+1)
      type_elem[n] = elemtype_tmp
      dim_elem[n]  = elemdim_tmp
      tag_elem[n]  = elemtag_tmp
      # それぞれのElementタイプに対して最大配列サイズを確定する。
      # およびそれぞれのElementタイプのローカルなインデクス-->グローバルなインデクスを取得する。
      #（他により良い実装はある気がする）
      if type_elem[n] == self.id_type_bar :
        m = type_elem[n]-1
        num_elembytype[m] = num_elembytype[m] + 1
        index_l2g_tmp[m].append(n)
      elif type_elem[n] == self.id_type_tri :
        m = type_elem[n]-1
        num_elembytype[m] = num_elembytype[m] + 1
        index_l2g_tmp[m].append(n)
      elif type_elem[n] == self.id_type_quad :
        m = type_elem[n]-1
        num_elembytype[m] = num_elembytype[m] + 1
        index_l2g_tmp[m].append(n)
      elif type_elem[n] == self.id_type_tet :
        m = type_elem[n]-1
        num_elembytype[m] = num_elembytype[m] + 1
        index_l2g_tmp[m].append(n)
      elif type_elem[n] == self.id_type_hex :
        m = type_elem[n]-1
        num_elembytype[m] = num_elembytype[m] + 1
        index_l2g_tmp[m].append(n)
      elif type_elem[n] == self.id_type_prism :
        m = type_elem[n]-1
        num_elembytype[m] = num_elembytype[m] + 1
        index_l2g_tmp[m].append(n)
      elif type_elem[n] == self.id_type_pyramid :
        m = type_elem[n]-1
        num_elembytype[m] = num_elembytype[m] + 1
        index_l2g_tmp[m].append(n)
      #print('Elements:', '-type',elemtype_tmp, '-dim.', elemdim_tmp, '-physics.', elemtag_tmp)

      # 全体の要素情報を記憶しておく
      nodetag_elem_list.append(nodetag_elem_tmp)
    #print( num_elembytype, index_l2g_tmp)


    # Getting Node ID on each element
    # 面積計算がアドホックなので３Dにするときは要注意
    print('Node ID on each element:')
    if num_dim == self.id_2d :
      elem2node_dict = {}
      volume_elem = np.zeros(num_elem).reshape(num_elem)
      coord_cellcenter_elem = np.zeros(3*num_elem).reshape(3,num_elem)
      dz = 1.0

      for i in range(0,self.num_type_elem):
        for n in range(0,num_elembytype[i]):
          elemtype_tmp, nodetag_elem_tmp, elemdim_tmp, elemtag_tmp = gmsh.model.mesh.getElement(index_l2g_tmp[i][n]+1)

          # 要素の並びが右回りか左回りか判定するためにここで面積(volume)計算も行っている
          # Volumeが負のときは右回りになるので並び替えはしない。正のときは左回りと判定して並び替えを行う
          if self.kind_nodebytype[i] == 'Triangle' :
            # Triangles
            m0 = int(nodetag_elem_tmp[0])-1
            m1 = int(nodetag_elem_tmp[1])-1
            m2 = int(nodetag_elem_tmp[2])-1
            vec01 = coord_node[m1,:]-coord_node[m0,:]
            vec02 = coord_node[m2,:]-coord_node[m0,:]
            veccross = np.cross(vec01,vec02)
            volume_tmp = 0.50*veccross*dz
            # Reordering
            #if any((x > 0 for x in volume_tmp)) :
            if volume_tmp[2] > 0.0 :
              nodetag_elem_tmp[0] = m0 + 1
              nodetag_elem_tmp[1] = m2 + 1
              nodetag_elem_tmp[2] = m1 + 1
            # Volume
            volume_elem[index_l2g_tmp[i][n]] = np.linalg.norm( volume_tmp )
            # Cell center coordinate
            coord_cellcenter_elem[:,index_l2g_tmp[i][n]] = (coord_node[m0,:]+coord_node[m1,:]+coord_node[m2,:])/3.0

          elif self.kind_nodebytype[i] == 'Quadrangle' :
            # Quadrnngle
            m0 = int(nodetag_elem_tmp[0])-1
            m1 = int(nodetag_elem_tmp[1])-1
            m2 = int(nodetag_elem_tmp[2])-1
            m3 = int(nodetag_elem_tmp[3])-1
            vec13 = coord_node[m1,:]-coord_node[m3,:]
            vec20 = coord_node[m2,:]-coord_node[m0,:]
            veccross = np.cross(vec13,vec20)
            volume_tmp = 0.50*veccross*dz
            # Reordering
            #if any((x > 0 for x in volume_tmp)) :
            if volume_tmp[2] > 0.0 :
              nodetag_elem_tmp[0] = m0 + 1
              nodetag_elem_tmp[1] = m3 + 1
              nodetag_elem_tmp[2] = m2 + 1
              nodetag_elem_tmp[3] = m1 + 1
            # Volume
            volume_elem[index_l2g_tmp[i][n]] = np.linalg.norm( volume_tmp )
            # Cell center coordinate
            coord_cellcenter_elem[:,index_l2g_tmp[i][n]] = (coord_node[m0,:]+coord_node[m1,:]+coord_node[m2,:]+coord_node[m3,:])*0.25
          else :
            pass

          # Elem to node indexes
          #elem2node_dict[ index_l2g_tmp[i][n] ] = nodetag_elem_tmp
          elem2node_dict[ index_l2g_tmp[i][n] ] = [int(i) for i in nodetag_elem_tmp]

    else :
      print('3D mesh is NOT implemented yet.')
      print('Program stopped.')
      exit()


    # Physical groupの情報を取得
    print('Getting physical group data...')
    physicalTags       = gmsh.model.getPhysicalGroups()
    physicalTags_array = np.array(physicalTags)
    num_physics  = len(physicalTags_array[:,0])
    dim_physics  = []
    tag_physics  = []
    name_physics = []
    for n in range(0,num_physics):
      dim_tmp  = physicalTags_array[n,0]
      tag_tmp  = physicalTags_array[n,1]
      name_tmp = gmsh.model.getPhysicalName(dim_tmp,tag_tmp)
      #
      dim_physics.append( dim_tmp )
      tag_physics.append( tag_tmp )
      name_physics.append( name_tmp )
      print('Physical group:', '-name', name_tmp, '-dim.', dim_tmp, '-tag.', tag_tmp)
    # Tag-->Nameの辞書作成
    tag2name_physics_dict = dict()
    for n in range(0,num_physics) :
      tag2name_physics_dict[tag_physics[n]] = name_physics[n]
    print('Dictionary of physical tag-->name:', tag2name_physics_dict)
    
    
    # Inner elementとBoundary elementを区別する。
    # 2Dの場合はInner elementタイプはtri(2) quad(3), Boundaryはbar(1)だけのはず
    # 3Dの場合は上記とはことなる（3Dはあとで実装する・・かも）
    print('Face and cell settings: ', str(num_dim)+'-Dimension case')
    if num_dim == self.id_2d :

      # ElementをCellに読み替える
      num_cell_bar  = num_elembytype[self.id_type_bar-1]
      num_cell_tri  = num_elembytype[self.id_type_tri-1]
      num_cell_quad = num_elembytype[self.id_type_quad-1]
      num_cell      = num_cell_tri + num_cell_quad

      print('Number of bar cell(?)    : ', num_cell_bar )
      print('Number of triangle cell  : ', num_cell_tri )
      print('Number of quadrangle cell: ', num_cell_quad)
      print('Number of cell           : ', num_cell)


      # Face-barの数はnum_cell_barと同じ。境界フェイズになる。オーバーラップもないはず
      num_face_bar = num_cell_bar
      print('Number of boudanry faces : ',num_face_bar)
      face2node_boundary = np.zeros(2*num_face_bar).reshape(2,num_face_bar).astype(int)
      face2tag_boundary  = np.zeros(num_face_bar).reshape(num_face_bar).astype(int)
      n_count_bar_tmp = 0
      for n in range(0, num_face_bar) :
        index_tmp = index_l2g_tmp[self.id_type_bar-1][n]
        m0 = elem2node_dict[ index_tmp ][0]-1
        m1 = elem2node_dict[ index_tmp ][1]-1
        face2node_boundary[0,n_count_bar_tmp] = m0
        face2node_boundary[1,n_count_bar_tmp] = m1
        face2tag_boundary[n_count_bar_tmp]    = tag_elem[ index_tmp ]
        #face2cell_boundary[0,n_count_bar_tmp] = index_l2g_tmp[self.id_type_bar-1][n]
        n_count_bar_tmp = n_count_bar_tmp + 1
        

      # フェイスの数(オーバーラップを許す)=(三角形Elementは3辺)+(四角形Elementは4辺)
      # オーバーラップ分はあとで取り除く必要がある
      num_face_overlap = num_cell_tri*3 + num_cell_quad*4
      print('Number of faces including overlapping faces: ',num_face_overlap)

      # 2Dの場合フェイスは線分で構成されているとする。線分には2つのノードがある。
      face2node_overlap = np.zeros(2*num_face_overlap).reshape(2,num_face_overlap).astype(int)
      face2cell_overlap = np.zeros(2*num_face_overlap).reshape(2,num_face_overlap).astype(int)
      cell2dim          = np.zeros(num_cell).reshape(num_cell).astype(int)
      cell2tag          = np.zeros(num_cell).reshape(num_cell).astype(int)
      volume_cell       = np.zeros(num_cell).reshape(num_cell)
      coord_cellcenter  = np.zeros(3*num_cell).reshape(3,num_cell)

      n_count_face_tmp = 0
      n_count_cell_tmp = 0
      # Triangle element
      for n in range(0, num_cell_tri) :
        index_tmp = index_l2g_tmp[self.id_type_tri-1][n]
        m0 = elem2node_dict[ index_tmp ][0]-1
        m1 = elem2node_dict[ index_tmp ][1]-1
        m2 = elem2node_dict[ index_tmp ][2]-1
        face2node_overlap[0,n_count_face_tmp  ] = m0
        face2node_overlap[1,n_count_face_tmp  ] = m1
        face2node_overlap[0,n_count_face_tmp+1] = m1
        face2node_overlap[1,n_count_face_tmp+1] = m2
        face2node_overlap[0,n_count_face_tmp+2] = m2
        face2node_overlap[1,n_count_face_tmp+2] = m0

        face2cell_overlap[0,n_count_face_tmp:n_count_face_tmp+3] = n_count_cell_tmp #index_tmp
        cell2dim[ n_count_cell_tmp ]         = dim_elem[ index_tmp ]
        cell2tag[ n_count_cell_tmp ]         = tag_elem[ index_tmp ]
        volume_cell[n_count_cell_tmp]        = volume_elem[ index_tmp ]
        coord_cellcenter[:,n_count_cell_tmp] = coord_cellcenter_elem[:, index_tmp ]

        n_count_face_tmp = n_count_face_tmp + 3
        n_count_cell_tmp = n_count_cell_tmp + 1

      # Quadrangle element
      for n in range(0, num_cell_quad) :
        index_tmp = index_l2g_tmp[self.id_type_quad-1][n]
        m0 = elem2node_dict[ index_tmp ][0]-1
        m1 = elem2node_dict[ index_tmp ][1]-1
        m2 = elem2node_dict[ index_tmp ][2]-1
        m3 = elem2node_dict[ index_tmp ][3]-1
        face2node_overlap[0,n_count_face_tmp  ] = m0
        face2node_overlap[1,n_count_face_tmp  ] = m1
        face2node_overlap[0,n_count_face_tmp+1] = m1
        face2node_overlap[1,n_count_face_tmp+1] = m2
        face2node_overlap[0,n_count_face_tmp+2] = m2
        face2node_overlap[1,n_count_face_tmp+2] = m3
        face2node_overlap[0,n_count_face_tmp+3] = m3
        face2node_overlap[1,n_count_face_tmp+3] = m0

        face2cell_overlap[0,n_count_face_tmp:n_count_face_tmp+4] = n_count_cell_tmp #index_tmp
        cell2dim[ n_count_cell_tmp ]         = dim_elem[ index_tmp ]
        cell2tag[ n_count_cell_tmp ]         = tag_elem[ index_tmp ]
        volume_cell[n_count_cell_tmp]        = volume_elem[ index_tmp ]
        coord_cellcenter[:,n_count_cell_tmp] = coord_cellcenter_elem[:, index_tmp ]

        n_count_face_tmp = n_count_face_tmp + 4
        n_count_cell_tmp = n_count_cell_tmp + 1

      print('Face to node:')
      print(face2node_overlap)

      # Check
      if num_face_overlap != n_count_face_tmp:
        print('Error. Check meshio.py')
        print('Program stopped.')
        exit()
      
      # Version check
      #version_current = python_version()
      #version_3_7_0 = "3.7.0" # 3.7未満ではdict型にインデクスがつかない（OrderedDict()ならつく？）
      #comp_version = self.compare_versions(version_current, version_3_7_0)

      # Dict型を使った高速検索を使うとき
      flag_search_dict = True

      if flag_search_dict:
        # Dict()を使うことで検索の計算量をO(n)にしている。
        # 重複したフェイスの同定-->Dict型で高速化する。
        print('Finding overlapping nodes (inner)...')
        flag_face_overlapped     = [False]*num_face_overlap
        flag_face_overlapped_opp = [False]*num_face_overlap
        face_opposite_overlapped = [-1]*num_face_overlap
        # 要素を辞書に作成する
        element_dict = {(face2node_overlap[0, n], face2node_overlap[1, n]): n for n in range(num_face_overlap)}
        # オーバーラップのインデクスを探索
        for n in range(0,num_face_overlap) :
          if flag_face_overlapped_opp[n] :
            continue
          node_s0 = face2node_overlap[0,n]
          node_s1 = face2node_overlap[1,n]
          # 辞書から直接インデックスを取得
          extracted_value = element_dict.get((node_s1, node_s0), -1)
          if extracted_value != -1:
            flag_face_overlapped[n] = True
            flag_face_overlapped_opp[extracted_value] = True
            face_opposite_overlapped[n] = extracted_value
      
        # Boundary faceと隣接セルの同定
        # 0-->隣接セル、1-->Boundary ID
        print('Finding overlapping nodes (boundary)...')
        face2cell_boundary = np.zeros(2*num_face_bar).reshape(2,num_face_bar).astype(int)

        element_dict = {(face2node_overlap[0, n], face2node_overlap[1, n]): n for n in range(num_face_overlap)}
        #face2cell_boundary = np.zeros((2, num_face_bar), dtype=int)
        for m in range(num_face_bar):
          node_n0 = face2node_boundary[0, m]
          node_n1 = face2node_boundary[1, m]
          #target_element = (node_n0, node_n1) #if node_n0 < node_n1 else (node_n1, node_n0)

          #n = element_dict.get((node_n0, node_n1), -1)
          n = element_dict.get((node_n0, node_n1) , element_dict.get((node_n1, node_n0), -1))
          if n != -1 and not flag_face_overlapped[n]:
            node_s0, node_s1 = face2node_overlap[0, n], face2node_overlap[1, n]
            face2cell_boundary[0, m] = face2cell_overlap[0, n]
            face2cell_boundary[1, m] = face2tag_boundary[m]
        
            if node_n0 != node_s0 or node_n1 != node_s1:
              face2node_boundary[0, m] = node_n1
              face2node_boundary[1, m] = node_n0


      elif not flag_search_dict :
        # １つ１つ調べ上げる。こちらの処理の計算量は最大でO(n^2)に注意。nは格子数
        print("Python version is lower than 3.7. The computational cost becomes high.")

        # 重複したフェイスの同定-->さかのぼって地道に調べる
        print('Finding overlapping nodes (inner)...')
        flag_face_overlapped     = [False]*num_face_overlap
        flag_face_overlapped_opp = [False]*num_face_overlap
        face_opposite_overlapped = [-1]*num_face_overlap
        # 要素を辞書に作成する
        for n in range(0,num_face_overlap) :
          node_s0 = face2node_overlap[0,n]
          node_s1 = face2node_overlap[1,n]
          if flag_face_overlapped_opp[n] :
            continue
          for m in range(n+1,num_face_overlap) :
            node_n0 = face2node_overlap[0,m]
            node_n1 = face2node_overlap[1,m]
            if node_n0 == node_s1 and node_n1 == node_s0 \
            or node_n0 == node_s0 and node_n1 == node_s1:
              # 重複が判明したフェイスにはTrue判定
              flag_face_overlapped[n] = True
              flag_face_overlapped_opp[m] = True
              # 反対側のFace idを記録しておく
              face_opposite_overlapped[n] = m
              #print('Overlap cell', face2cell_overlap[0,n]+1, face2cell_overlap[0,m]+1, '-nodes:',node_s0+1, node_s1+1, node_n0+1, node_n1+1)
              break

        # Boundary faceと隣接セルの同定
        # 0-->隣接セル、1-->Boundary ID
        print('Finding overlapping nodes (boundary)...')
        face2cell_boundary = np.zeros(2*num_face_bar).reshape(2,num_face_bar).astype(int)
      #for n in range(0,num_face_overlap) :
      #  # 境界フェイスは重複はしない
      #  if not flag_face_overlapped[n]:
      #    node_s0 = face2node_overlap[0,n]
      #    node_s1 = face2node_overlap[1,n]
      #    for m in range(0, num_face_bar) :
      #      node_n0 = face2node_boundary[0,m]
      #      node_n1 = face2node_boundary[1,m]
      #      if node_n0 == node_s1 and node_n1 == node_s0 \
      #      or node_n0 == node_s0 and node_n1 == node_s1:
      #        face2cell_boundary[0,m] = face2cell_overlap[0,n]
      #        face2cell_boundary[1,m] = face2tag_boundary[m]
      #        # Inner faceとBoundary faceでNodeの方向が一致してない場合は、Boundary faceのノードを入れ替える
      #        if node_n0 == node_s1 and node_n1 == node_s0: 
      #          face2node_boundary[0,m] = node_n1
      #          face2node_boundary[1,m] = node_n0
      #        print('Cell',face2cell_boundary[0,m]+1,'-boundary attribute:',face2cell_boundary[1,m], '-nodes:',node_s0+1, node_s1+1, node_n0+1, node_n1+1)
        for m in range(0, num_face_bar) :
          node_n0 = face2node_boundary[0,m]
          node_n1 = face2node_boundary[1,m]
          for n in range(0,num_face_overlap) :
            # 境界フェイスは重複はしない
            if not flag_face_overlapped[n]:
              node_s0 = face2node_overlap[0,n]
              node_s1 = face2node_overlap[1,n]
              if node_n0 == node_s1 and node_n1 == node_s0 \
              or node_n0 == node_s0 and node_n1 == node_s1:
                face2cell_boundary[0,m] = face2cell_overlap[0,n]
                face2cell_boundary[1,m] = face2tag_boundary[m]
                # Inner faceとBoundary faceでNodeの方向が一致してない場合は、Boundary faceのノードを入れ替える
                if node_n0 == node_s1 and node_n1 == node_s0: 
                  face2node_boundary[0,m] = node_n1
                  face2node_boundary[1,m] = node_n0
                break
          #print('Cell',face2cell_boundary[0,m]+1,'-boundary attribute:',face2cell_boundary[1,m], '-nodes:',node_s0+1, node_s1+1, node_n0+1, node_n1+1)

      else :
        print("Please check meshdata.py")
        print("Program stopped.")
        exit()


      # 重複したフェイスの削除
      # --重複していないフェイスの数(Inner face+Boundary face)
      num_face_merge = sum( flag_tmp==False for flag_tmp in flag_face_overlapped) 
      print('Number of faces including inner and boundary faces: ', num_face_merge)
      # -重複していたフェイスの数=Inner faceの数
      num_face_inner = sum( flag_tmp==True for flag_tmp in flag_face_overlapped) 
      print('Number of inner faces (=number of faces that already overlapped: ', num_face_inner)
      # -重複していないフェイスの数から重複していたフェイスの数を差し引くと、Boundary faceの数となる
      num_face_boundary = num_face_merge - num_face_inner
      print('Number of boundary faces: ', num_face_boundary)

      face2node_inner = np.zeros(2*num_face_inner).reshape(2,num_face_inner).astype(int)
      face2cell_inner = np.zeros(2*num_face_inner).reshape(2,num_face_inner).astype(int)
      n_count_face_tmp = 0
      index_inner_face =[]
      for n in range(0,num_face_overlap) :
        if flag_face_overlapped[n] :
          # Inner cell
          face2node_inner[0,n_count_face_tmp] = face2node_overlap[0,n]
          face2node_inner[1,n_count_face_tmp] = face2node_overlap[1,n]
          face2cell_inner[0,n_count_face_tmp] = face2cell_overlap[0,n]
          face2cell_inner[1,n_count_face_tmp] = face2cell_overlap[0,face_opposite_overlapped[n]]
         # print('Self cell: ',face2cell_inner[0,n_count_face_tmp]+1, '   ','Opposite cell: ', face2cell_inner[1,n_count_face_tmp]+1)
          n_count_face_tmp = n_count_face_tmp+1     

      #for n in range(0,num_face_inner):
      #  print(n,'Cell',face2cell_inner[0,n]+1, face2cell_inner[1,n]+1, 'Nodes',face2node_inner[0,n]+1, face2node_inner[1,n]+1)
      #for n in range(0,num_face_boundary):
      #  print(n,'Cell',face2cell_boundary[0,n]+1, face2cell_boundary[1,n]+1, 'Nodes',face2node_boundary[0,n]+1, face2node_boundary[1,n]+1)
      #exit()

      # Check
      if num_face_inner != n_count_face_tmp:
        print('Error. Check meshio.py: (number of inner face)')
        print('Program stopped.')
        exit()
      if num_face_bar != num_face_boundary:
        print('Error. Check meshio.py: (number of boundary face)')
        print('Program stopped.')
        exit()

    else :
      print('3D mesh is NOT implemented yet.')
      print('Program stopped.')
      exit()

    meshnode_dict   = { 'num_dim':num_dim,   \
                        'num_node':num_node, \
                        'coord_node':coord_node }
    meshelem_dict   = { 'num_elem':num_elem, \
                        'type_elem':type_elem, \
                        'dim_elem':dim_elem, \
                        'tag_elem':tag_elem, \
                        'elem2node_dict':elem2node_dict, \
                        'num_elembytype':num_elembytype, \
                        'elem_index_l2g':index_l2g_tmp,  \
                        'num_type_elem':self.num_type_elem, \
                        'num_nodebytype':self.num_nodebytype }
    facecell_list   = [num_face_inner, num_face_boundary, num_cell, face2node_inner, face2node_boundary, face2cell_inner, face2cell_boundary, face2tag_boundary, tag2name_physics_dict, cell2dim, cell2tag]
    celltmp_list    = [volume_cell, coord_cellcenter]

    # Clear all the model data
    gmsh.clear()
    
    # Finalize
    gmsh.finalize()

    return meshnode_dict, meshelem_dict, facecell_list, celltmp_list


  def set_virtualcell_boundary(self, config, facecell_list):

    # Virtual cells identification on boundary

    print('Setting virtual-cell on boundary...')

    num_face_boundary  = facecell_list[1]
    face2cell_boundary = facecell_list[6]
    face2tag_boundary  = facecell_list[7]
    tag2name_physics   = facecell_list[8]

    # Freestream, symmetric, outlet (gradient-free): virtual (virtualcell_bd=1)
    # Wall: non-virtual (virtualcell_bd=0)
    virtualcell_boundary = np.zeros(num_face_boundary).reshape(num_face_boundary).astype(int)
    for n_face in range(0,num_face_boundary):
      bd_name  = tag2name_physics[ face2cell_boundary[1,n_face] ]
      bd_kind  = config['boundary_conditions'][bd_name]['kind']
      if bd_kind == self.KIND_BOUNDARY_FREESTREAM :
        virtualcell_boundary[n_face] = 1
      elif bd_kind == self.KIND_BOUNDARY_TOTAL_PRESS_TEMP :
        virtualcell_boundary[n_face] = 1
      elif bd_kind == self.KIND_BOUNDARY_MASSFLOWRATE :
        virtualcell_boundary[n_face] = 1
      elif bd_kind == self.KIND_BOUNDARY_SYMMETRY :
        virtualcell_boundary[n_face] = 1
      elif bd_kind == self.KIND_BOUNDARY_OUTLET :
        virtualcell_boundary[n_face] = 1
      elif bd_kind  == self.KIND_BOUNDARY_OUTLET_PRESS_FIXED :
        virtualcell_boundary[n_face] = 1
      elif bd_kind == self.KIND_BOUNDARY_WALL_FIX :
        virtualcell_boundary[n_face] = 0
      elif bd_kind == self.KIND_BOUNDARY_WALL_SLIP :
        virtualcell_boundary[n_face] = 1
      elif bd_kind == self.KIND_BOUNDARY_AMBIENT :
        virtualcell_boundary[n_face] = 1
      else:
        print( 'No boundary ID', 'N_Face:',n_face)
        print( 'Check boundary condition. Program stopped')
        exit()

    return virtualcell_boundary


  def set_geometry(self, facecell_list, virtualcell_boundary):

    print('Setting geometry data...')

    num_face_inner     = facecell_list[0]
    num_face_boundary  = facecell_list[1]
    num_cell           = facecell_list[2]
    face2node_inner    = facecell_list[3]
    face2node_boundary = facecell_list[4]
    face2cell_inner    = facecell_list[5]
    face2cell_boundary = facecell_list[6]
    face2tag_boundary  = facecell_list[7]
    tag2name_physics   = facecell_list[8]
    cell2dim           = facecell_list[9]
    cell2tag           = facecell_list[10]

    geom_dict = { 'num_face_inner':num_face_inner, \
                  'num_face_boundary':num_face_boundary, \
                  'num_cell':num_cell, \
                  'face2node_inner':face2node_inner, \
                  'face2node_boundary':face2node_boundary, \
                  'face2cell_inner':face2cell_inner, \
                  'face2cell_boundary':face2cell_boundary, \
                  'cell2tag':cell2tag, \
                  'virtualcell_boundary':virtualcell_boundary,  \
                  'face2tag_boundary':face2tag_boundary, \
                  'tag2name_physics':tag2name_physics }

    return geom_dict


  def set_metrics(self, meshnode_dict, meshelem_dict, celltmp_list, geom_dict):

    # Metrics
    # Calculate area vectors
    print( 'Setting metrics...' )

    num_dim    = meshnode_dict['num_dim']
    num_node   = meshnode_dict['num_node']
    coord_node = meshnode_dict['coord_node']

    num_face_inner     = geom_dict['num_face_inner']
    num_face_boundary  = geom_dict['num_face_boundary']
    num_cell           = geom_dict['num_cell']
    face2node_inner    = geom_dict['face2node_inner']
    face2node_boundary = geom_dict['face2node_boundary']
    face2cell_inner    = geom_dict['face2cell_inner']
    face2cell_boundary = geom_dict['face2cell_boundary']

    volume_cell        = celltmp_list[0]
    coord_cellcenter   = celltmp_list[1]

    dz = 1.0

    print('--Calculate area and normal vectors...')

    # Area vector
    area_vec          = np.zeros(10*num_face_inner).reshape(10,num_face_inner)
    area_vec_boundary = np.zeros(10*num_face_boundary).reshape(10,num_face_boundary)

    #　内向きを正としている（Nodeは右回りが正）
    # --Inner faces
    for n_face in range(0,num_face_inner):
      im = face2node_inner[0,n_face]
      ip = face2node_inner[1,n_face]
      dx = coord_node[ip,0] - coord_node[im,0]
      dy = coord_node[ip,1] - coord_node[im,1]
      dl = np.sqrt( dx**2 + dy**2 )
      area_vec[0,n_face] = dl*dz
      area_vec[1,n_face] = dy/dl
      area_vec[2,n_face] =-dx/dl
      area_vec[3,n_face] = 0.0
      area_vec[4,n_face] = dx/dl
      area_vec[5,n_face] = dy/dl
      area_vec[6,n_face] = 0.0
      area_vec[7,n_face] = 0.0
      area_vec[8,n_face] = 0.0
      area_vec[9,n_face] = 0.0

    # --Boundary faces
    for n_face in range(0,num_face_boundary):
      im = face2node_boundary[0,n_face]
      ip = face2node_boundary[1,n_face]
      dx = coord_node[ip,0] - coord_node[im,0]
      dy = coord_node[ip,1] - coord_node[im,1]
      dl = np.sqrt( dx**2 + dy**2 )
      area_vec_boundary[0,n_face] = dl*dz
      area_vec_boundary[1,n_face] = dy/dl
      area_vec_boundary[2,n_face] =-dx/dl
      area_vec_boundary[3,n_face] = 0.0
      area_vec_boundary[4,n_face] = dx/dl
      area_vec_boundary[5,n_face] = dy/dl
      area_vec_boundary[6,n_face] = 0.0
      area_vec_boundary[7,n_face] = 0.0
      area_vec_boundary[8,n_face] = 0.0
      area_vec_boundary[9,n_face] = 0.0
      #print(n_face, area_vec_boundary[0,n_face], area_vec_boundary[1,n_face], area_vec_boundary[2,n_face])


    # Lengths between face center and cell center
    print('--Calculate lengths...')
    # lenght[0,:]-->自分自身を構成するセル側への距離
    # lenght[1,:]-->隣のセル側への距離
    length          = np.zeros(2*num_face_inner).reshape(2,num_face_inner)
    # lenght_bd[:]-->自分自身を構成するセル側への距離
    length_boundary = np.zeros(num_face_boundary).reshape(num_face_boundary)
    # --Inner faces
    for n_face in range(0,num_face_inner):
      # Face center
      im = face2node_inner[0,n_face]
      ip = face2node_inner[1,n_face]
      coord_face = 0.50*( coord_node[im,:]+coord_node[ip,:] )

      # Cell center
      n_cell     = face2cell_inner[0,n_face]
      coord_cell = coord_cellcenter[:,n_cell]
      # Lenght between cell center and face center
      length[0,n_face] = np.linalg.norm(coord_cell-coord_face)

      # Counter cell center
      n_cell     = face2cell_inner[1,n_face]
      coord_cell = coord_cellcenter[:,n_cell]
      # Lenght between cell center and face center
      length[1,n_face] = np.linalg.norm(coord_cell-coord_face)

    # --Boundary faces
    for n_face in range(0,num_face_boundary):
      # Face center
      im = face2node_boundary[0,n_face]
      ip = face2node_boundary[1,n_face]
      coord_face = 0.50*( coord_node[im,:]+coord_node[ip,:] )

      # Cell center
      n_cell     = face2cell_boundary[0,n_face]
      coord_cell = coord_cellcenter[:,n_cell]
      # Lenght between cell center and face center
      length_boundary[n_face] = np.linalg.norm(coord_cell-coord_face)
      
    metrics_dict = { 'area_vec_inner':area_vec, \
                     'area_vec_boundary':area_vec_boundary, \
                     'length_inner':length, \
                     'length_boundary':length_boundary, \
                     'volume_cell':volume_cell, \
                     'coord_cellcenter': coord_cellcenter }

    return metrics_dict


  def domain_partition(self, config):
    return