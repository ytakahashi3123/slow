#!/usr/bin/env python3

# Program to verify the boundary contribution of the Green-Gauss gradient
#
# get_gradient の境界ループは内部面用の length を境界面番号で引いていた。
# 重みは dl_n = dl_s*vcell なので dl_s が約分され結果は正しかったが、
# length は 2 x num_face_inner なので num_face_boundary > num_face_inner の
# 格子で IndexError になる。チュートリアルの格子はどれも内部面のほうが多く、
# 既存のスナップショット回帰では踏めなかった。
#
# ここでは境界面のほうが多い最小の格子（正方形セル 2 個）を組み、
# 境界ループが動くことと、線形場に対して勾配が厳密に再現されることを固定する。

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

import numpy as np
import pytest

from slow.gradient.gradient import gradient


NUM_PRIMITIV = 6
NUM_SPATIAL  = 3

# 線形場 f_m(x,y) = OFFSET[m] + SLOPE_X[m]*x + SLOPE_Y[m]*y
OFFSET  = np.array([ 1.0e-4,  600.0,  120.0,  0.0,  300.0,  8.6e0])
SLOPE_X = np.array([ 2.0e-5,  -30.0,   11.0,  0.0,   17.0, -1.3e0])
SLOPE_Y = np.array([-1.0e-5,   45.0,   -7.0,  0.0,  -23.0,  0.7e0])


def linear_field(coord):
  # 与えた座標での線形場の値（原始変数の本数だけ並べる）

  return OFFSET + SLOPE_X*coord[0] + SLOPE_Y*coord[1]


def make_two_cell_mesh(virtualcell_boundary):
  """
  一辺 1 の正方形セル 2 個（[0,1]x[0,1] と [1,2]x[0,1]、奥行き 1）。

  内部面 1 枚に対して境界面 6 枚なので num_face_boundary > num_face_inner になる。
  法線の向きはソルバの規約に合わせる: area_vec[1:4] はセル a に対して内向き
  （内部面は face2cell[1] から face2cell[0] を向き、境界面もセル側を向く）。
  """

  coord_cell = np.array([[0.5, 1.5],
                         [0.5, 0.5],
                         [0.0, 0.0]])

  # 内部面: 共有辺 x=1。セル a=0 の内向き、すなわち -x 方向
  area_vec_inner = np.zeros((10, 1))
  area_vec_inner[0,0] = 1.0
  area_vec_inner[1,0] = -1.0
  length_inner = np.array([[0.5], [0.5]])
  face2cell_inner = np.array([[0], [1]])

  # 境界面: (セル番号, 面中心, セルから見た外向き法線)
  faces_bd = [ (0, (0.0, 0.5), (-1.0,  0.0)),
               (0, (0.5, 0.0), ( 0.0, -1.0)),
               (0, (0.5, 1.0), ( 0.0,  1.0)),
               (1, (1.5, 0.0), ( 0.0, -1.0)),
               (1, (1.5, 1.0), ( 0.0,  1.0)),
               (1, (2.0, 0.5), ( 1.0,  0.0)) ]

  num_face_bd  = len(faces_bd)
  area_vec_bd  = np.zeros((10, num_face_bd))
  face2cell_bd = np.zeros((2, num_face_bd), dtype=int)
  length_bd    = np.zeros(num_face_bd)
  coord_face_bd = np.zeros((3, num_face_bd))

  for n_face, (n_cell, coord_face, vec_out) in enumerate(faces_bd):
    area_vec_bd[0,n_face] = 1.0
    # 格納するのはセルに対する内向き法線なので外向きを反転する
    area_vec_bd[1,n_face] = -vec_out[0]
    area_vec_bd[2,n_face] = -vec_out[1]
    face2cell_bd[0,n_face] = n_cell
    coord_face_bd[0,n_face] = coord_face[0]
    coord_face_bd[1,n_face] = coord_face[1]
    length_bd[n_face] = np.linalg.norm( coord_face_bd[0:2,n_face] - coord_cell[0:2,n_cell] )

  geom_dict = { 'num_face_inner': 1,                                  \
                'num_face_boundary': num_face_bd,                     \
                'num_cell': 2,                                        \
                'face2cell_inner': face2cell_inner,                   \
                'face2cell_boundary': face2cell_bd,                   \
                'virtualcell_boundary': np.asarray(virtualcell_boundary) }

  metrics_dict = { 'area_vec_inner': area_vec_inner,       \
                   'area_vec_boundary': area_vec_bd,       \
                   'length_inner': length_inner,           \
                   'length_boundary': length_bd,           \
                   'volume_cell': np.array([1.0, 1.0]) }

  return geom_dict, metrics_dict, coord_cell, coord_face_bd


def make_state(geom_dict, coord_cell, coord_face_bd):
  """
  線形場をセル中心と境界面に設定する。

  境界面の値は、境界ループの重み（仮想セルの有無で決まる）を通したあとに
  面中心の値になるように選ぶ。こうすると Green-Gauss が線形場に対して厳密になる。
    --仮想セルあり (vcell=1): var_face = 0.5*(セル値 + 与えた値) なので、
      面を挟んでセル中心を鏡像にした位置の値を与える
    --仮想セルなし (vcell=0): var_face = 与えた値 なので、面中心の値を与える
  """

  virtualcell_bd = geom_dict['virtualcell_boundary']
  face2cell_bd   = geom_dict['face2cell_boundary']

  var_primitiv = np.zeros((NUM_PRIMITIV, geom_dict['num_cell']))
  for n_cell in range(0, geom_dict['num_cell']):
    var_primitiv[:,n_cell] = linear_field(coord_cell[:,n_cell])

  var_primitiv_bd = np.zeros((NUM_PRIMITIV, geom_dict['num_face_boundary']))
  for n_face in range(0, geom_dict['num_face_boundary']):
    coord_face = coord_face_bd[:,n_face]
    if virtualcell_bd[n_face] == 1:
      n_cell = face2cell_bd[0,n_face]
      coord_ghost = 2.0*coord_face - coord_cell[:,n_cell]
      var_primitiv_bd[:,n_face] = linear_field(coord_ghost)
    else:
      var_primitiv_bd[:,n_face] = linear_field(coord_face)

  return var_primitiv, var_primitiv_bd


def run_gradient(geom_dict, metrics_dict, var_primitiv, var_primitiv_bd):
  # get_gradient は config を参照しないので、辞書は dimension_dict だけ与えれば足りる

  gradient_obj = gradient()
  var_gradient = np.zeros((NUM_SPATIAL, NUM_PRIMITIV, geom_dict['num_cell']))

  return gradient_obj.get_gradient({}, {'num_primitive': NUM_PRIMITIV}, \
                                   geom_dict, metrics_dict, \
                                   var_primitiv, var_primitiv_bd, var_gradient)


VIRTUALCELL_CASES = [
  pytest.param([1, 1, 1, 1, 1, 1], id='virtual-all'),
  pytest.param([0, 0, 0, 0, 0, 0], id='virtual-none'),
  pytest.param([1, 0, 1, 0, 1, 0], id='virtual-mixed'),
]


def test_the_mesh_really_has_more_boundary_faces_than_inner_faces():
  # この検査の前提。ここが崩れると回帰テストの意味が無くなる

  geom_dict, _, _, _ = make_two_cell_mesh([1, 1, 1, 1, 1, 1])

  assert geom_dict['num_face_boundary'] > geom_dict['num_face_inner']


def test_the_mesh_is_geometrically_closed():
  # 各セルについて外向き法線の面積重み和が 0 になること（法線の向きの取り違え防止）

  geom_dict, metrics_dict, _, _ = make_two_cell_mesh([1, 1, 1, 1, 1, 1])

  closure = np.zeros((3, geom_dict['num_cell']))
  for n_face in range(0, geom_dict['num_face_inner']):
    n_a = geom_dict['face2cell_inner'][0,n_face]
    n_b = geom_dict['face2cell_inner'][1,n_face]
    vec = metrics_dict['area_vec_inner'][1:4,n_face]*metrics_dict['area_vec_inner'][0,n_face]
    closure[:,n_a] -= vec
    closure[:,n_b] += vec
  for n_face in range(0, geom_dict['num_face_boundary']):
    n_a = geom_dict['face2cell_boundary'][0,n_face]
    closure[:,n_a] -= metrics_dict['area_vec_boundary'][1:4,n_face] \
                      *metrics_dict['area_vec_boundary'][0,n_face]

  np.testing.assert_allclose(closure, 0.0, atol=1.0e-14)


@pytest.mark.parametrize('virtualcell_boundary', VIRTUALCELL_CASES)
def test_boundary_loop_handles_more_boundary_faces_than_inner_faces(virtualcell_boundary):
  # 回帰: 境界ループが内部面用の length を引いていると IndexError になる

  geom_dict, metrics_dict, coord_cell, coord_face_bd = make_two_cell_mesh(virtualcell_boundary)
  var_primitiv, var_primitiv_bd = make_state(geom_dict, coord_cell, coord_face_bd)

  var_gradient = run_gradient(geom_dict, metrics_dict, var_primitiv, var_primitiv_bd)

  assert np.all( np.isfinite(var_gradient) )


@pytest.mark.parametrize('virtualcell_boundary', VIRTUALCELL_CASES)
def test_gradient_is_exact_for_a_linear_field(virtualcell_boundary):
  # Green-Gauss は線形場に対して厳密。境界の寄与が狂うとここで落ちる

  geom_dict, metrics_dict, coord_cell, coord_face_bd = make_two_cell_mesh(virtualcell_boundary)
  var_primitiv, var_primitiv_bd = make_state(geom_dict, coord_cell, coord_face_bd)

  var_gradient = run_gradient(geom_dict, metrics_dict, var_primitiv, var_primitiv_bd)

  for n_cell in range(0, geom_dict['num_cell']):
    np.testing.assert_allclose(var_gradient[0,:,n_cell], SLOPE_X, rtol=1.0e-12, atol=1.0e-14)
    np.testing.assert_allclose(var_gradient[1,:,n_cell], SLOPE_Y, rtol=1.0e-12, atol=1.0e-14)
  # 2 次元実装なので z 成分は計算されない
  np.testing.assert_allclose(var_gradient[2,:,:], 0.0, atol=0.0)


def test_boundary_weights_depend_only_on_the_virtual_cell_flag():
  """
  境界の重みは dl_n = dl_s*vcell なので dl_s が約分され、距離の値には依らない。
  この性質のために内部面用の length を引いていても結果が正しく、誤りが露見しなかった。
  性質そのものを固定して、距離を使う形に「直した」ときに気付けるようにする。
  """

  virtualcell_boundary = [1, 0, 1, 0, 1, 0]
  geom_dict, metrics_dict, coord_cell, coord_face_bd = make_two_cell_mesh(virtualcell_boundary)
  var_primitiv, var_primitiv_bd = make_state(geom_dict, coord_cell, coord_face_bd)

  reference = run_gradient(geom_dict, metrics_dict, var_primitiv, var_primitiv_bd).copy()

  # 境界面の距離だけを何倍にしても結果は変わらない
  metrics_dict['length_boundary'] = metrics_dict['length_boundary']*7.3
  scaled = run_gradient(geom_dict, metrics_dict, var_primitiv, var_primitiv_bd)

  np.testing.assert_allclose(scaled, reference, rtol=0.0, atol=0.0)
