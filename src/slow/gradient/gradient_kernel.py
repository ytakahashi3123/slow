#!/usr/bin/env python3

# Program to accumulate the Green-Gauss gradient over faces and cells

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/10

from slow.general.jit import kernel


@kernel
def accumulate_gradient(num_face, num_face_bd, num_cell, num_primitiv,
                        face2cell, face2cell_bd, virtualcell_bd,
                        area_vec, area_vec_bd, length, length_bd, volume,
                        var_primitiv, var_primitiv_bd, var_gradient):
  """
  Green-Gauss gradient: var_gradient = (1/V) * sum over faces of ( var_face * n_out * S )

  面ループの本体。numba は dict を受け取れないので、値の取り出しは呼び出し側で行い、
  ここは numpy 配列とスカラーだけを見る。

  ループの形は元のままで、原始変数についての内側ループを明示的に書いてある
  （6 要素の numpy スライスはスカラー演算より遅く、numba も掛けられないため）。

  面の法線 area_vec[1:4] はセル a に対して内向きなので、a には引き b には足す。
  """

  var_gradient[:,:,:] = 0.0

  # --Inner loop
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
    for m in range(0,num_primitiv):
      var_face = fact_p*var_primitiv[m,n_cell_self] + fact_m*var_primitiv[m,n_cell_neig]
      var_gradient[0,m,n_cell_self] = var_gradient[0,m,n_cell_self] - var_face*vec_x
      var_gradient[1,m,n_cell_self] = var_gradient[1,m,n_cell_self] - var_face*vec_y
      var_gradient[0,m,n_cell_neig] = var_gradient[0,m,n_cell_neig] + var_face*vec_x
      var_gradient[1,m,n_cell_neig] = var_gradient[1,m,n_cell_neig] + var_face*vec_y

  # --Boundary loop
  for n_face in range(0,num_face_bd):
    n_cell_self = face2cell_bd[0,n_face]
    area        = area_vec_bd[0,n_face]
    vec_x       = area*area_vec_bd[1,n_face]
    vec_y       = area*area_vec_bd[2,n_face]
    # 境界面なので境界面用の距離を使う（内部面用の length は
    #  2 x num_face_inner なので num_face_boundary > num_face_inner で添字が溢れる）
    # なお重みは dl_n = dl_s*vcell なので dl_s が約分され、仮想セルの有無だけで決まる
    dl_s   = length_bd[n_face]
    dl_n   = length_bd[n_face]*float(virtualcell_bd[n_face])
    fact_m = dl_s/(dl_s+dl_n)
    fact_p = dl_n/(dl_s+dl_n)
    for m in range(0,num_primitiv):
      var_face = fact_p*var_primitiv[m,n_cell_self] + fact_m*var_primitiv_bd[m,n_face]
      var_gradient[0,m,n_cell_self] = var_gradient[0,m,n_cell_self] - var_face*vec_x
      var_gradient[1,m,n_cell_self] = var_gradient[1,m,n_cell_self] - var_face*vec_y

  # Divide by the cell volume
  for n_cell in range(0,num_cell):
    inv_volume = 1.0/volume[n_cell]
    for m in range(0,num_primitiv):
      var_gradient[0,m,n_cell] = var_gradient[0,m,n_cell]*inv_volume
      var_gradient[1,m,n_cell] = var_gradient[1,m,n_cell]*inv_volume

  return var_gradient


@kernel
def find_neighbour_maxmin(num_face, num_face_bd, num_primitiv,
                          face2cell, face2cell_bd,
                          var_primitiv, var_primitiv_bd, var_neig_maxmin):
  """
  Maximum and minimum of the primitive variables over the neighbouring cells

  var_neig_maxmin[0] に最大値、[1] に最小値を入れる。自セルの値も含める
  （制限関数の分子 var_max - var_self が負にならないようにするため）。

  面ループの本体。max / min は 2 個の比較なので、元の numpy 版と厳密に同じ値になる。
  """

  var_neig_maxmin[0,:,:] = var_primitiv[:,:]
  var_neig_maxmin[1,:,:] = var_primitiv[:,:]

  for n_face in range(0,num_face):
    n_cell_self = face2cell[0,n_face]
    n_cell_neig = face2cell[1,n_face]
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

  return var_neig_maxmin


@kernel
def apply_minmod_limiter(num_face, num_face_bd, num_primitiv,
                         face2cell, face2cell_bd, virtualcell_bd,
                         area_vec, area_vec_bd, length, length_bd,
                         var_primitiv, var_gradient, var_neig_maxmin, var_limiter):
  """
  Minmod slope limiter: phi = min over faces of clip( (var_max - var_self)/grad_face, 0, 1 )

  面ループの本体。セル中心から面までの外挿量 grad_face に対して、隣接セルの最大／最小を
  超えない範囲に収まる係数を求め、そのセルの全ての面について最小値を採る。

  grad_face に 1e-20 を足す（引く）のは 0 除算を避けるためで、符号は保つ。
  """

  var_limiter[:,:] = 1.0

  # Inner faces
  for n_face in range(0,num_face):
    n_cell_self = face2cell[0,n_face]
    n_cell_neig = face2cell[1,n_face]
    dl_s   = length[0,n_face]
    dl_n   = length[1,n_face]
    vec_x  = area_vec[1,n_face]*(dl_s+dl_n)
    vec_y  = area_vec[2,n_face]*(dl_s+dl_n)
    vec_z  = area_vec[3,n_face]*(dl_s+dl_n)

    for m in range(0,num_primitiv):
      # --selfside cell
      grad_face =-(var_gradient[0,m,n_cell_self]*vec_x+var_gradient[1,m,n_cell_self]*vec_y+var_gradient[2,m,n_cell_self]*vec_z)
      if grad_face >= 0.0:
        grad_face = grad_face + 1.e-20
        del_face  = var_neig_maxmin[0,m,n_cell_self] - var_primitiv[m,n_cell_self]
      else :
        grad_face = grad_face - 1.e-20
        del_face  = var_neig_maxmin[1,m,n_cell_self] - var_primitiv[m,n_cell_self]
      del_face  = max( 0.0, min(1.0, del_face/grad_face) )
      var_limiter[m,n_cell_self] = min(var_limiter[m,n_cell_self], del_face)

      # --neigboring side cell
      grad_face = (var_gradient[0,m,n_cell_neig]*vec_x+var_gradient[1,m,n_cell_neig]*vec_y+var_gradient[2,m,n_cell_neig]*vec_z)
      if grad_face >= 0.0:
        grad_face = grad_face + 1.e-20
        del_face  = var_neig_maxmin[0,m,n_cell_neig] - var_primitiv[m,n_cell_neig]
      else :
        grad_face = grad_face - 1.e-20
        del_face  = var_neig_maxmin[1,m,n_cell_neig] - var_primitiv[m,n_cell_neig]
      del_face  = max( 0.0, min(1.0, del_face/grad_face) )
      var_limiter[m,n_cell_neig] = min(var_limiter[m,n_cell_neig], del_face)

  # Boundary faces
  for n_face in range(0,num_face_bd):
    n_cell_self = face2cell_bd[0,n_face]
    dl_s   = length_bd[n_face]
    dl_n   = length_bd[n_face]*float(virtualcell_bd[n_face])
    vec_x  = area_vec_bd[1,n_face]*(dl_s+dl_n)
    vec_y  = area_vec_bd[2,n_face]*(dl_s+dl_n)
    vec_z  = area_vec_bd[3,n_face]*(dl_s+dl_n)

    for m in range(0,num_primitiv):
      grad_face =-(var_gradient[0,m,n_cell_self]*vec_x+var_gradient[1,m,n_cell_self]*vec_y+var_gradient[2,m,n_cell_self]*vec_z)
      if grad_face >= 0.0:
        grad_face = grad_face + 1.e-20
        del_face  = var_neig_maxmin[0,m,n_cell_self] - var_primitiv[m,n_cell_self]
      else :
        grad_face = grad_face - 1.e-20
        del_face  = var_neig_maxmin[1,m,n_cell_self] - var_primitiv[m,n_cell_self]
      del_face  = max( 0.0, min(1.0, del_face/grad_face) )
      var_limiter[m,n_cell_self] = min(var_limiter[m,n_cell_self], del_face)

  return var_limiter
