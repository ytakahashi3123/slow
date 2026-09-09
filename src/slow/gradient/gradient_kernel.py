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
