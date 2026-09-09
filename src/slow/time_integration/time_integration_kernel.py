#!/usr/bin/env python3

# Program to hold the face and cell loops of the time-integration routines

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/10

import numpy as np

from slow.general import thermodynamics
from slow.general.jit import kernel


@kernel
def accumulate_character_time(num_face, num_face_bd, num_cell, num_primitiv,
                              face2cell, face2cell_bd, virtualcell_bd,
                              area_vec, area_vec_bd, length, length_bd, volume,
                              specfic_heat_ratio, var_primitiv, var_primitiv_bd,
                              viscosity, viscosity_bd, character_time):
  """
  Characteristic time of each cell: V / max over faces of ( lambda*S )

  局所時間刻みのもとになる量。lambda=|u.n|+c+2*mu/(rho*d) を面上の平均状態で評価し、
  セルに接する面についての最大値を採る。

  2 次元実装なので法線の z 成分は 0 として扱う（元の実装のまま）。
  """

  character_time[:] = 0.0

  prim = np.zeros(num_primitiv)

  for n_face in range(0,num_face):
    area = area_vec[0,n_face]
    vecx = area_vec[1,n_face]
    vecy = area_vec[2,n_face]
    vecz = 0.0

    leng_a = length[0,n_face]
    leng_b = length[1,n_face]

    n_cell_a = face2cell[0,n_face]
    n_cell_b = face2cell[1,n_face]

    # Values on cell interface
    for m in range(0,num_primitiv):
      prim[m] = 0.50*( var_primitiv[m,n_cell_a] + var_primitiv[m,n_cell_b] )

    lenght_tmp = leng_a + leng_b
    visc_tmp   = 0.50*( viscosity[n_cell_a] + viscosity[n_cell_b] )
    eigenvalue = thermodynamics.get_max_eigenvalue(specfic_heat_ratio, prim[0], prim[1], \
                                                   prim[2], prim[3], prim[5],            \
                                                   visc_tmp, lenght_tmp, vecx, vecy, vecz)
    character_time[n_cell_a] = max(character_time[n_cell_a], eigenvalue*area)
    character_time[n_cell_b] = max(character_time[n_cell_b], eigenvalue*area)

  for n_face in range(0,num_face_bd):
    area = area_vec_bd[0,n_face]
    vecx = area_vec_bd[1,n_face]
    vecy = area_vec_bd[2,n_face]
    vecz = 0.0

    vcell_bd = virtualcell_bd[n_face]

    leng_a = length_bd[n_face]
    leng_b = float(vcell_bd)*leng_a

    n_cell_a = face2cell_bd[0,n_face]

    # 仮想セルがあれば平均、無ければ境界値をそのまま使う
    for m in range(0,num_primitiv):
      prim[m] = float(vcell_bd)*0.50*( var_primitiv[m,n_cell_a] + var_primitiv_bd[m,n_face] ) \
              + float(1-vcell_bd)*var_primitiv_bd[m,n_face]

    lenght_tmp = leng_a + leng_b
    visc_tmp   = float(vcell_bd)*0.50*( viscosity[n_cell_a] + viscosity_bd[n_face] ) \
               + float(1-vcell_bd)*viscosity_bd[n_face]
    eigenvalue = thermodynamics.get_max_eigenvalue(specfic_heat_ratio, prim[0], prim[1], \
                                                   prim[2], prim[3], prim[5],            \
                                                   visc_tmp, lenght_tmp, vecx, vecy, vecz)
    character_time[n_cell_a] = max(character_time[n_cell_a], eigenvalue*area)

  for n_cell in range(0,num_cell):
    character_time[n_cell] = volume[n_cell]/character_time[n_cell]

  return character_time


@kernel
def update_primitive_from_conservative(num_cell, gas_constant, specific_heat_volum,
                                       var_conserv, var_primitiv):
  """
  Primitive variables from the conservative ones

  Q=(rho, rho*u, rho*v, rho*w, E) から (rho, u, v, w, T, p) を作る。
  T = ( E - 0.5*rho*U^2 )/(rho*Cv)、p = rho*R*T。

  負の密度・温度・圧力の検査は呼び出し側で行う（そこでログを出して止めるため）。
  """

  for n_cell in range(0,num_cell):
    dens = var_conserv[0,n_cell]
    uvel = var_conserv[1,n_cell]/dens
    vvel = var_conserv[2,n_cell]/dens
    wvel = var_conserv[3,n_cell]/dens
    temp = ( var_conserv[4,n_cell]                                           \
           - 0.50*dens*(uvel*uvel + vvel*vvel + wvel*wvel) )/(dens*specific_heat_volum)

    var_primitiv[0,n_cell] = dens
    var_primitiv[1,n_cell] = uvel
    var_primitiv[2,n_cell] = vvel
    var_primitiv[3,n_cell] = wvel
    var_primitiv[4,n_cell] = temp
    var_primitiv[5,n_cell] = dens*gas_constant*temp

  return var_primitiv
