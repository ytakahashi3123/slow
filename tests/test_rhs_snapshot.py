#!/usr/bin/env python3

# Snapshot regression test for the residual (RHS) evaluation
#
# メッシュ読み込みから境界条件・輸送係数・空間勾配・数値流束までを実際に通し、
# 得られた残差を tests/data/rhs_snapshot.npz の基準値と突き合わせる。
# 陰解法 (LU-SGS) は残差の計算に関与しないので、このテストは時間積分側の
# 変更からは独立している。リファクタリングで挙動が変わっていないことの確認に使う。
#
# 基準値を作り直すとき（残差の値を意図的に変えたとき）だけ、以下を実行する:
#   .venv/bin/python tests/test_rhs_snapshot.py
# 実行すると tests/data/rhs_snapshot.npz を上書きするので、
# 差分をレビューしたうえでコミットすること。

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

import pathlib

import numpy as np
import pytest

from slow.boundary.boundary import boundary
from slow.flowfield.flowfield import flowfield
from slow.gradient.gradient import gradient
from slow.meshdata.meshdata import meshdata
from slow.orbital.orbital import orbital
from slow.rhs.rhs import rhs


DATA_DIR      = pathlib.Path(__file__).parent/'data'
FILE_CONFIG   = DATA_DIR/'config_rhs_snapshot.yml'
FILE_MESH     = DATA_DIR/'chimera.msh'
FILE_SNAPSHOT = DATA_DIR/'rhs_snapshot.npz'

# 検査する移流スキーム
ADVECTION_SCHEMES = ['slau2', 'haenel']


def set_perturbed_flowfield(var_primitiv, var_conserv, metrics_dict, gas_property_dict, get_total_energy):
  """
  一様な初期条件では勾配も制限関数もほとんど働かないため、セル中心座標の
  滑らかな関数で原始変数に摂動を与える。全ての項に非自明な値を通すのが目的。

  メッシュ寸法に依存しないよう、座標はバウンディングボックスで正規化する。
  """

  gas_constant        = gas_property_dict['gas_constant']
  specific_heat_volum = gas_property_dict['specific_heat_volume']

  coord = metrics_dict['coord_cellcenter']
  span  = np.where( coord.max(axis=1) - coord.min(axis=1) > 0.0, \
                    coord.max(axis=1) - coord.min(axis=1), 1.0 )
  xhat  = (coord[0,:] - coord[0,:].min())/span[0]
  yhat  = (coord[1,:] - coord[1,:].min())/span[1]

  dens0 = var_primitiv[0,:].copy()
  uvel0 = var_primitiv[1,:].copy()
  vvel0 = var_primitiv[2,:].copy()
  temp0 = var_primitiv[4,:].copy()

  # 密度と温度は正のまま保つよう振幅を 1 未満に抑える
  var_primitiv[0,:] = dens0*( 1.0 + 0.30*np.sin(2.0*np.pi*xhat)*np.cos(3.0*np.pi*yhat) )
  var_primitiv[1,:] = uvel0*( 1.0 + 0.40*np.sin(1.0*np.pi*xhat + 2.0*np.pi*yhat) )
  var_primitiv[2,:] = uvel0*0.25*np.cos(2.0*np.pi*xhat)*np.sin(1.0*np.pi*yhat) + vvel0
  var_primitiv[3,:] = 0.0
  var_primitiv[4,:] = temp0*( 1.0 + 0.20*np.cos(3.0*np.pi*xhat)*np.sin(2.0*np.pi*yhat) )
  var_primitiv[5,:] = var_primitiv[0,:]*gas_constant*var_primitiv[4,:]

  var_conserv[0,:] = var_primitiv[0,:]
  var_conserv[1,:] = var_primitiv[0,:]*var_primitiv[1,:]
  var_conserv[2,:] = var_primitiv[0,:]*var_primitiv[2,:]
  var_conserv[3,:] = var_primitiv[0,:]*var_primitiv[3,:]
  for n_cell in range(0, var_conserv.shape[1]):
    var_conserv[4,n_cell] = get_total_energy(var_primitiv[0,n_cell], specific_heat_volum, \
                                             var_primitiv[4,n_cell], var_primitiv[1:4,n_cell])

  return var_primitiv, var_conserv


def evaluate_rhs(advection_scheme):
  """メッシュ読み込みから残差の評価までを一度だけ通し、途中の変数もまとめて返す"""

  orbital_obj   = orbital()
  meshdata_obj  = meshdata()
  flowfield_obj = flowfield()
  boundary_obj  = boundary()
  gradient_obj  = gradient()
  rhs_obj       = rhs()

  config = orbital_obj.read_config_yaml(str(FILE_CONFIG))
  # メッシュはテストデータの絶対パスで与える（実行時のカレントディレクトリに依存させない）
  config['meshdata_io']['filename_mesh']            = str(FILE_MESH)
  config['numericalflux_setting']['advection_scheme'] = advection_scheme

  dimension_dict = orbital_obj.set_dimension()

  meshnode_dict, meshelem_dict, geom_dict, metrics_dict = meshdata_obj.set_mesh_routine(config)

  gas_property_dict = flowfield_obj.set_gas_properties(config)

  transport_coefficient_dict, \
  var_primitiv, var_primitiv_bd, \
  var_conserv, var_conserv_prev = flowfield_obj.define_variables(config, dimension_dict, geom_dict)

  var_primitiv, var_conserv, var_conserv_prev, _ = \
    flowfield_obj.initialize_flowfield(config, dimension_dict, geom_dict, metrics_dict, \
                                       meshnode_dict, meshelem_dict, gas_property_dict, \
                                       var_primitiv, var_conserv, var_conserv_prev)

  var_primitiv, var_conserv = set_perturbed_flowfield(var_primitiv, var_conserv, metrics_dict, \
                                                      gas_property_dict, orbital_obj.get_total_energy)

  var_gradient, var_limiter, var_neig_maxmin = gradient_obj.initialize_gradient(config, dimension_dict, geom_dict)
  var_rhs = rhs_obj.initialize_rhs(config, dimension_dict, geom_dict)

  var_primitiv_bd = boundary_obj.boundary_condition(config, geom_dict, metrics_dict, \
                                                    gas_property_dict, var_primitiv, var_conserv, \
                                                    var_primitiv_bd)

  transport_coefficient_dict = flowfield_obj.set_transport_coefficients(config, geom_dict, \
                                                                        gas_property_dict, var_primitiv, var_primitiv_bd, \
                                                                        transport_coefficient_dict)

  var_gradient, var_limiter = gradient_obj.gradient_routine(config, dimension_dict, geom_dict, \
                                                            metrics_dict, var_primitiv, var_primitiv_bd, \
                                                            var_gradient, var_neig_maxmin, var_limiter)

  var_rhs = rhs_obj.rhs_routine(config, dimension_dict, geom_dict, metrics_dict, \
                                gas_property_dict, transport_coefficient_dict, \
                                var_primitiv, var_primitiv_bd, var_gradient, var_limiter, \
                                var_rhs)

  # 段階ごとに比較できるよう、中間結果も併せて返す
  return { 'area_vec_inner': metrics_dict['area_vec_inner'],       \
           'volume_cell':    metrics_dict['volume_cell'],          \
           'viscosity':      transport_coefficient_dict['viscosity'], \
           'primitiv_bd':    var_primitiv_bd,                      \
           'gradient':       var_gradient,                         \
           'limiter':        var_limiter,                          \
           'rhs':            var_rhs }


@pytest.fixture(scope='module')
def snapshot():
  if not FILE_SNAPSHOT.exists():
    pytest.fail(f'基準値が見つからない: {FILE_SNAPSHOT}\n'
                f'  作り直すには: python tests/test_rhs_snapshot.py')
  return np.load(FILE_SNAPSHOT)


@pytest.mark.parametrize('advection_scheme', ADVECTION_SCHEMES)
def test_rhs_matches_snapshot(snapshot, advection_scheme):
  # メッシュから残差までの計算結果が基準値と一致すること。
  # リファクタリングでは完全一致するはずなので許容差は倍精度の丸め程度に取る

  result = evaluate_rhs(advection_scheme)

  for name, value in result.items():
    key = f'{advection_scheme}__{name}'
    assert key in snapshot, f'基準値に {key} がない。基準値の作り直しが必要'
    np.testing.assert_allclose(value, snapshot[key], rtol=1.0e-12, \
                               atol=1.0e-12*np.abs(snapshot[key]).max(), \
                               err_msg=f'{name} ({advection_scheme}) が基準値と一致しない')


def test_snapshot_exercises_every_term(snapshot):
  # 基準値が「ほぼゼロの一様場」になっていないことの確認。
  # 摂動を与え損なうとテストが何も検査しなくなるため、それを防ぐ

  for advection_scheme in ADVECTION_SCHEMES:
    for name in ('gradient', 'rhs'):
      value = snapshot[f'{advection_scheme}__{name}']
      assert np.abs(value).max() > 0.0, f'{name} ({advection_scheme}) が全てゼロ'

    # 制限関数が 1.0 のまま（=どこも制限されていない）だと MUSCL の検査にならない
    limiter = snapshot[f'{advection_scheme}__limiter']
    assert limiter.min() < 1.0, f'制限関数 ({advection_scheme}) がどこにも効いていない'


def main():
  # 基準値の作り直し。意図的に残差を変えたときだけ実行する

  snapshot = {}
  for advection_scheme in ADVECTION_SCHEMES:
    for name, value in evaluate_rhs(advection_scheme).items():
      snapshot[f'{advection_scheme}__{name}'] = value

  np.savez(FILE_SNAPSHOT, **snapshot)
  print(f'Writing snapshot data : {FILE_SNAPSHOT}')
  for key in sorted(snapshot):
    print(f'  {key:32s} shape={snapshot[key].shape}')

  return


if __name__ == '__main__':

  main()
