#!/usr/bin/env python3

# Program to verify the LU-SGS sweep against an explicitly assembled implicit operator
#
# lusgs_diagonal / lusgs_sweep は D+L+U を面ごとに組み立てて行列を持たないため、
# セル番号・法線の向き・掃引順の取り違えが黙って「近似解法」に化ける。
# ここでは小さな人工メッシュに対して D, L, U を密行列として独立に組み立て、
# スイープの結果が (D+L)*D^-1*(D+U)*dq = b を厳密に満たすことを検査する。
# 非対角ブロックの流束ヤコビアンは複素ステップ微分で作るので、
# flux_jacobian.py の実装とは独立になっている。

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

import pathlib

import numpy as np
import pytest

from slow.time_integration import lusgs_diagonal
from slow.time_integration import lusgs_sweep
from slow.time_integration import roe_dissipation


NUM_CONSERV = 5

SPECIFIC_HEAT_RATIO = 1.40
GAS_CONSTANT        = 8.3144598/28.8e-3
SPECIFIC_HEAT_VOLUM = GAS_CONSTANT/(SPECIFIC_HEAT_RATIO-1.0)

TIMESTEP_OUTER = 2.0e-6


# ---------------------------------------------------------------- 参照実装（独立）

def normal_flux(var_conserv, vec):
  # 保存変数 Q=(rho, rho*u, rho*v, rho*w, E) に対する法線流束 F(Q).n

  dens, momx, momy, momz, energy = var_conserv
  uvel, vvel, wvel = momx/dens, momy/dens, momz/dens

  cvel = uvel*vec[0] + vvel*vec[1] + wvel*vec[2]
  pres = (SPECIFIC_HEAT_RATIO-1.0)*( energy - 0.50*dens*(uvel**2 + vvel**2 + wvel**2) )

  return np.array([ dens*cvel,                    \
                    momx*cvel + vec[0]*pres,      \
                    momy*cvel + vec[1]*pres,      \
                    momz*cvel + vec[2]*pres,      \
                    (energy + pres)*cvel ])


def jacobian_by_complex_step(var_conserv, vec):
  # 複素ステップ微分による厳密な A(n)=d(F.n)/dQ

  step     = 1.0e-20
  jacobian = np.zeros((NUM_CONSERV, NUM_CONSERV))
  for n in range(0, NUM_CONSERV):
    var_perturbed    = np.array(var_conserv, dtype=complex)
    var_perturbed[n] = var_perturbed[n] + step*1.0j
    jacobian[:,n]    = np.imag( normal_flux(var_perturbed, vec) )/step

  return jacobian


def conservative_from_primitive(prim):
  # (rho,u,v,w,T,p) -> (rho, rho*u, rho*v, rho*w, E)

  dens = prim[0]
  vel  = prim[1:4]
  energy = dens*SPECIFIC_HEAT_VOLUM*prim[4] + 0.50*dens*np.dot(vel, vel)

  return np.array([dens, dens*vel[0], dens*vel[1], dens*vel[2], energy])


def max_eigenvalue(prim, viscosity, length, vec):
  # lambda = |u.n| + c + 2*mu/(rho*d)

  cvel = prim[1]*vec[0] + prim[2]*vec[1] + prim[3]*vec[2]
  sos  = np.sqrt( SPECIFIC_HEAT_RATIO*prim[5]/prim[0] )

  return abs(cvel) + sos + 2.0*viscosity/(prim[0]*length)


def absolute_jacobian(prim_a, prim_b, vec):
  # |A_Roe|。roe_dissipation は test_roe_dissipation.py で P|Lambda|P^-1 == A が固定済み

  absjacobian = np.zeros((NUM_CONSERV, NUM_CONSERV))
  pmatrix     = np.zeros((NUM_CONSERV, NUM_CONSERV))
  roe_dissipation.get_face_dissipation(SPECIFIC_HEAT_RATIO, SPECIFIC_HEAT_VOLUM, \
                                       prim_a, prim_b, vec[0], vec[1], vec[2],   \
                                       absjacobian, pmatrix)

  return absjacobian


# ---------------------------------------------------------------- 人工メッシュ

def make_case(kind_steady_mode='steady', lusgs_beta=1.0,
              kind_backward_difference='2nd_backward_diff'):
  """
  4 セル・内部面 5 枚・境界面 3 枚の人工メッシュ。
  内部面は掃引が依存する不変条件を満たすように並べる:
    face2cell_inner[0,n] < face2cell_inner[1,n] かつ第 0 行が単調非減少
  """

  num_cell    = 4
  face2cell   = np.array([[0, 0, 1, 1, 2],
                          [1, 2, 2, 3, 3]])
  num_face    = face2cell.shape[1]

  face2cell_bd   = np.array([[0, 1, 3],
                             [1, 2, 3]])
  virtualcell_bd = np.array([1, 0, 1])
  num_face_bd    = face2cell_bd.shape[1]

  # 面法線は単位ベクトル。2 次元実装なので z 成分は 0
  angles   = np.array([0.3, 1.1, 2.0, 2.9, 4.1])
  area_vec = np.zeros((10, num_face))
  area_vec[0,:] = np.array([1.3e-3, 0.9e-3, 1.7e-3, 1.1e-3, 2.1e-3])
  area_vec[1,:] = np.cos(angles)
  area_vec[2,:] = np.sin(angles)
  area_vec[3,:] = 0.0

  angles_bd   = np.array([0.7, 2.4, 5.0])
  area_vec_bd = np.zeros((10, num_face_bd))
  area_vec_bd[0,:] = np.array([1.5e-3, 1.2e-3, 0.8e-3])
  area_vec_bd[1,:] = np.cos(angles_bd)
  area_vec_bd[2,:] = np.sin(angles_bd)
  area_vec_bd[3,:] = 0.0

  length    = np.array([[1.1e-3, 0.8e-3, 1.4e-3, 0.9e-3, 1.6e-3],
                        [1.3e-3, 1.0e-3, 1.2e-3, 1.5e-3, 0.7e-3]])
  length_bd = np.array([1.0e-3, 1.2e-3, 0.9e-3])
  volume    = np.array([2.1e-9, 3.4e-9, 1.7e-9, 2.8e-9])

  # 一様場では非対角が効かないので、セルごとに十分ばらつかせる
  var_primitiv = np.zeros((6, num_cell))
  var_primitiv[0,:] = np.array([1.0e-4, 1.4e-4, 0.7e-4, 1.9e-4])
  var_primitiv[1,:] = np.array([ 600.0,  520.0,  680.0,  450.0])
  var_primitiv[2,:] = np.array([ 120.0, -80.0,    40.0, -150.0])
  var_primitiv[3,:] = 0.0
  var_primitiv[4,:] = np.array([ 300.0,  340.0,  270.0,  410.0])
  var_primitiv[5,:] = var_primitiv[0,:]*GAS_CONSTANT*var_primitiv[4,:]

  var_primitiv_bd = np.zeros((6, num_face_bd))
  var_primitiv_bd[0,:] = np.array([1.1e-4, 0.9e-4, 1.6e-4])
  var_primitiv_bd[1,:] = np.array([ 580.0,  640.0,  470.0])
  var_primitiv_bd[2,:] = np.array([  90.0, -110.0,   60.0])
  var_primitiv_bd[3,:] = 0.0
  var_primitiv_bd[4,:] = np.array([ 310.0,  290.0,  380.0])
  var_primitiv_bd[5,:] = var_primitiv_bd[0,:]*GAS_CONSTANT*var_primitiv_bd[4,:]

  var_conserv = np.zeros((NUM_CONSERV, num_cell))
  for n_cell in range(0, num_cell):
    var_conserv[:,n_cell] = conservative_from_primitive(var_primitiv[:,n_cell])

  # 非定常項に非自明な値を通すため、過去の解は現在と少しずらす
  var_conserv_prev = np.zeros((2, NUM_CONSERV, num_cell))
  var_conserv_prev[0,:,:] = var_conserv*0.97
  var_conserv_prev[1,:,:] = var_conserv*0.94

  var_rhs = np.zeros((NUM_CONSERV, num_cell))
  var_rhs[:,:] = np.outer( np.array([1.0e-4, 2.0e-1, -1.5e-1, 0.0, 5.0e2]), \
                           np.array([1.0, -0.7, 1.3, -0.4]) )

  var_dt = np.array([1.0e-7, 1.4e-7, 0.8e-7, 1.9e-7])

  viscosity    = np.array([1.8e-5, 2.0e-5, 1.6e-5, 2.3e-5])
  viscosity_bd = np.array([1.9e-5, 1.7e-5, 2.1e-5])

  config = { 'time_integration': { 'kind_steady_mode': kind_steady_mode,                   \
                                   'kind_backward_difference': kind_backward_difference,   \
                                   'lusgs_beta': lusgs_beta,                               \
                                   'timestep_outer': TIMESTEP_OUTER } }

  return { 'config':         config,                                                  \
           'dimension_dict': {'num_conservative': NUM_CONSERV},                       \
           'geom_dict':      {'num_face_inner': num_face,                             \
                              'num_face_boundary': num_face_bd,                       \
                              'num_cell': num_cell,                                   \
                              'face2cell_inner': face2cell,                           \
                              'face2cell_boundary': face2cell_bd,                     \
                              'virtualcell_boundary': virtualcell_bd},                \
           'metrics_dict':   {'area_vec_inner': area_vec,                             \
                              'area_vec_boundary': area_vec_bd,                       \
                              'length_inner': length,                                 \
                              'length_boundary': length_bd,                           \
                              'volume_cell': volume},                                 \
           'gas_property_dict': {'specfic_heat_ratio': SPECIFIC_HEAT_RATIO,           \
                                 'specific_heat_volume': SPECIFIC_HEAT_VOLUM},        \
           'transport_coefficient_dict': {'viscosity': viscosity,                     \
                                          'viscosity_boundary': viscosity_bd},        \
           'var_primitiv':     var_primitiv,                                          \
           'var_primitiv_bd':  var_primitiv_bd,                                       \
           'var_conserv':      var_conserv,                                           \
           'var_conserv_prev': var_conserv_prev,                                      \
           'var_rhs':          var_rhs,                                               \
           'var_dt':           var_dt }


def time_term(case, n_cell):
  """対角の時間項・右辺の非定常項・面の寄与にかかる係数（テスト側で独立に組み立てる）"""

  config = case['config']['time_integration']
  volume = case['metrics_dict']['volume_cell'][n_cell]
  var_dt = case['var_dt'][n_cell]

  if config['kind_steady_mode'] == 'steady':
    return volume/var_dt, np.zeros(NUM_CONSERV), 0.50*config['lusgs_beta']

  if config['kind_backward_difference'] == '2nd_backward_diff':
    diag_time = 1.50*volume/TIMESTEP_OUTER + volume/var_dt
    dq_unst   = ( 1.50*case['var_conserv'][:,n_cell]           \
                - 2.00*case['var_conserv_prev'][0,:,n_cell]    \
                + 0.50*case['var_conserv_prev'][1,:,n_cell] )*volume/TIMESTEP_OUTER
  else:
    diag_time = volume/TIMESTEP_OUTER + volume/var_dt
    dq_unst   = ( case['var_conserv'][:,n_cell]                \
                - case['var_conserv_prev'][0,:,n_cell] )*volume/TIMESTEP_OUTER

  return diag_time, dq_unst, 0.50


def rhs_vector(case):
  """陰解系の右辺 b = -RHS - 非定常項 を 1 本のベクトルに並べる"""

  num_cell = case['geom_dict']['num_cell']
  rhs_vec  = np.zeros(NUM_CONSERV*num_cell)
  for n_cell in range(0, num_cell):
    _, dq_unst, _ = time_term(case, n_cell)
    rhs_vec[NUM_CONSERV*n_cell:NUM_CONSERV*(n_cell+1)] = -case['var_rhs'][:,n_cell] - dq_unst

  return rhs_vec


def assemble_operator(case, kind_dissipation):
  """
  D, L, U を密行列として独立に組み立てる。

  面 (a,b) の法線 n はセル b の外向き（セル a の内向き）。
    D_a += face_scale*( A(u_a,-n) + Diss )*S,  D_b += face_scale*( A(u_b,+n) + Diss )*S
    L[b,a] = 0.5*( A(u_a,+n) - Diss )*S        （b>a なので下三角）
    U[a,b] = 0.5*( A(u_b,-n) - Diss )*S        （a<b なので上三角）
  scalar 散逸では Diss=lambda*I とし、対角には A の項を含めない（従来の実装に合わせる）。
  """

  num_cell = case['geom_dict']['num_cell']
  num_face = case['geom_dict']['num_face_inner']
  face2cell = case['geom_dict']['face2cell_inner']
  area_vec  = case['metrics_dict']['area_vec_inner']
  length    = case['metrics_dict']['length_inner']
  prim      = case['var_primitiv']
  viscosity = case['transport_coefficient_dict']['viscosity']

  size    = NUM_CONSERV*num_cell
  diag    = np.zeros((size, size))
  lower   = np.zeros((size, size))
  upper   = np.zeros((size, size))
  identity = np.identity(NUM_CONSERV)

  def block(matrix, row_cell, col_cell):
    return matrix[NUM_CONSERV*row_cell:NUM_CONSERV*(row_cell+1),
                  NUM_CONSERV*col_cell:NUM_CONSERV*(col_cell+1)]

  # 内部面
  for n_face in range(0, num_face):
    area  = area_vec[0,n_face]
    vec   = area_vec[1:4,n_face]
    n_a   = face2cell[0,n_face]
    n_b   = face2cell[1,n_face]

    cons_a = conservative_from_primitive(prim[:,n_a])
    cons_b = conservative_from_primitive(prim[:,n_b])

    if kind_dissipation == 'scalar':
      leng = length[0,n_face] + length[1,n_face]
      eigen = max( max_eigenvalue(prim[:,n_a], viscosity[n_a], leng, vec), \
                   max_eigenvalue(prim[:,n_b], viscosity[n_b], leng, vec) )
      diss_a = eigen*identity
      diss_b = eigen*identity
      # scalar 版の対角は lambda*S のみを積む（A の項は面の総和で相殺する）
      block(diag, n_a, n_a)[:,:] += eigen*area*identity
      block(diag, n_b, n_b)[:,:] += eigen*area*identity
    else:
      diss = absolute_jacobian(prim[:,n_a], prim[:,n_b], vec)
      diss_a = diss
      diss_b = diss
      block(diag, n_a, n_a)[:,:] += ( jacobian_by_complex_step(cons_a, -vec) + diss )*area
      block(diag, n_b, n_b)[:,:] += ( jacobian_by_complex_step(cons_b,  vec) + diss )*area

    # 非対角: 前進は b の式に現れる dq_a、後退は a の式に現れる dq_b
    block(lower, n_b, n_a)[:,:] = 0.50*( jacobian_by_complex_step(cons_a,  vec) - diss_b )*area
    block(upper, n_a, n_b)[:,:] = 0.50*( jacobian_by_complex_step(cons_b, -vec) - diss_a )*area

  # 境界面（対角のみに効く）
  num_face_bd    = case['geom_dict']['num_face_boundary']
  face2cell_bd   = case['geom_dict']['face2cell_boundary']
  virtualcell_bd = case['geom_dict']['virtualcell_boundary']
  area_vec_bd    = case['metrics_dict']['area_vec_boundary']
  length_bd      = case['metrics_dict']['length_boundary']
  prim_bd        = case['var_primitiv_bd']
  viscosity_bd   = case['transport_coefficient_dict']['viscosity_boundary']

  for n_face in range(0, num_face_bd):
    area  = area_vec_bd[0,n_face]
    vec   = area_vec_bd[1:4,n_face]
    n_a   = face2cell_bd[0,n_face]
    vcell = float(virtualcell_bd[n_face])

    if kind_dissipation == 'scalar':
      prim_face = vcell*0.50*( prim[:,n_a] + prim_bd[:,n_face] ) + (1.0-vcell)*prim_bd[:,n_face]
      visc_face = vcell*0.50*( viscosity[n_a] + viscosity_bd[n_face] ) + (1.0-vcell)*viscosity_bd[n_face]
      leng      = length_bd[n_face]*( 1.0 + vcell )
      eigen     = max_eigenvalue(prim_face, visc_face, leng, vec)
      block(diag, n_a, n_a)[:,:] += eigen*area*identity
    else:
      diss   = absolute_jacobian(prim[:,n_a], prim_bd[:,n_face], vec)
      cons_a = conservative_from_primitive(prim[:,n_a])
      block(diag, n_a, n_a)[:,:] += ( jacobian_by_complex_step(cons_a, -vec) + diss )*area

  # 時間項と面の寄与の係数
  for n_cell in range(0, num_cell):
    diag_time, _, face_scale = time_term(case, n_cell)
    block(diag, n_cell, n_cell)[:,:] = diag_time*identity \
                                     + face_scale*block(diag, n_cell, n_cell)

  return diag, lower, upper


def run_solver(case, kind_dissipation):
  """get_diagonal + sweep_jacobian を実際に通す"""

  case['config']['time_integration']['kind_lusgs_dissipation'] = kind_dissipation

  num_cell = case['geom_dict']['num_cell']
  if kind_dissipation == 'matrix':
    var_diagonal = np.zeros((NUM_CONSERV, NUM_CONSERV, num_cell))
  else:
    var_diagonal = np.zeros(num_cell)
  var_dq = np.zeros((NUM_CONSERV, num_cell))

  var_diagonal, var_dq = lusgs_diagonal.get_diagonal(
      case['config'], case['dimension_dict'], case['geom_dict'], case['metrics_dict'],
      case['gas_property_dict'], case['transport_coefficient_dict'],
      case['var_primitiv'], case['var_primitiv_bd'],
      case['var_conserv'], case['var_conserv_prev'],
      case['var_rhs'], case['var_dt'], var_diagonal, var_dq)

  dq_initial = var_dq.copy()

  var_dq = lusgs_sweep.sweep_jacobian(
      case['config'], case['dimension_dict'], case['geom_dict'], case['metrics_dict'],
      case['gas_property_dict'], case['transport_coefficient_dict'],
      case['var_primitiv'], case['var_conserv'], var_diagonal, var_dq)

  return var_diagonal, dq_initial, var_dq


def flatten(var_dq):
  # (num_conserv, num_cell) -> セルごとに 5 成分ずつ並べた 1 本のベクトル

  return var_dq.T.reshape(-1)


MODES = [
  pytest.param('steady',   1.00, id='steady-beta1'),
  pytest.param('steady',   1.01, id='steady-beta1.01'),
  pytest.param('unsteady', 1.01, id='unsteady-bdf2'),
]


# ---------------------------------------------------------------- 検査

@pytest.mark.parametrize('kind_dissipation', ['scalar', 'matrix'])
@pytest.mark.parametrize('kind_steady_mode, lusgs_beta', MODES)
def test_diagonal_matches_independent_assembly(kind_dissipation, kind_steady_mode, lusgs_beta):
  # 対角 D が独立に組んだものと一致すること（matrix 版は逆行列で保持される）

  case = make_case(kind_steady_mode=kind_steady_mode, lusgs_beta=lusgs_beta)
  var_diagonal, _, _ = run_solver(case, kind_dissipation)
  diag, _, _ = assemble_operator(case, kind_dissipation)

  num_cell = case['geom_dict']['num_cell']
  for n_cell in range(0, num_cell):
    expected = diag[NUM_CONSERV*n_cell:NUM_CONSERV*(n_cell+1),
                    NUM_CONSERV*n_cell:NUM_CONSERV*(n_cell+1)]
    if kind_dissipation == 'matrix':
      expected_inv = np.linalg.inv(expected)
      # 成分の大きさが数桁にわたるので、許容差はブロック全体の大きさに対して取る
      np.testing.assert_allclose(var_diagonal[:,:,n_cell], expected_inv, rtol=1.0e-9, \
                                 atol=1.0e-12*np.abs(expected_inv).max())
    else:
      np.testing.assert_allclose(var_diagonal[n_cell], expected[0,0], rtol=1.0e-12)


@pytest.mark.parametrize('kind_dissipation', ['scalar', 'matrix'])
@pytest.mark.parametrize('kind_steady_mode, lusgs_beta', MODES)
def test_dq_is_initialized_with_the_inverse_diagonal(kind_dissipation, kind_steady_mode, lusgs_beta):
  # スイープ前の dq が D^-1*b であること

  case = make_case(kind_steady_mode=kind_steady_mode, lusgs_beta=lusgs_beta)
  _, dq_initial, _ = run_solver(case, kind_dissipation)
  diag, _, _ = assemble_operator(case, kind_dissipation)

  rhs_vec = rhs_vector(case)
  np.testing.assert_allclose(diag @ flatten(dq_initial), rhs_vec, rtol=1.0e-9, \
                             atol=1.0e-12*np.abs(rhs_vec).max())


@pytest.mark.parametrize('kind_dissipation', ['scalar', 'matrix'])
@pytest.mark.parametrize('kind_steady_mode, lusgs_beta', MODES)
def test_sweep_solves_the_factorized_system(kind_dissipation, kind_steady_mode, lusgs_beta):
  # 本題: スイープの結果が (D+L)*D^-1*(D+U)*dq = b を厳密に満たすこと。
  # セル番号・法線の向き・掃引順のいずれかを取り違えると成立しない

  case = make_case(kind_steady_mode=kind_steady_mode, lusgs_beta=lusgs_beta)
  _, _, var_dq = run_solver(case, kind_dissipation)
  diag, lower, upper = assemble_operator(case, kind_dissipation)

  operator = (diag + lower) @ np.linalg.inv(diag) @ (diag + upper)

  rhs_vec = rhs_vector(case)
  np.testing.assert_allclose(operator @ flatten(var_dq), rhs_vec, rtol=1.0e-8, \
                             atol=1.0e-11*np.abs(rhs_vec).max())


@pytest.mark.parametrize('kind_dissipation', ['scalar', 'matrix'])
def test_sweep_is_not_trivially_the_diagonal_solve(kind_dissipation):
  # 上の一致が「非対角がゼロで自明に成り立っている」のではないことの確認

  case = make_case()
  _, dq_initial, var_dq = run_solver(case, kind_dissipation)

  change = np.abs(flatten(var_dq) - flatten(dq_initial)).max()
  scale  = np.abs(flatten(dq_initial)).max()

  assert change > 1.0e-3*scale, '非対角の寄与が効いていない'


# ---------------------------------------------------------------- 実メッシュ上の不変条件

DATA_DIR    = pathlib.Path(__file__).parent/'data'
FILE_CONFIG = DATA_DIR/'config_rhs_snapshot.yml'
FILE_MESH   = DATA_DIR/'chimera.msh'


@pytest.fixture(scope='module')
def frozen_mesh():
  """凍結した 25 セルのメッシュを読み込む（tests/data の設定と格子）"""

  from slow.meshdata.meshdata import meshdata
  from slow.orbital.orbital import orbital

  orbital_obj  = orbital()
  meshdata_obj = meshdata()

  config = orbital_obj.read_config_yaml(str(FILE_CONFIG))
  config['meshdata_io']['filename_mesh'] = str(FILE_MESH)

  _, _, geom_dict, metrics_dict = meshdata_obj.set_mesh_routine(config)

  return geom_dict, metrics_dict


def test_face_ordering_makes_the_sweep_an_exact_triangular_solve(frozen_mesh):
  """
  スイープが黙って依存している面の並びを固定する。
    face2cell_inner[0,n] < face2cell_inner[1,n] かつ第 0 行が単調非減少
  これが成り立つ限り、前進スイープは厳密な下三角 Gauss-Seidel 解法になる。
  格子の番号付けを変える変更（領域分割など）ではここが最初に壊れる。
  """

  geom_dict, _ = frozen_mesh
  face2cell = geom_dict['face2cell_inner']

  assert np.all( face2cell[0,:] < face2cell[1,:] ), \
         '内部面の第 0 セル番号が第 1 セル番号より小さくない'
  assert np.all( np.diff(face2cell[0,:]) >= 0 ), \
         '内部面の第 0 セル番号が単調非減少でない'


def test_outward_normals_close_each_cell(frozen_mesh):
  """
  各セルについて外向き法線ベクトルの面積重み和がゼロになること。

  面の法線 area_vec[1:4] はセル a に対して内向き（rhs も lusgs も
  セル a の外向きを -n として扱う）。境界面も同じ向きであることを含めて確認する。
  この向きを取り違えると陰解演算子の A^+ と A^- が入れ替わる。
  """

  geom_dict, metrics_dict = frozen_mesh

  num_cell     = geom_dict['num_cell']
  num_face     = geom_dict['num_face_inner']
  num_face_bd  = geom_dict['num_face_boundary']
  face2cell    = geom_dict['face2cell_inner']
  face2cell_bd = geom_dict['face2cell_boundary']
  area_vec     = metrics_dict['area_vec_inner']
  area_vec_bd  = metrics_dict['area_vec_boundary']

  closure = np.zeros((3, num_cell))
  for n_face in range(0, num_face):
    n_a = face2cell[0,n_face]
    n_b = face2cell[1,n_face]
    closure[:,n_a] -= area_vec[1:4,n_face]*area_vec[0,n_face]
    closure[:,n_b] += area_vec[1:4,n_face]*area_vec[0,n_face]

  for n_face in range(0, num_face_bd):
    n_a = face2cell_bd[0,n_face]
    closure[:,n_a] -= area_vec_bd[1:4,n_face]*area_vec_bd[0,n_face]

  np.testing.assert_allclose(closure, 0.0, atol=1.0e-12*area_vec[0,:].max())
