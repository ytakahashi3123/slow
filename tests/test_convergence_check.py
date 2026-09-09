#!/usr/bin/env python3

# Program to verify the convergence criteria of the inner and outer loops
#
# 相対／絶対の判定は内側ループと外側ループで同じ向きでなければならない。
# 外側ループは分岐が逆になっており、flag_convergence_relative_outerloop: True で
# 絶対値を、False で相対値を見ていた。残差はエネルギーで 1e9 の大きさなので、
# 既定の設定 (relative: True, 1e-8) では実質いつまでも収束しない状態だった。

# Author: Y.Takahashi, Hokkaido University
# Date; 2026/09/09

import numpy as np
import pytest

from slow.orbital.orbital import orbital


NUM_CONSERV = 5

# エネルギーの残差の大きさ（デバッグ格子の実測に合わせる）
RESIDUAL_INIT = 5.93e9


def make_config(flag_relative_inner=True, criterion_inner=1.0e-8,
                flag_relative_outer=True, criterion_outer=1.0e-8):
  return { 'time_integration': { 'flag_convergence_relative_innerloop': flag_relative_inner, \
                                 'criterion_convergence_innerloop': criterion_inner,         \
                                 'flag_convergence_relative_outerloop': flag_relative_outer, \
                                 'criterion_convergence_outerloop': criterion_outer } }


def make_orbital(sum_rhs_init=RESIDUAL_INIT, sum_dq_init=1.0e3):
  # check_convergence_* は初期残差を自分のインスタンス属性から読む
  orbital_obj = orbital()
  orbital_obj.sum_rhs_init = np.full(NUM_CONSERV, sum_rhs_init)
  orbital_obj.sum_dq_init  = np.full(NUM_CONSERV, sum_dq_init)

  return orbital_obj


def residual(value):
  # 判定に使われるのはエネルギー成分 [4] のみ
  return np.array([0.0, 0.0, 0.0, 0.0, value])


# ---------------------------------------------------------------- 外側ループ

def test_outer_relative_criterion_compares_against_the_initial_residual():
  # relative: True なら初期残差との比で判定する

  orbital_obj = make_orbital()
  config      = make_config(flag_relative_outer=True, criterion_outer=1.0e-8)

  # 比が基準を下回る
  assert orbital_obj.check_convergence_outer(config, False, residual(RESIDUAL_INIT*1.0e-9))
  # 比が基準を上回る
  assert not orbital_obj.check_convergence_outer(config, False, residual(RESIDUAL_INIT*1.0e-7))


def test_outer_relative_criterion_is_not_the_absolute_residual():
  # 相対と絶対で判定が分かれる値を与え、相対指定のときに絶対値を見ていないことを確かめる

  orbital_obj = make_orbital(sum_rhs_init=1.0e-4)
  config      = make_config(flag_relative_outer=True, criterion_outer=1.0e-8)

  # 絶対値 1e-9 は基準 1e-8 を下回るが、初期残差 1e-4 との比は 1e-5 で基準に届かない
  assert not orbital_obj.check_convergence_outer(config, False, residual(1.0e-9))


def test_outer_absolute_criterion_compares_the_residual_itself():
  # relative: False なら残差そのもので判定する

  orbital_obj = make_orbital()
  config      = make_config(flag_relative_outer=False, criterion_outer=1.0e-8)

  assert orbital_obj.check_convergence_outer(config, False, residual(1.0e-9))
  # 初期残差との比では収束しているが、絶対値では基準に届かない
  assert not orbital_obj.check_convergence_outer(config, False, residual(RESIDUAL_INIT*1.0e-12))


def test_outer_check_never_clears_a_flag_that_is_already_set():
  # 既に立っているフラグを倒さない（呼び出し側は False を渡してくる）

  orbital_obj = make_orbital()
  config      = make_config(flag_relative_outer=True, criterion_outer=1.0e-8)

  assert orbital_obj.check_convergence_outer(config, True, residual(RESIDUAL_INIT))


# ---------------------------------------------------------------- 内側ループ

def test_inner_relative_criterion_compares_against_the_initial_deltaq():
  orbital_obj = make_orbital(sum_dq_init=1.0e3)
  config      = make_config(flag_relative_inner=True, criterion_inner=1.0e-8)

  assert orbital_obj.check_convergence_inner(config, False, residual(1.0e3*1.0e-9))
  assert not orbital_obj.check_convergence_inner(config, False, residual(1.0e3*1.0e-7))


def test_inner_absolute_criterion_compares_the_deltaq_itself():
  orbital_obj = make_orbital(sum_dq_init=1.0e3)
  config      = make_config(flag_relative_inner=False, criterion_inner=1.0e-8)

  assert orbital_obj.check_convergence_inner(config, False, residual(1.0e-9))
  assert not orbital_obj.check_convergence_inner(config, False, residual(1.0e-7))


# ---------------------------------------------------------------- 両者の向きがそろっていること

@pytest.mark.parametrize('flag_relative', [True, False])
def test_inner_and_outer_use_the_same_sense_of_relative(flag_relative):
  """
  同じ初期値・同じ基準・同じ現在値を与えたとき、内側と外側で判定が一致すること。
  片方だけ相対と絶対が入れ替わっていると、この検査で落ちる。
  """

  init      = 1.0e6
  criterion = 1.0e-8
  orbital_obj = make_orbital(sum_rhs_init=init, sum_dq_init=init)
  config = make_config(flag_relative_inner=flag_relative, criterion_inner=criterion,
                       flag_relative_outer=flag_relative, criterion_outer=criterion)

  for value in (init, init*1.0e-4, init*1.0e-9, 1.0e-9, 1.0e-7):
    assert orbital_obj.check_convergence_outer(config, False, residual(value)) \
        == orbital_obj.check_convergence_inner(config, False, residual(value)), \
           f'残差 {value} で内側と外側の判定が食い違う'
