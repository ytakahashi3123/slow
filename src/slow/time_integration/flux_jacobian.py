#!/usr/bin/env python3

# Program to calculate the convective flux Jacobian on a cell interface

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/31


from slow.general.jit import kernel


@kernel
def set_flux_jacobian(jacobian, specfic_heat_ratio, eigenvalue, cvel, uvel, vvel, wvel, enth, vecx, vecy, vecz):
  """
  Set the split convective flux Jacobian A(n) - lambda*I used in the LU-SGS sweep.

  A(n)=d(F.n)/dQ は保存変数 Q=(rho, rho*u, rho*v, rho*w, E) に対する法線流束のヤコビ行列。
  面ループ内から呼ばれるため、確保済みの jacobian(5,5) をその場で上書きする。

  --jacobian: (num_conservative, num_conservative) の配列。戻り値ではなくこれを書き換える
  --eigenvalue: 対角から差し引く最大固有値 lambda
  --cvel: 反変速度 u.n
  --enth: 全エンタルピー H=(E+p)/rho（orbital.get_enthalpy の戻り値）
  --vecx, vecy, vecz: 面の単位法線ベクトル

  A(n) は法線について奇関数（A(-n)=-A(n)）である。後退スイープはこの性質を使い、
  法線を反転させることで自セルの外向き法線に対する A^- を得ている。
  """

  # Set derivatives
  # -- Energy derivatives
  pett =  specfic_heat_ratio - 1.0
  # -- Dynamic pressure
  q    =  0.50*( uvel*uvel + vvel*vvel + wvel*wvel )
  # -- Total density derivatives: dp/drho|_(rho*u,E) = (gamma-1)*q
  prho =  pett * q


  # Jacobian matrix (Left side)
  jacobian[0,0] = 0.0 - eigenvalue
  jacobian[0,1] = vecx
  jacobian[0,2] = vecy
  jacobian[0,3] = vecz
  jacobian[0,4] = 0.0

  jacobian[1,0] =  vecx*prho            - cvel*uvel
  jacobian[1,1] = -vecx*(pett-1.0)*uvel + cvel - eigenvalue
  jacobian[1,2] = -vecx*pett*vvel       + vecy*uvel
  jacobian[1,3] = -vecx*pett*wvel       + vecz*uvel
  jacobian[1,4] =  vecx*pett

  jacobian[2,0] =  vecy*prho            - cvel*vvel
  jacobian[2,1] = -vecy*pett*uvel       + vecx*vvel
  jacobian[2,2] = -vecy*(pett-1.0)*vvel + cvel - eigenvalue
  jacobian[2,3] = -vecy*pett*wvel       + vecz*vvel
  jacobian[2,4] =  vecy*pett

  jacobian[3,0] =  vecz*prho            - cvel*wvel
  jacobian[3,1] = -vecz*pett*uvel       + vecx*wvel
  jacobian[3,2] = -vecz*pett*vvel       + vecy*wvel
  jacobian[3,3] = -vecz*(pett-1.0)*wvel + cvel - eigenvalue
  jacobian[3,4] =  vecz*pett

  jacobian[4,0] =  cvel*(prho-enth)
  jacobian[4,1] = -cvel*pett*uvel + vecx*enth
  jacobian[4,2] = -cvel*pett*vvel + vecy*enth
  jacobian[4,3] = -cvel*pett*wvel + vecz*enth
  jacobian[4,4] =  cvel*(pett+1.0) - eigenvalue

  return
