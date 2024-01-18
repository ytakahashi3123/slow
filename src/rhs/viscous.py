#!/usr/bin/env python3

# Program to calculate viscous flux on interface

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/25

import numpy as np
from orbital.orbital import orbital

@orbital.time_measurement_decorated
def flux_viscous(config, dimension_dict, geom_dict, metrics_dict, gas_property_dict, transport_coefficient_dict, var_primitiv, var_primitiv_bd, var_gradient, var_rhs):

  def get_stress_tensor():
    # Stress tensors: \mu du_i/dx_j
    tau_tmp  = 1.0/3.0*(grad_face[0,1] + grad_face[1,2] + grad_face[2,3])
    tau[0,0] = 2.0*viscosty_face*( grad_face[0,1] - tau_tmp )
    tau[0,1] =     viscosty_face*( grad_face[0,2] + grad_face[1,1] )
    tau[0,2] =     viscosty_face*( grad_face[0,3] + grad_face[2,1] )
    tau[1,0] = tau[0,1]
    tau[1,1] = 2.0*viscosty_face*( grad_face[1,2] - tau_tmp )
    tau[1,2] =     viscosty_face*( grad_face[1,3] + grad_face[2,2] )
    tau[2,0] = tau[0,2]
    tau[2,1] = tau[1,2]
    tau[2,2] = 2.0*viscosty_face*( grad_face[2,3] - tau_tmp )

    # Viscous stress Work
    uvel_face      = 0.50*(var_a[1]+var_b[1])
    vvel_face      = 0.50*(var_a[2]+var_b[2])
    wvel_face      = 0.50*(var_a[3]+var_b[3])
    work_stress[0] = tau[0,0]*uvel_face + tau[0,1]*vvel_face + tau[0,2]*wvel_face
    work_stress[1] = tau[1,0]*uvel_face + tau[1,1]*vvel_face + tau[1,2]*wvel_face
    work_stress[2] = tau[2,0]*uvel_face + tau[2,1]*vvel_face + tau[2,2]*wvel_face

    return

  def get_heat_flux():
    # Heat flux: \lambda dT/dx
    heat_flux[0] = thermal_cond_face*grad_face[0,4]
    heat_flux[1] = thermal_cond_face*grad_face[1,4]
    heat_flux[2] = thermal_cond_face*grad_face[2,4]

    return


  # Main routine
  
  # Input parameters
  num_conserv  = dimension_dict['num_conservative']
  num_primitiv = dimension_dict['num_primitive']

  num_face     = geom_dict['num_face_inner']
  num_face_bd  = geom_dict['num_face_boundary']
  num_cell     = geom_dict['num_cell']
  face2cell    = geom_dict['face2cell_inner'] 
  face2cell_bd = geom_dict['face2cell_boundary']

  area_vec    = metrics_dict['area_vec_inner']
  area_vec_bd = metrics_dict['area_vec_boundary']
  length      = metrics_dict['length_inner']
  length_bd   = metrics_dict['length_boundary']
  volume      = metrics_dict['volume_cell']

  specfic_heat_ratio  = gas_property_dict['specfic_heat_ratio']
  specific_heat_volum = gas_property_dict['specific_heat_volume']

  viscosity           = transport_coefficient_dict['viscosity']
  thermal_cond        = transport_coefficient_dict['thermal_conductivity']
  viscosity_bd        = transport_coefficient_dict['viscosity_boundary']
  thermal_cond_bd     = transport_coefficient_dict['thermal_conductivity_boundary']

  # Initialize
  grad_face   = np.zeros(3*3).reshape(3,3)
  tau         = np.zeros(3*3).reshape(3,3)
  work_stress = np.zeros(3).reshape(3)
  heat_flux   = np.zeros(3).reshape(3)
  flux_rhs    = np.linspace(0, 0, num_conserv, dtype=float)
  

  # Flux calculation
  # --Inner loop
  for n_face in range(0,num_face):

    # Face area vector
    area   = area_vec[0,n_face]
    vecx   = area_vec[1,n_face]
    vecy   = area_vec[2,n_face]
    vecz   = area_vec[3,n_face]

    # Spatial gradient variables
    # Calculate values on face, interpolated from left (self) and right sides (neigboring) cell
    # --from self cell side
    n_cell_a = face2cell[0,n_face]
    var_a    = var_primitiv[:,n_cell_a]
    # --from neigboring cell
    n_cell_b = face2cell[1,n_face]
    var_b    = var_primitiv[:,n_cell_b]


    # Gradient variables on face
    # 単純な平均で評価している
    # grad_face[0,0]: drho/dx,  grad_face[1,0]: drho/dy,  grad_face[2,0]: drho/dz
    # grad_face[0,1]: du/dx,    grad_face[1,1]: du/dy,    grad_face[2,1]: du/dz
    # grad_face[0,2]: dv/dx,    grad_face[1,2]: dv/dy,    grad_face[2,2]: dv/dz
    # grad_face[0,3]: dw/dx,    grad_face[1,3]: dw/dy,    grad_face[2,3]: dw/dz
    # grad_face[0,4]: dT/dx,    grad_face[1,4]: dT/dy,    grad_face[2,4]: dT/dz
    # grad_face[0,5]: dp/dx,    grad_face[1,5]: dp/dy,    grad_face[2,5]: dp/dz
    grad_face = 0.50*(var_gradient[:,:,n_cell_a] + var_gradient[:,:,n_cell_b])


    # Transport coefficients on face
    # 熱伝導率は相加平均で評価したが、相乗平均でも良いかも
    viscosty_face     = 0.50*( viscosity[n_cell_a] + viscosity[n_cell_b] )
    thermal_cond_face = 0.50*( thermal_cond[n_cell_a] + thermal_cond[n_cell_b] )


    # Stress tensors-->tau, work_stress
    get_stress_tensor()


    # Heat flux-->heat_flux
    get_heat_flux()


    # Calculate fluxes
    flux_rhs[0] = 0.0
    flux_rhs[1] = - ( tau[0,0]*vecx + tau[0,1]*vecy + tau[0,2]*vecz )
    flux_rhs[2] = - ( tau[1,0]*vecx + tau[1,1]*vecy + tau[1,2]*vecz )
    flux_rhs[3] = - ( tau[2,0]*vecx + tau[2,1]*vecy + tau[2,2]*vecz )
    flux_rhs[4] = - ( (work_stress[0]+heat_flux[0])*vecx + (work_stress[1]+heat_flux[1])*vecy + (work_stress[2]+heat_flux[2])*vecz )

    var_rhs[:,n_cell_a] = var_rhs[:,n_cell_a] - flux_rhs[:]*area
    var_rhs[:,n_cell_b] = var_rhs[:,n_cell_b] + flux_rhs[:]*area


  # --Boundary loop
  for n_face in range(0,num_face_bd):

    # Face area vector
    area = area_vec_bd[0,n_face]
    vecx = area_vec_bd[1,n_face]
    vecy = area_vec_bd[2,n_face]
    vecz = area_vec_bd[3,n_face]

    # Calculate values on face, interpolated from left (self) and right sides (neigboring) cell
    # --from self cell side
    n_cell_a = face2cell_bd[0,n_face]
    var_a    = var_primitiv[:,n_cell_a]
    # --from neigboring cell
    var_b    = var_primitiv_bd[:,n_face]


    # Gradient variables on face
    # 隣接セルの勾配をそのまま使う
    # grad_face[0,0]: drho/dx,  grad_face[1,0]: drho/dy,  grad_face[2,0]: drho/dz
    # grad_face[0,1]: du/dx,    grad_face[1,1]: du/dy,    grad_face[2,1]: du/dz
    # grad_face[0,2]: dv/dx,    grad_face[1,2]: dv/dy,    grad_face[2,2]: dv/dz
    # grad_face[0,3]: dw/dx,    grad_face[1,3]: dw/dy,    grad_face[2,3]: dw/dz
    # grad_face[0,4]: dT/dx,    grad_face[1,4]: dT/dy,    grad_face[2,4]: dT/dz
    # grad_face[0,5]: dp/dx,    grad_face[1,5]: dp/dy,    grad_face[2,5]: dp/dz
    grad_face = var_gradient[:,:,n_cell_a] 


    # Transport coefficients on face
    # 熱伝導率は相加平均で評価したが、相乗平均でも良いかも
    viscosty_face     = 0.50*( viscosity[n_cell_a] + viscosity_bd[n_face] )
    thermal_cond_face = 0.50*( thermal_cond[n_cell_a] + thermal_cond_bd[n_face] )


    # Stress tensors-->tau, work_stress
    get_stress_tensor()


    # Heat flux-->heat_flux
    get_heat_flux()


    # Calculate fluxes
    flux_rhs[0] = 0.0
    flux_rhs[1] = - ( tau[0,0]*vecx + tau[0,1]*vecy + tau[0,2]*vecz )
    flux_rhs[2] = - ( tau[1,0]*vecx + tau[1,1]*vecy + tau[1,2]*vecz )
    flux_rhs[3] = - ( tau[2,0]*vecx + tau[2,1]*vecy + tau[2,2]*vecz )
    flux_rhs[4] = - ( (work_stress[0]+heat_flux[0])*vecx + (work_stress[1]+heat_flux[1])*vecy + (work_stress[2]+heat_flux[2])*vecz )

    var_rhs[:,n_cell_a] = var_rhs[:,n_cell_a] - flux_rhs[:]*area


  return var_rhs

