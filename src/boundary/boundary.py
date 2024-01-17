#!/usr/bin/env python3

# Module to set boundary condition

# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/11

import numpy as np
from orbital.orbital import orbital

class boundary(orbital):


  def __init__(self):
    print("Calling class: boundary")

    # not used
    self.ID_BOUNDARY_FREESTREAM   = 1
    self.ID_BOUNDARY_AXISYMMETRIC = 6
    self.ID_BOUNDARY_OUTLET       = 10
    self.ID_BOUNDARY_WALL_FIX     = 20

    # Boundary name-->orbital/orbital.py
    
    return


  def boundary_condition(self, config, geom_dict, metrics_dict, gas_property_dict, var_primitiv, var_conserv, var_primitiv_bd):
    # Boundary condition

    # 以下はLoop内でやる意味はないのでinitなどでやったほうがよいかも
    # Freestream
    #key_freestream = config['boundary_conditions'].get( self.KIND_BOUNDARY_FREESTREAM )
    #if key_freestream != None :
    #  density_freestream       = config['boundary_conditions'][self.KIND_BOUNDARY_FREESTREAM]['density']
    #  velocity_freestream      = config['boundary_conditions'][self.KIND_BOUNDARY_FREESTREAM]['velocity']
    #  temperature_freestream   = config['boundary_conditions'][self.KIND_BOUNDARY_FREESTREAM]['temperature']

    # Inlet given by total pressure and temperature
    #key_total_press_temp = config['boundary_conditions'].get( self.KIND_BOUNDARY_TOTAL_PRESS_TEMP )
    #if key_total_press_temp != None :
    #  pressure_total_press_temp    = config['boundary_conditions'][self.KIND_BOUNDARY_TOTAL_PRESS_TEMP]['pressure']
    #  temperature_total_press_temp = config['boundary_conditions'][self.KIND_BOUNDARY_TOTAL_PRESS_TEMP]['temperature']

    # Wall condition
    #key_wall = config['boundary_conditions'].get( self.KIND_BOUNDARY_WALL_FIX )
    #if key_wall != None :
    #  temperature_wall = config['boundary_conditions'][self.KIND_BOUNDARY_WALL_FIX]['temperature']

    # Slip wall
    #key_wall_slip = config['boundary_conditions'].get( self.KIND_BOUNDARY_WALL_SLIP )

    # Symmetry condition
    #key_symmetry = config['boundary_conditions'].get( self.KIND_BOUNDARY_SYMMETRY )

    # Outlet condition
    #key_outlet = config['boundary_conditions'].get( self.KIND_BOUNDARY_OUTLET )

    # Outlet by fixed pressure 
    #key_outlet_press_fixed = config['boundary_conditions'].get( self.KIND_BOUNDARY_OUTLET_PRESS_FIXED )
    #if key_outlet_press_fixed != None :
    #  pressure_outlet = config['boundary_conditions'][self.KIND_BOUNDARY_OUTLET_PRESS_FIXED]['pressure']

    # Ambient condition
    #key_ambient = config['boundary_conditions'].get( self.KIND_BOUNDARY_AMBIENT )
    #if key_ambient != None :
    #  pressure_ambient    = config['boundary_conditions'][self.KIND_BOUNDARY_AMBIENT]['pressure']
    #  temperature_ambient = config['boundary_conditions'][self.KIND_BOUNDARY_AMBIENT]['temperature']


    num_face_bd      = geom_dict['num_face_boundary']
    face2cell_bd     = geom_dict['face2cell_boundary']
    tag2name_physics = geom_dict['tag2name_physics']

    area_vec_bd = metrics_dict['area_vec_boundary']

    specfic_heat_ratio  = gas_property_dict['specfic_heat_ratio']
    gas_constant        = gas_property_dict['gas_constant']
    specific_heat_volum = gas_property_dict['specific_heat_volume']


    #--Primitive variables
    for n_face in range(0,num_face_bd):
      n_cell   = face2cell_bd[0,n_face]
      bd_name  = tag2name_physics[ face2cell_bd[1,n_face] ]
      bd_kind  = config['boundary_conditions'][bd_name]['kind']

      vecx   = area_vec_bd[1,n_face]
      vecy   = area_vec_bd[2,n_face]
      vecz   = area_vec_bd[3,n_face]
      vect1x = area_vec_bd[4,n_face]
      vect1y = area_vec_bd[5,n_face]
      vect1z = area_vec_bd[6,n_face]

      if( bd_kind == self.KIND_BOUNDARY_FREESTREAM ):
        # Freestream condition gives density, velocity, and temperature in control file
        density_freestream       = config['boundary_conditions'][bd_name]['density']
        velocity_freestream      = config['boundary_conditions'][bd_name]['velocity']
        temperature_freestream   = config['boundary_conditions'][bd_name]['temperature']

        var_primitiv_bd[0,n_face] = density_freestream
        var_primitiv_bd[1,n_face] = velocity_freestream[0]
        var_primitiv_bd[2,n_face] = velocity_freestream[1]
        var_primitiv_bd[3,n_face] = velocity_freestream[2]
        var_primitiv_bd[4,n_face] = temperature_freestream
        var_primitiv_bd[5,n_face] = self.get_pressure_eos(density_freestream, gas_constant, temperature_freestream)


      elif( bd_kind == self.KIND_BOUNDARY_TOTAL_PRESS_TEMP ):
        # Inlet parameters given by total pressure and temperature
        pressure_total_press_temp    = config['boundary_conditions'][bd_name]['pressure']
        temperature_total_press_temp = config['boundary_conditions'][bd_name]['temperature']

        gamma        =  specfic_heat_ratio-1.0
        gamma_ratio  = (specfic_heat_ratio-1.0)/specfic_heat_ratio
        pressure     = min( var_primitiv[5,n_cell], pressure_total_press_temp )
        mach_number2 = 2.0/gamma*( (pressure_total_press_temp/pressure)**gamma_ratio-1.0 )
        temperature  = temperature_total_press_temp/(1.0+0.50*gamma*mach_number2)
        density      = self.get_density_eos(pressure, gas_constant, temperature)
        velocity     = np.sqrt(mach_number2*specfic_heat_ratio*pressure/density)

        var_primitiv_bd[0,n_face] = density
        var_primitiv_bd[1,n_face] = velocity*vecx
        var_primitiv_bd[2,n_face] = velocity*vecy
        var_primitiv_bd[3,n_face] = velocity*vecz
        var_primitiv_bd[4,n_face] = temperature
        var_primitiv_bd[5,n_face] = pressure


      elif( bd_kind == self.KIND_BOUNDARY_SYMMETRY or bd_kind == self.KIND_BOUNDARY_WALL_SLIP ):
        # Symmetry or slip wall condition
        uvel_tmp =-( var_primitiv[1,n_cell]*vecx   + var_primitiv[2,n_cell]*vecy   + var_primitiv[3,n_cell]*vecz   )
        vvel_tmp = ( var_primitiv[1,n_cell]*vect1x + var_primitiv[2,n_cell]*vect1y + var_primitiv[3,n_cell]*vect1z )
        wvel_tmp =  0.0
        var_primitiv_bd[0,n_face] = var_primitiv[0,n_cell]
        var_primitiv_bd[1,n_face] = uvel_tmp*vecx + vvel_tmp*vect1x
        var_primitiv_bd[2,n_face] = uvel_tmp*vecy + vvel_tmp*vect1y
        var_primitiv_bd[3,n_face] = uvel_tmp*vecz + vvel_tmp*vect1z
        var_primitiv_bd[4,n_face] = var_primitiv[4,n_cell]
        var_primitiv_bd[5,n_face] = var_primitiv[5,n_cell]


      elif( bd_kind == self.KIND_BOUNDARY_OUTLET ):
        # Gradient-free (Zeroth extrapolation)
        uvel_tmp =  min( 0.0, ( var_primitiv[1,n_cell]*vecx   + var_primitiv[2,n_cell]*vecy   + var_primitiv[3,n_cell]*vecz   ) )
        vvel_tmp =            ( var_primitiv[1,n_cell]*vect1x + var_primitiv[2,n_cell]*vect1y + var_primitiv[3,n_cell]*vect1z )
        wvel_tmp =             0.0

        var_primitiv_bd[0,n_face] = var_primitiv[0,n_cell]
        var_primitiv_bd[1,n_face] = uvel_tmp*vecx + vvel_tmp*vect1x
        var_primitiv_bd[2,n_face] = uvel_tmp*vecy + vvel_tmp*vect1y
        var_primitiv_bd[3,n_face] = uvel_tmp*vecz + vvel_tmp*vect1z
        var_primitiv_bd[4,n_face] = var_primitiv[4,n_cell]
        var_primitiv_bd[5,n_face] = var_primitiv[5,n_cell]


      elif( bd_kind == self.KIND_BOUNDARY_OUTLET_PRESS_FIXED ):
        # Outlet by static pressure fixed
        pressure_outlet = config['boundary_conditions'][bd_name]['pressure']

        density     = var_primitiv[0,n_cell]
        pressure    = pressure_outlet
        temperature = pressure/( gas_constant*density )

        uvel_tmp =  min( 0.0, ( var_primitiv[1,n_cell]*vecx   + var_primitiv[2,n_cell]*vecy   + var_primitiv[3,n_cell]*vecz   ) )
        vvel_tmp =            ( var_primitiv[1,n_cell]*vect1x + var_primitiv[2,n_cell]*vect1y + var_primitiv[3,n_cell]*vect1z )
        wvel_tmp =             0.0

        var_primitiv_bd[0,n_face] = density
        var_primitiv_bd[1,n_face] = uvel_tmp*vecx + vvel_tmp*vect1x
        var_primitiv_bd[2,n_face] = uvel_tmp*vecy + vvel_tmp*vect1y
        var_primitiv_bd[3,n_face] = uvel_tmp*vecz + vvel_tmp*vect1z
        var_primitiv_bd[4,n_face] = temperature
        var_primitiv_bd[5,n_face] = pressure


      elif( bd_kind == self.KIND_BOUNDARY_AMBIENT ):
        # Ambient condition
        pressure_ambient    = config['boundary_conditions'][bd_name]['pressure']
        temperature_ambient = config['boundary_conditions'][bd_name]['temperature']

        pressure    = pressure_ambient
        temperature = temperature_ambient
        density     = self.get_density_eos(pressure, gas_constant, temperature)
        velocity    = 0.0
        
        var_primitiv_bd[0,n_face] = density
        var_primitiv_bd[1,n_face] = velocity*vecx
        var_primitiv_bd[2,n_face] = velocity*vecy
        var_primitiv_bd[3,n_face] = velocity*vecz
        var_primitiv_bd[4,n_face] = temperature
        var_primitiv_bd[5,n_face] = pressure 


      elif( bd_kind == self.KIND_BOUNDARY_WALL_FIX ):
        # Wall condition
        temperature_wall = config['boundary_conditions'][bd_name]['temperature']

        pressure_tmp     = var_primitiv[5,n_cell]
        gas_constant_tmp = gas_constant
        temperature_tmp  = temperature_wall
        density_tmp      = self.get_density_eos(pressure_tmp, gas_constant_tmp, temperature_tmp)
        velocity_x_tmp   = 0.0
        velocity_y_tmp   = 0.0
        velocity_z_tmp   = 0.0

        var_primitiv_bd[0,n_face] = density_tmp
        var_primitiv_bd[1,n_face] = velocity_x_tmp
        var_primitiv_bd[2,n_face] = velocity_y_tmp
        var_primitiv_bd[3,n_face] = velocity_z_tmp
        var_primitiv_bd[4,n_face] = temperature_tmp
        var_primitiv_bd[5,n_face] = pressure_tmp 

      else:
        print( 'No boundary ID', 'N_Face:',n_face, 'N_Cell:',n_cell)
        continue

    return var_primitiv_bd
    

