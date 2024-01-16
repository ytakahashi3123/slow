#!/usr/bin/env python3

# Module to update solution in time marching
# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/31

import numpy as np
from orbital.orbital import orbital

def update_solution(config, geom_dict, var_conserv, var_dq):

  # Main routine
  
  # Input parameters
  num_cell = geom_dict['num_cell']

  # Update conservative variables
  for n_cell in range(0,num_cell):
     var_conserv[:,n_cell] = var_conserv[:,n_cell] + var_dq[:,n_cell]

  return var_conserv