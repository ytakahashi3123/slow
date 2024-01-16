#!/usr/bin/env python3

# Module to update solution by Euler explicit scheme
# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/31

import numpy as np
from orbital.orbital import orbital

def explicit_euler(config, geom_dict, metrics_dict, var_dt, var_rhs, var_conserv):

  # Main routine
  
  # Input parameters
  num_cell = geom_dict['num_cell']
  volume   = metrics_dict['volume_cell']

  # Update conservative variables
  for n_cell in range(0,num_cell):
    var_conserv[:,n_cell] = var_conserv[:,n_cell] - var_rhs[:,n_cell]*var_dt[n_cell]/volume[n_cell]

  return var_conserv
