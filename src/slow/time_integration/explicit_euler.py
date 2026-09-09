#!/usr/bin/env python3

# Module to update solution by Euler explicit scheme
# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/31

def explicit_euler(config, geom_dict, metrics_dict, var_dt, var_rhs, var_conserv):

  # Main routine

  # Input parameters
  volume   = metrics_dict['volume_cell']

  # Update conservative variables
  # var_rhs は V*dQ/dt の符号を反転した量なので、引いて進める
  var_conserv -= var_rhs*var_dt/volume

  return var_conserv
