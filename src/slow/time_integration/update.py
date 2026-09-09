#!/usr/bin/env python3

# Module to update solution in time marching
# Author: Y.Takahashi, Hokkaido University
# Date; 2022/03/31

def update_solution(config, geom_dict, var_conserv, var_dq):

  # Main routine

  # Update conservative variables
  # 全セル・全成分について同じ足し算なので、セルループにする理由が無い
  var_conserv += var_dq

  return var_conserv
