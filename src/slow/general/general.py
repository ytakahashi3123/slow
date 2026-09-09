#!/usr/bin/env python3

# Author: Y.Takahashi, Hokkaido University
# Date: 2022/03/31

import logging

logger = logging.getLogger(__name__)


class general:

  def __init__(self):
    logger.info('Calling class: general')

# FUnctions
  def argument(self, filename_default):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', action='store', type=str, default=filename_default)
    args = parser.parse_args()
    return args


  def read_config_yaml(self, file_control):

    import yaml as yaml
    import sys as sys
    #import pprint as pprint

    logger.info('Reading control file...: %s', file_control)

    try:
      with open(file_control) as file:
        config = yaml.safe_load(file)
#        pprint.pprint(config)
    except Exception as e:
      logger.error('Exception occurred while loading YAML...')
      logger.error('%s', e)
      sys.exit(1)

    return config


  def make_directory(self, dir_path):
  
    import os as os
    import shutil as shutil

    if not os.path.exists(dir_path):
      os.mkdir(dir_path)
    
    return


  def make_directory_rm(self, dir_path):
  
    import os as os
    import shutil as shutil

    if not os.path.exists(dir_path):
      os.mkdir(dir_path)
    else:
      shutil.rmtree(dir_path)
      os.mkdir(dir_path)

    return
    

  def check_file_exist(self, dir_path):
  
    import os as os

    if os.path.exists(dir_path):
      flag_file_exist = True
    else: 
      flag_file_exist = False

    return flag_file_exist


  def split_file(self, filename,addfile, splitchar):
    """
    特定の文字の前に'_***'を加える.
    特定文字列が２つ以上ある場合は未対応
    """
    splitchar_tmp   = splitchar
    filename_split  = filename.rsplit(splitchar_tmp, 1)
    filename_result = filename_split[0]+addfile+splitchar+filename_split[1]

    return filename_result
    

  def insert_suffix(self, filename, suffix, splitchar):
    parts = filename.split(splitchar)
    if len(parts) == 2:
      new_filename = f"{parts[0]}{suffix}.{parts[1]}"
      return new_filename
    else:
      # ファイル名が拡張子を含まない場合の処理
      return filename + suffix


  def get_file_extension(self, filename):
    # ドットで分割し、最後の要素が拡張子となる
    parts = filename.split(".")
    if len(parts) > 1:
      return parts[-1].lower()
    else:
      # ドットが含まれていない場合は拡張子が存在しない
      return None


  def zero_pad_number(self, number, length=3):
    """
    Zero-pad the given number to the specified length.

    Parameters:
    - number: The number to be zero-padded.
    - length: The desired length of the resulting string, including zero-padding.

    Returns:
    - Zero-padded string representation of the number.
    """
    return f"{number:0{length}d}"

    