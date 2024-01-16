#!/bin/bash

LD=../src/slow.py
LOG=log_slow

#$LD 2>&1 | tee $LOG
#python3 $LD > $LOG
python3.9 $LD > $LOG
