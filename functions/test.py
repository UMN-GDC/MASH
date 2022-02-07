#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 13:06:04 2022

@author: christian
"""

# args.py
import argparse

p = argparse.ArgumentParser()

# accept two lists of arguments
# like -a 1 2 3 4 -b 1 2 3
p.add_argument('-a', nargs="+", type=int)
p.add_argument('-b', nargs="+", type=int)
args = p.parse_args()

for ai in args.a :
    print(ai)
    