#!/usr/bin/python3
import os
import sys

def err(m):
    print("Error: " + m); sys.exit(1)

def run(cmd):
    print(cmd); return os.system(cmd)

def exist(p): return os.path.exists(p)

p = '/usr/local/bin/'
f = ['cloude_decom', 'cloude_view']

for i in f:
    x = p + i
    if exist(x):
        run('sudo rm -f ' + x)
