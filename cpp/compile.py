#!/usr/bin/python3
'''tested on Ubuntu v. 20. Support to be added for: other linux, Mac OS, Windows'''
import os
import sys
import platform

def err(m):
    print("Error: " + m); sys.exit(1)

def run(cmd):
    print(cmd); return os.system(cmd)

''' detect operating system '''
WIN = os.name == 'nt'
MAC = platform.system() == 'Darwin'
LNX = platform.system() == 'Linux'

# ubuntu 20 dependencies
if len(os.popen('g++ 2>&1').read().strip().split('fatal')) < 1:
    if LNX:
        run('sudo apt install g++')
if not os.path.exists('/usr/include/GL/glut.h'):
    if LNX:
        run('sudo apt install freeglut3-dev')

def build(n):
    run('rm -f cpp/' + n)
    run('cd cpp; ./make_' + n + '.sh')
    run('chmod 755 cpp/' + n)
    run('sudo mv cpp/' + n + ' /usr/local/bin/' + n)  # because /usr/bin has "integrity protection" on mac

# build both programs
for n in ['cloude_decom', 'cloude_view']:
    build(n)
