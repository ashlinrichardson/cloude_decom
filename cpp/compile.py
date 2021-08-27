#!/usr/bin/python3
'''tested on Ubuntu v. 20 and MacOS. Support to be added for: Windows'''
import os
import sys
import platform

def err(m):
    print("Error: " + m); sys.exit(1)

def run(cmd):
    print(cmd); return os.system(cmd)

def exist(p): return os.path.exists(p)

def no_cpp():  # returns True if g++ not available
    return len(os.popen('g++ 2>&1').read().strip().split('no input')) < 1

''' detect operating system '''
WIN = os.name == 'nt'
MAC = platform.system() == 'Darwin'
LNX = platform.system() == 'Linux'

if LNX:  # ubuntu dependencies
    if no_cpp():
        run('sudo apt install g++')
    if not exist('/usr/include/GL/glut.h'):
        run('sudo apt install freeglut3-dev')

if MAC and no_cpp():
    err('please install xcode and try again.. (e.g google "how to install xcode on mac os"')

def build(n):
    if not exist('/usr/local/bin'):
        run('sudo mkdir -p -m 775 /usr/local/bin')
    run('rm -f cpp/' + n)
    run('cd cpp; ./make_' + n + '.sh')
    run('chmod 755 cpp/' + n)
    run('sudo mv cpp/' + n + ' /usr/local/bin/' + n)  # because /usr/bin has "integrity protection" on mac

# build both programs
for n in ['cloude_decom', 'cloude_view']:
    build(n)
