#!/usr/bin/python3
'''tested on Ubuntu v. 20. Support to be added for: other linux, Mac OS, Windows'''
import os
import sys

def err(m):
    print("Error: " + m)
    sys.exit(1)

def run(cmd):
    print(cmd)
    return os.system(cmd)

# ubuntu 20 dependencies
run('sudo apt install g++ freeglut3-dev')

def build(n):
    run('rm -f ' + n) 
    run('./make_' + n + '.sh')
    run('chmod 755 ' + n)
    run('sudo mv ' + n + ' /usr/bin/' + n)

for n in ['cloude_decom', 'cloude_view']:
    build(n)
