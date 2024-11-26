'''misc.py: helper functions for reading data binary files and headers
'''
import numpy as np
import sys
import os

def err(m):
    print('Error: ', m)
    sys.exit(1)

'''read PolSARPro config.txt file to get image dimensions!
'''
def read_config(cfg_file):
    lines = [x.strip().lower() for x in open(cfg_file).readlines()]
    
    if lines[0] != 'nrow':
        err("line 1 != 'Nrow'")
    if lines[3] != 'ncol':
        err("line 4 != 'Ncol'") 
    if lines[6] != 'polarcase':
        err("line 7 != 'PolarCase'")
    if lines[9] != 'polartype':
        err("line 10 != 'PolarType'")
    nrow = int(lines[1])
    ncol = int(lines[4]) 
    
    return {'nrow': nrow, 'ncol': ncol}

'''read ENVI-format header file to get image dimensions!
'''
def read_hdr(hdr):
    samples, lines, bands = 0, 0, 0
    print('+r', hdr)
    for line in [x.strip() for x in open(hdr).readlines()]:
        print(line)
        line = line.strip()
        words = line.split('=')
        if len(words) == 2:
            f, g = words[0].strip(), words[1].strip()
            if f == 'samples':
                samples = g
            if f == 'lines':
                lines = g
            if f == 'bands':
                bands = g
    return {'ncol': int(samples), 'nrow': int(lines), 'nband': int(bands)}


def exist(f):
    return os.path.exists(f)


def hdr_fn(bin_fn):  # return filename for hdr file, given binfile name
    hfn = bin_fn[:-4] + '.hdr'
    if not exist(hfn):
        hfn2 = bin_fn + '.hdr'
        if not exist(hfn2):
            err("header not found at:" + hfn + " or: " + hfn2)
        return hfn2
    return hfn


# use numpy to read a floating-point data file (4 bytes per float, byte order 0)
def read_float(fn):
    print("+r", fn)
    return np.fromfile(fn, dtype = np.float32) # "float32") # '<f4')


def read_binary(fn):
    hdr = hdr_fn(fn) # read header and print parameters
    x = read_hdr(hdr)
    samples, lines, bands = x['ncol'], x['nrow'], x['nband']
    samples, lines, bands = int(samples), int(lines), int(bands)
    print("\tsamples", samples, "lines", lines, "bands", bands)
    data = read_float(fn)
    return samples, lines, bands, data


if __name__ == '__main__':
    x = read_config('../T3/config.txt')
    y = read_hdr('../T3/T11.hdr')
    
    if x['nrow'] != y['nrow'] != 624:
        err('unexpected number of rows')

    if x['ncol'] != y['ncol'] != 1920:
        err('unexpected number of cols')
    
    if y['nband'] != 1:
        err('unexpected number of bands')

    print('test passed')
