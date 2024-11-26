'''misc.py: helper functions for reading data binary files and headers
'''
import sys

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
    # print('+r', hdr)
    for line in open(hdr).readlines():
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
