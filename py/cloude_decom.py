'''cloude_decom.py

instructions: 

# to run basic example:
    python3 py/cloude_decom.py T3 

# to run basic example faster ( because of using a smaller dataset ):
    python3 py/cloude_decom.py T3_small

# special rgb for targeting:
    python3 py/cloude_decom.py T3 --special_rgb

# specify target location (x,y) coordinate to "cancel" a specific ship:
    python3 py/cloude_decom.py T3  --x=616 --y=393 --no_gui

# cancel another ship:
    python3 py/cloude_decom.py T3  --x=779 --y=718 --no_gui

# when running on a polygon target in shapefile, gui is suppressed ( opt.bin, etc. are produced in the T3 folder )
    python3 py/cloude_decom.py T3 --shapefile=T3/shapefiles/water.shp

Todo:
(DONE) simple gui (point select)
(DONE) polygon select ( shapefile or input to GUI )
(TODO) write files on ESC or close plot

To install dependencies:
    python3 -m pip install matplotlib numpy fiona shapely rasterio

To upgrade pip:
    python3 -m pip install --upgrade pip
'''
import warnings; warnings.filterwarnings("ignore", message="Unable to import Axes3D")
from misc import read_config, read_binary, write_binary, write_hdr, parfor, err, exist
from matplotlib.backend_bases import MouseEvent
from matplotlib.path import Path
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import rasterio
import cProfile
import shapely
import pickle
import ctypes
import shutil
import fiona
import cmath
import copy
import math
import time
import sys
import os
from rasterio.features import geometry_mask
from shapely.geometry import shape
args = sys.argv
sep = os.path.sep

in_dir = None
xp, yp = None, None  # target pixel, "graphics" convention (+x right, +y down)
special_rgb = '--special_rgb' in args  # use special (r,g,b) = (dn, vn, sn) visualization ( alphas ) 
no_gui = '--no_gui' in args  # option to suppress gui window (just run at specific target)
args_new, args_special = [], []

shapefile_mask = None  # mask produced from shapefile provided at command-line

for arg in args:
    if arg[:2] == '--':
        args_special += [arg]
        
        # check for shapefile argument
        w = arg.strip('--').split('=')
        if w[0] == 'shapefile':
            shapefile = w[1]
            print("+r", shapefile)

            shapefile_shapes = None
            try:
                shapefile_shapes = [shape(feature['geometry']) for feature in fiona.open(shapefile)]       
            except:
                err("fiona failed to open shapefile: " + shapefile)
    
            raster_t11 = os.path.normpath(args[1]) + os.path.sep + 'T11.bin'
            try:
                print("+r", raster_t11)
                src = rasterio.open(raster_t11)  # raster_crs = src.crs
                transform, width, height = src.transform, src.width, src.height
    
                shapefile_mask = geometry_mask(shapefile_shapes,
                                               transform=transform,
                                               invert=True,
                                               out_shape=(height, width))

                # convert to image coords (x,y) "graphics convention"
                shapefile_mask = [(y, x) for x, y in zip(*np.where(shapefile_mask))]
            except:
                err("rasterio failed to open raster file: " + raster_t11)

        if w[0] == 'x':
            xp = int(w[1])
        if w[0] == 'y':
            yp = int(w[1])
    else:
        args_new += [arg]
args = args_new


F = None
M_PI = math.pi
nrow, ncol = None, None  # image dimensions
eps = np.finfo(np.float64).eps  # "machine epsilon" for 64-bit floating point number


t11_p, t22_p, t33_p, t12_r_p, t12_i_p, t13_r_p, t13_i_p, t23_r_p, t23_i_p =\
    None, None, None, None, None, None, None, None, None
t11c, t22c, t33c, t12c, t13c, t23c =\
    None, None, None, None, None, None
v1_v, v2_v, v3_v, e1_v, e2_v, e3_v =\
    None, None, None, None, None, None
opt, hv, pwr, sopt, aopt, popt =\
    None, None, None, None, None, None
dn, vn, sn = None, None, None

class vec3:
    def __init__(self, a, b, c):
        self.a, self.b, self.c = a, b, c

    def norm(self):
        return math.sqrt(abs(self.a)*abs(self.a) + abs(self.b)*abs(self.b) + abs(self.c)*abs(self.c));

    def normalize(self):
        norm = self.norm()  # give vector l2 length of 1 
        self.a, self.b, self.c = self.a / norm,\
                                 self.b / norm,\
                                 self.c / norm

    def __truediv__(self, a):  # define division by scalar ( X / a ) operator
        return vec3(self.a / a,
                    self.b / a,
                    self.c / a)
    
    def __sub__(self, b):
        return vec3(self.a - b.b,
                    self.b - b.b,
                    self.c - b.c)
    
    def __str__(self):  # tostring()
        return ' '.join([str(x) for x in [self.a, self.b, self.c]])


def solve_cubic(a, b, c, d):
    _2t13, _2t23, sqrt3 = 2. ** 0.3333333333333333,\
                          2. ** 0.6666666666666666,\
                          math.sqrt(3.)
    
    t1, t2 = b*(-2.*b*b + 9.*a*c) - 27.*a*a*d, 3.*a*c -b*b
    t0 = (t1 +  (4.*(t2*t2*t2) + (t1*t1)) ** 0.5)
    t3 = t0 ** 0.333333333333333333333333
    aX6, bX2, X2 = 6. * a * t3, -2. * b * t3, t3 * t3

    return vec3((bX2 + _2t23*X2 - 2.*_2t13*t2)/aX6,
                (2.*bX2 + _2t13*(2.*(1. + 1j * sqrt3) * t2 + 1j * _2t13*(1j + sqrt3)*X2 ))/(2.*aX6),
                (2.*bX2 + _2t13*(2.*(1. - 1j * sqrt3) * t2 - _2t13*(1.+ 1j * sqrt3)*X2 ))/(2.*aX6))


class herm3:
    def __init__(self, A, B, C, D, E, F):
        self.a, self.b, self.c, self.d, self.e, self.f = A, B, C, D, E, F

    def solve_characteristic(self):
        a, b, c, d, e, f = self.a, self.b, self.c, self.d, self.e, self.f

        _A, _B = -1. + 0j,\
                 a + d + f

        ec = self.e.conjugate()
        cc, bc = self.c.conjugate(), self.b.conjugate()
        _C = (-(a*d) - a*f - d*f + b* bc + c* cc + e* ec)
        _D = d*(a*f - c* cc) + e*(b* cc - a* ec) + bc *(-(b*f) + c* ec)
        
        return solve_cubic(_A, _B, _C, _D)

    def __mul__(self, other):
        if isinstance(other, vec3):
            A, B = self, other
            return vec3(A.a*B.a + A.b*B.b + A.c*B.c,
                        A.d*B.b + A.e*B.c + B.a* (A.b).conjugate(), 
                        A.f*B.c + B.a*(A.c).conjugate() + B.b* (A.e).conjugate())
        else:
            err("herm3 * vec3 operation defined only")
    
    def __str__(self):
        return( ' '.join([str(x) for x in [self.a, self.b, self.c]]) + '\n' +
                ' '.join([str(x) for x in [0, self.d, self.e]]) + '\n' + 
                ' '.join([str(x) for x in [0, 0, self.f]]))
        

def eigv(A, _lambda):  # herm3<cf> &A, cf & lambda){
    ''' syms a lambda b y c z d y e z
        solve( '(a-lambda)+b*y+c*z', 'conj(b) + y*(d-lambda) +e*z')  '''
    Abc = (A.b).conjugate()
    return vec3(1. + 0j, # cf(1.,0.),
                -((A.a)*(A.e)-_lambda*(A.e)-(A.c)* Abc )/((A.b)*(A.e)-(A.d)*(A.c)+_lambda*(A.c)),
                (-(A.b)* Abc -_lambda*(A.a)+(A.d)*(A.a)-(A.d)*_lambda+(_lambda*_lambda))/((A.b)*(A.e)-(A.d)*(A.c)+_lambda*(A.c)))


def eig(A):
    lambdas = A.solve_characteristic()

    e1, e2, e3 = eigv(A, lambdas.a),\
                 eigv(A, lambdas.b),\
                 eigv(A, lambdas.c)

    l1, l2, l3 = lambdas.a, lambdas.b, lambdas.c  

    e1.normalize()
    e2.normalize()
    e3.normalize()

    X = [[abs(l1), e1],
         [abs(l2), e2],
         [abs(l3), e3]]

    X.sort(reverse=True)  # sort eigenvectors by eigenvalue ( decreasing order )
    L = [X[i][0] for i in range(len(X))]

    L = vec3(L[0], L[1], L[2])
    E1, E2, E3 = X[0][1], X[1][1], X[2][1]
    d1, d2, d3 = (A*E1)/(L.a) - E1, (A*E2)/(L.b) - E2, (A*E3)/(L.c) - E3;

    diff = d1.norm() + d2.norm() + d3.norm()
    return [L, E1, E2, E3]


def read_T3(d):
    d = os.path.normpath(d)
    t11f = d + os.path.sep + 'T11.bin'
    t11h = d + os.path.sep + 'T11.hdr'
    if not os.path.exists(t11f):
        err('could not find file: ' + t11f)
    
    global t11_p, t22_p, t33_p, t12_r_p, t12_i_p, t13_r_p, t13_i_p, t23_r_p, t23_i_p 
    sep = os.path.sep
    t11_p, t22_p = read_binary(d + sep + 'T11.bin')[3], read_binary(d + sep + 'T22.bin')[3]
    t33_p = read_binary(d + sep + 'T33.bin')[3]
    t12_r_p, t12_i_p = read_binary(d + sep + 'T12_real.bin')[3], read_binary(d + sep + 'T12_imag.bin')[3]
    t13_r_p, t13_i_p = read_binary(d + sep + 'T13_real.bin')[3], read_binary(d + sep + 'T13_imag.bin')[3]
    t23_r_p, t23_i_p = read_binary(d + sep + 'T23_real.bin')[3], read_binary(d + sep + 'T23_imag.bin')[3]


def lamcloude(a, b=None, c=None, z1=None, z2=None, z3=None):

    if b is None:
        X = a
        a, b, c, z1, z2, z3 = X

    if math.isnan(a.real):
        return([float('nan') for i in range(6)])

    p = 1./3.
    tra = (a + b + c) / 3.
    z1p, z2p, z3p = z1.conjugate(),\
                    z2.conjugate(),\
                    z3.conjugate()
    fac0 = z1 * z1p + z2 * z2p + z3 * z3p

    deta = a * b * c - c * z1 * z1p - b * z2 * z2p + z1 * z2p * z3 + z1p * z2 * z3p - a * z3 * z3p
    s2 = a * a - a * b + b * b - a * c - b * c + c * c + 3. * fac0
    s1 = a * b + a * c + b * c - fac0

    fac1 = 27. * deta - 27. * s1 * tra + 54. * (tra ** 3.)
    tr3 = fac1 + cmath.sqrt((fac1 ** 2.)- 4. * (s2 ** 3.))
    j1s3 = 1j * math.sqrt(3.)
    fac2, fac3 = 1. + j1s3,\
                 1. - j1s3

    p2p = 2. ** p
    ptr3p = tr3 ** p
    p22p = 2. ** ( 2 * p ) 
    e1 = (tra + ptr3p / (3. * p2p) + (s2 * p2p + eps) / (3. * ptr3p + eps)).real
    e2 = (tra - (fac2 * s2) / (3. * ptr3p * p22p + eps) - (fac3 * ptr3p) / (6. * p2p + eps)).real
    e3 = (tra - (fac3 * s2) / (3. * p22p * ptr3p + eps) - (fac2 * ptr3p) / (6. * p2p + eps)).real
    
    if(e1 < e3):  # sort the eigenvalues
        tmp = e1; e1 = e3; e3 = tmp  
    if(e1 < e2):
        tmp = e1; e1 = e2; d2 = tmp 
    if(e2 < e3):
        tmp = e2; e2 = e3; e3 = tmp

    if(e1 < e3 or e1 < e2 or e2 < e3):
        print("Warning: not sorted (%e, %e, %e)\n", e1, e2, e3);

    v2 = ((a - e1) * z3 - z1p * z2) / ((b - e1) * z2 - z3 * z1 + eps)  # dominant eigenvector
    v3 = (e1 - a - z1 * v2) / (z2 + eps)
    v1 = 1.

    av1, av2, av3 = abs(v1), abs(v2), abs(v3)
    n = math.sqrt(av1 * av1 + av2 * av2 + av3 * av3) + eps

    v1, v2, v3 = v1 / n, v2 / n, v3 / n  # normalised components as output
    return [e1, e2, e3, v1, v2, v3] 


def lamcloude_vectorised(a, b, c, z1, z2, z3):
    e1 = np.full(npx, np.nan + 0j) # e1 = np.array([math.nan + 0j for i in range(npx)])
    e2 = np.copy(e1) # e2 = copy.deepcopy(e1)
    e3 = np.copy(e1) # e3 = copy.deepcopy(e1)
    v1 = np.copy(e1) # v1 = copy.deepcopy(e1)
    v2 = np.copy(e1) #  copy.deepcopy(e1)
    v3 = np.copy(e1) #copy.deepcopy(e1)

    print("prepare data..")
    inputs = [[a[i], b[i], c[i], z1[i], z2[i], z3[i]] for i in range(npx)]
    print("run lamcloude..")
    results = parfor(lamcloude, inputs)

    for i in range(npx):
        e1[i] = results[i][0]
        e2[i] = results[i][1]
        e3[i] = results[i][2]
        v1[i] = results[i][3]
        v2[i] = results[i][4]
        v3[i] = results[i][5]
    return [e1, e2, e3, v1, v2, v3]
    '''
    for i in range(npx):
        X = lamcloude(a[i], b[i], c[i], z1[i], z2[i], z3[i])
        e1[i] = X[0]
        e2[i] = X[1]
        e3[i] = X[2]
        v1[i] = X[3]
        v2[i] = X[4]
        v3[i] = X[5]
    return [e1, e2, e3, v1, v2, v3]
    '''

def rank1_t3(e1, v1, v2, v3):  #  generate T3 rank 1
    e1v1, e1v2 = e1 * v1, e1 * v2
    v2c, v3c = v2.conjugate(), v3.conjugate()
    '''
    [t11c, t12c, t13c, t22c, t23c, t33c]   '''
    return [e1v1 * v1.conjugate(), 
            e1v1 * v2c,
            e1v1 * v3c,
            e1v2 * v2c,
            e1v2 * v3c,
            e1 * v3 * v3c]


def rank1_t3_vectorised(e1, v1, v2, v3):  #  generate T3 rank 1 ( numpy vectorised version )
    e1v1, e1v2 = e1 * v1, e1 * v2
    v2c, v3c = np.conjugate(v2), np.conjugate(v3) # v3.conjugate()
    '''
    [t11c, t12c, t13c, t22c, t23c, t33c]   '''
    return [e1v1 * np.conjugate(v1), # v1.conjugate(),
            e1v1 * v2c, 
            e1v1 * v3c, 
            e1v2 * v2c, 
            e1v2 * v3c, 
            e1 * v3 * v3c]


def write_out(variable_name):
    global in_dir
    dd = os.path.normpath(args[1]) + os.path.sep
    cmd =  ('write_binary(' + variable_name + '.tolist(), "' +
            dd + variable_name + '.bin"); write_hdr("' +
            dd + variable_name + '.hdr", ncol, nrow, 1, ["' +
            variable_name + '.bin"])')
    # print(cmd)
    exec(cmd)
    
    # copy the T11.hdr file over the header for this file, so the map info / CRS are retained
    t11h = os.path.normpath(os.path.abspath(in_dir)) + os.path.sep + 'T11.hdr'
    shutil.copy(t11h,
                dd + variable_name + '.hdr')

def decom(o2d1, o2d2, o2d3, o3d1, o3d2, o3d3, o2d1c, o2d2c, o2d3c, o3d1c, o3d2c, o3d3c):   # calculate decom for pixel at linear index "i"
    print("decom..")
    ''' # aliases
      t12c = z1
      t13c = z2
      t23c = z3
      t11c = a
      t22c = b
      t33c = c
    '''
    # project data onto null channels // null_vecs=[o2d o3d];
    z1 = o2d1c * v1_v + o2d2c * v2_v + o2d3c * v3_v
    z2 = o3d1c * v1_v + o3d2c * v2_v + o3d3c * v3_v
    
    # find optimum weights
    npct12c = np.conjugate(t12c)
    popt = np.angle(t13c * npct12c) * 180. / M_PI   # cmath.phase(z2 * z1.conjugate()) * 180. / M_PI
    za = (t12c * npct12c - t13c * np.conjugate(t13c)) + 2j * np.abs(t12c) * np.abs( t13c) #  abs(z1) * abs(z2)
    aopt = np.angle(za) * 90. / M_PI # cmath.phase(za) * 90. / M_PI
    ar = aopt * M_PI / 180.
    br = popt * M_PI / 180.
    
    # optimum weight vector
    npcar = np.cos(ar)
    npsar = np.sin(ar)
    npe1br = np.exp(1j * br)
    w1 = np.conjugate(npcar * o2d1 + npsar * npe1br * o3d1) # .conjugate()
    w2 = np.conjugate(npcar * o2d2 + npsar * npe1br * o3d2) # .conjugate()
    w3 = np.conjugate(npcar * o2d3 + npsar * npe1br * o3d3) # .conjugate()
    
    # find optimum subspace signal
    zopt = w1 * v1_v + w2 * v2_v + w3 * v3_v
    ip = np.abs(zopt * np.conjugate(zopt))  # zopt.conjugate
    ip_eps = ip + eps
    sopt = 10. * np.log(ip_eps) / math.log(10.)  # optimum normalised power
    
    sp = t11c + t22c + t33c  # span power
    abs_sp = np.abs(sp)
    pwr = 10. * np.log(np.abs(sp)) / math.log(10.)  # span channel
    
    sm = np.abs(t33c) #  math.fabs(t33)
    hv = 10. * np.log(sm) / math.log(10.)
    sm2 = sopt + pwr
    
    opt = 10. ** ( sm2 / 10.) # pow(10., sm2 / 10.)  # linear opt channel

    return [ opt, hv, pwr, sopt, aopt, popt]


def nullspace_vectors(xp, yp, mask=None):
    # print('nullspace_vectors', xp, yp)
    t11, t22, t33, t12_r, t12_i, t13_r, t13_i, t23_r, t23_i = 0, 0, 0, 0, 0, 0, 0, 0, 0

    if (xp is None or yp is None) and mask is None:
        err("must specify single point (xp, yp) or list of coords (mask)")

    n_use, n_all = 0, 0
    x_bar, y_bar = 0, 0
    x_bar_all = 0
    y_bar_all = 0
    if mask is None:
        x_bar = xp
        y_bar = yp
        y_bar_all = xp
        y_bar_all = yp
        i = yp * ncol + xp
        if not math.isnan(t11_p[i]):
            t11 = t11_p[i]
            t22 = t22_p[i]
            t33 = t33_p[i]
            t12_r = t12_r_p[i]
            t12_i = t12_i_p[i]
            t13_r = t13_r_p[i]
            t13_i = t13_i_p[i]
            t23_r = t23_r_p[i]
            t23_i = t23_i_p[i]
            n_use += 1
    else:
        # print("mask_pixels", mask)
        print("|mask pixels|=", len(mask))
        for (x,y) in mask:
            i = y * ncol + x
            if not math.isnan(t11_p[i]):
                t11 += t11_p[i]
                t22 += t22_p[i]
                t33 += t33_p[i]
                t12_r += t12_r_p[i]
                t12_i += t12_i_p[i]
                t13_r += t13_r_p[i]
                t13_i += t13_i_p[i]
                t23_r += t23_r_p[i]
                t23_i += t23_i_p[i]
                n_use += 1
                x_bar += x
                y_bar += y
        x_bar_all += x
        y_bar_all += y
        n_all += 1

    print("Target centroid (x,y) before excluding NaN:", x_bar_all/n_use, y_bar_all/n_use)
    if n_use < 1:
        print("Error: no valid data area selected")
        sys.exit(1)
    
    '''(ws > 1){
      for(int di = yp - dw; di <= yp + dw; di++){
        if(di >=0 && di < nrow){
          for(int dj = xp - dw; dj <= xp + dw; dj ++){
            if(dj >=0 && dj < ncol){
              int j = di * ncol + dj;
    
              t11 += (double)t11_p[j]; t22 += (double)t22_p[j];
              t33 += (double)t33_p[j];
              t12_r += (double)t12_r_p[j]; t12_i += (double)t12_i_p[j];
              t13_r += (double)t13_r_p[j]; t13_i += (double)t13_i_p[j];
              t23_r += (double)t23_r_p[j]; t23_i += (double)t23_i_p[j];
    
              n_use ++;
            }
          }
        }
      }
    }
    '''
    # printf("N_USE: %f\n", n_use);
    t11 /= n_use; t22 /= n_use; t33 /= n_use
    t12_r /= n_use; t12_i /= n_use
    t13_r /= n_use; t13_i /= n_use
    t23_r /= n_use; t23_i /= n_use
    print("Target T11:", t11)
    print("Target T22:", t22)
    print("Target T33:", t33)
    x_bar /= n_use; y_bar /= n_use;
    print("Target centroid (x,y)=", x_bar, y_bar)
    
    a = t11; b = t22; c = t33
    z1 = t12_r + t12_i * 1j
    z2 = t13_r + t13_i * 1j
    z3 = t23_r + t23_i * 1j
    
    # /* avoid 0 elements.. conditioning */
    eps2 = (a + b + c) * (1.0e-9) + eps
    F = (float(nrow + ncol) / 2.) + 0j
    z1 = z1 + eps2 * F
    z2 = z2 + eps2 * F
    z3 = z3 + eps2 * F
    a = a + eps2 * F
    b = b + eps2 * F
    c = c + eps2 * F
    
    T = herm3(a, z1, z2, b, z3, c)
    
    L, E1, E2, E3 = eig(T)
    '''print("L", L)
    print("E2", E2)
    print("E3", E3)'''
    
    o2d1 = E2.a
    o2d2 = E2.b
    o2d3 = E2.c
    
    o3d1 = E3.a
    o3d2 = E3.b
    o3d3 = E3.c
    
    o2d1c, o2d2c, o2d3c, o3d1c, o3d2c, o3d3c = o2d1.conjugate(), o2d2.conjugate(),\
                                               o2d3.conjugate(), o3d1.conjugate(),\
                                               o3d2.conjugate(), o3d3.conjugate()
    return [o2d1, o2d2, o2d3, o3d1, o3d2, o3d3,
            o2d1c, o2d2c, o2d3c, o3d1c, o3d2c, o3d3c]


# program start
in_dir = os.path.normpath(args[1])
x = read_config(in_dir + os.path.sep + 'config.txt')
nrow, ncol = x['nrow'], x['ncol']
F = (float(nrow) + float(ncol)) / 2.
npx = nrow * ncol
read_T3(in_dir)  # load the T3 matrix data

# initialize output variables
# out_r, out_g, out_b, out_e1, out_e2, out_e3, out_opt, out_sopt, out_v1 = [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)]


# xp, yp = int(args[3]) if len(args) > 3 else None,\
#         int(args[2]) if len(args) > 2 else None    # cloude_decom T3 313 798 # col/ row for previous test target

# check specified col, row indices ( respectively ) are "in bounds"
if xp is not None and (xp < 0 or xp > ncol):
        err("x coord out of bounds")
if yp is not None and (yp < 0 or yp > nrow):
        err("y coord out of bounds")

print("complex arrays..")
N = len(t11_p)
t11c = np.array(t11_p) + 0j
t22c = np.array(t22_p) + 0j
t33c = np.array(t33_p) + 0j
t12c = np.array(t12_r_p) + np.array(t12_i_p) * 1j
t13c = np.array(t13_r_p) + np.array(t13_i_p) * 1j
t23c = np.array(t23_r_p) + np.array(t23_i_p) * 1j
''' # aliases
t12c = z1
t13c = z2
t23c = z3
t11c = a
t22c = b
t33c = c
'''
# /* avoid 0 elements.. conditioning */
print("conditioning..")
eps2 = (t11c + t22c + t33c ) * (1.0e-9) + eps
t12c += eps2 * F 
t13c += eps2 * F 
t23c += eps2 * F
t11c += eps2 * F 
t22c += eps2 * F
t33c += eps2 * F 

pickle_filename = os.path.normpath(args[1]) + os.path.sep + 'cloude_decom.pkl'
if os.path.exists(pickle_filename):
    print("unpickling..")
    [t11c, t22c, t33c, t12c, t13c, t23c, v1_v, v2_v, v3_v, e1_v, e2_v, e3_v, dn, vn, sn] =\
        pickle.load(open(pickle_filename, 'rb'))
else:
    print("lamcloude..")
    [e1_v, e2_v, e3_V, v1_v, v2_v, v3_v] =\
        lamcloude_vectorised(t11c, t22c, t33c, t12c, t13c, t23c)
    
    print("rank1 t3..")
    [t11c, t12c, t13c, t22c, t23c, t33c] =\
        rank1_t3_vectorised(e1_v, v1_v, v2_v, v3_v)
   
    ''' 
    T11c = np.abs(t11c)
    T22c = np.abs(t22c)
    T33c = np.abs(t33c)

    for x in ['T11c', 'T22c', 'T33c']:
        write_out(x)

    '''

    # generate alpha etc. eigenvector parameters : ) 
    alpha = np.arccos(np.abs(v1_v))
    phi = np.angle(t12c)
    theta = np.angle((t22c - t33c) + 2. * 1j * t23c.real) / 4.  

    # generate RGB colour composite from multiple eigenvector angles
    dn = alpha * 2. / M_PI # alpha angle in red channel                # out: red channel
    theta2 = theta + (theta > M_PI / 4.) * (M_PI / 4. - theta)
    theta2 = theta2 + (theta2 < -M_PI / 4.) * (-M_PI / 4. - theta2)
    vn = (theta2 + M_PI / 4.) * 2. / M_PI   # az slope is green        # out: green channel
    sn = np.abs(phi) / M_PI  # mag of Pauli phase is blue (180 is Bragg)  # out: blue channel

    # special RGB encoding: (r,g,b) = (dn, vn, sn)
    for x in ['alpha', 'phi', 'theta', 'dn', 'theta2', 'vn', 'sn']:
        write_out(x)

    print("pickling..")
    pickle.dump([t11c, t22c, t33c,
                 t12c, t13c, t23c,
                 v1_v, v2_v, v3_v,
                 e1_v, e2_v, e3_v,
                 dn, vn, sn],
                open(pickle_filename, 'wb'))


print("null vectors..")
if (xp is not None and yp is not None) or shapefile_mask is not None:
    [o2d1, o2d2, o2d3, 
     o3d1, o3d2, o3d3,
     o2d1c, o2d2c, o2d3c,
     o3d1c, o3d2c, o3d3c] = nullspace_vectors(xp if shapefile_mask is None else xp,
                                              yp if shapefile_mask is None else yp,
                                              shapefile_mask)

    [opt, hv, pwr, sopt, aopt, popt]  =\
        decom(o2d1, o2d2, o2d3,
              o3d1, o3d2, o3d3,
              o2d1c, o2d2c, o2d3c,
              o3d1c, o3d2c, o3d3c)

    for x in ['opt']: # 'hv', 'pwr', 'sopt', 'aopt', 'popt']:
        write_out(x)

    if no_gui or (shapefile_mask is not None):
        sys.exit(0)

    # shapefile option goes here:


def naninf_list(x):
    Y = []
    X = list(x.ravel().tolist()) if type(x) != list else x
    for i in X:
        if not (math.isnan(i) or math.isinf(i)):
            Y.append(i)
    return Y

def scale(rgb_i):
    values = naninf_list(rgb_i) # values.reshape(np.prod(values.shape)).tolist()
    values.sort()
    n_pct = 1. # percent for stretch value
    frac = n_pct / 100.
    rgb_min, rgb_max = values[int(math.floor(float(len(values))*frac))],\
                       values[int(math.floor(float(len(values))*(1. - frac)))]
    rng = rgb_max - rgb_min  # apply restored or derived scaling
    rgb_i = (rgb_i - rgb_min) / (rng if rng != 0. else 1.)
    rgb_i[rgb_i < 0.] = 0.  # clip
    rgb_i[rgb_i > 1.] = 1.
    return rgb_i


# default visualization
print("scaling rgb..")
rgb, rgb_pkl = None, None
if special_rgb:
    rgb_pkl = os.path.normpath(args[1]) + os.path.sep + 'rgb_special.pkl'
    if exist(rgb_pkl):
        [rgb] = pickle.load(open(rgb_pkl), 'rb')
    else:
        rgb = np.zeros((nrow, ncol, 3))
        rgb[:, :, 0] = scale(np.array(dn)).reshape((nrow, ncol))
        rgb[:, :, 1] = scale(np.array(vn)).reshape((nrow, ncol))
        rgb[:, :, 2] = scale(np.array(sn)).reshape((nrow, ncol))
        pickle.dump([rgb],
                    open(rgb_pkl, 'wb'))
else:
    rgb_pkl = os.path.normpath(args[1]) + os.path.sep + 'rgb.pkl'
    if exist(rgb_pkl):
        [rgb] = pickle.load(open(rgb_pkl, 'rb'))
    else:
        rgb = np.zeros((nrow, ncol, 3))
        rgb[:, :, 0] = scale(np.array(t22_p)).reshape((nrow, ncol))
        rgb[:, :, 1] = scale(np.array(t33_p)).reshape((nrow, ncol))
        rgb[:, :, 2] = scale(np.array(t11_p)).reshape((nrow, ncol))
        pickle.dump([rgb],
                    open(rgb_pkl, 'wb'))
print("done scaling.")

# gui stuff
vertices = []  # (x,y) coordinates in our polygon ( once finalized )
drawing_poly = False  # we are in the process of drawing a polygon
decom_plotted = False
fig, ax = plt.subplots()

im = ax.imshow(rgb)
# plt.title('JAXA ALOS-1 data over SanFransisco, USA')
plt.xlabel('(R,G,B)=(T22, T33, T11)')
plt.tight_layout()
line, = ax.plot([], [], 'r-', lw=2)  # dynamically updated line
polygon_lines = []

def data_to_pixel(x, y):
    pixel_coords = ((x + 0.5, y + 0.5))
    return np.floor(pixel_coords).astype(int)


def get_polygon_mask(polygon_coords, image_width, image_height):
    """
    Create binary mask: points under poly have value True
    :param polygon_coords: List of (x, y) tuples representing the polygon's vertices
    :param image_width: Width of image (# of columns)
    :param image_height: Height of image (# of rows)
    :return: A binary mask (2D array): True == point is under poly
    """
    polygon = Path(polygon_coords)
    x, y = np.meshgrid(np.arange(image_width), np.arange(image_height))
    points = np.vstack((x.ravel(), y.ravel())).T  # Shape: (num_points, 2)
    inside_mask = polygon.contains_points(points)
    inside_image = inside_mask.reshape(image_height, image_width)
    return inside_image

# update dynamic line segment
def update_line(event):
    global line
    if len(vertices) > 0:
        x, y = data_to_pixel(event.xdata, event.ydata)  # current mouse position
        line.set_data([vertices[-1][0], x],
                      [vertices[-1][1], y])  # update the line
        fig.canvas.draw()

def on_press(event):  # called when point is clicked
    global decom_plotted
    global drawing_poly  # state variable
    global vertices, polygon_lines, line

    # print("on_press(): now release mouse button over target location")
    if event.button == 1:
        if not drawing_poly:  # left mouse click
            print("imshow(rgb) 1")
            plt.cla()
            ax.imshow(rgb)
            plt.xlabel('(R,G,B)=(T22, T33, T11)')
            plt.draw()  # Redraw the canvas
            decom_plotted = False
        else:
            if len(vertices) > 1:
                new_line, = ax.plot([vertices[-1][0], vertices[0][0]],
                                    [vertices[-1][1], vertices[0][1]], 'g-', lw=2)
                polygon_lines += [new_line]

                line.set_data([], [])  # clear the red line data
                fig.canvas.draw()
                print("Polygon completed!")
                print(vertices)
                mask = get_polygon_mask(vertices, ncol, nrow)
                y_ind, x_ind = np.where(mask)  # row/col indices under mask
                mask_points = list(zip(x_ind, y_ind))  
                print(str(len(mask_points)), 'points under mask')
                print(mask_points)

                [o2d1, o2d2, o2d3, 
                 o3d1, o3d2, o3d3,
                 o2d1c, o2d2c, o2d3c,
                 o3d1c, o3d2c, o3d3c] = nullspace_vectors(None, None, mask_points)

                global opt, hv, pwr, sopt, aopt, popt
                [opt, hv, pwr, sopt, aopt, popt]  =\
                    decom(o2d1, o2d2, o2d3,
                          o3d1, o3d2, o3d3,
                          o2d1c, o2d2c, o2d3c,
                          o3d1c, o3d2c, o3d3c)

                for x in ['opt']: # , 'hv', 'pwr', 'sopt', 'aopt', 'popt']:
                    write_out(x)

                print("imshow(opt) poly")
                ax.imshow(scale(opt).reshape((nrow, ncol)),
                          cmap='gray',
                          vmin=0,
                          vmax=1)  # Update the image
                plt.xlabel('opt.bin')
                plt.draw()
                decom_plotted = True

            vertices = []  # Reset vertices after finalizing polygon
            # drawing_poly = False

    else:
        pass  # no action on middle button

def on_release(event):
    global opt
    global hv
    global pwr
    global sopt
    global aopt
    global popt
    global drawing_poly
    global decom_plotted
    global vertices, polygon_lines, line
    x, y = event.xdata, event.ydata

    if event.button == 1 and drawing_poly:
        drawing_poly = False
        return

    if event.button == 1 and (x is not None) and (y is not None):  # ensure click within axes
        x = math.floor(x + 0.5)
        y = math.floor(y + 0.5)
        print('on_release(): run decom at (x,y)=' + str(x) + ',' + str(y))

        [o2d1, o2d2, o2d3,
         o3d1, o3d2, o3d3,
         o2d1c, o2d2c, o2d3c,
         o3d1c, o3d2c, o3d3c] = nullspace_vectors(x, y)

        [opt, hv, pwr, sopt, aopt, popt]  = decom(o2d1, o2d2, o2d3,
                                                  o3d1, o3d2, o3d3,
                                                  o2d1c, o2d2c, o2d3c,
                                                  o3d1c, o3d2c, o3d3c)

        for x in ['opt']: # 'hv', 'pwr', 'sopt', 'aopt', 'popt']:
            write_out(x)
    
        print("imshow(opt) 2")
        plt.cla()
        ax.imshow(scale(opt).reshape((nrow, ncol)),
                  cmap='gray',
                  vmin=0,
                  vmax=1)  # Update the image
        plt.xlabel('opt.bin')
        plt.draw()
        decom_plotted = True

    elif event.button == 3:  # Right mouse btn.. start or continue drawing poly
        if event.xdata is not None and event.ydata is not None:
            if not drawing_poly:
                vertices = []
                for line in polygon_lines:
                    line.remove()
                polygon_lines = []
                line, = ax.plot([], [], 'r-', lw=2)

            if decom_plotted: # or not drawing_poly:
                print("imshow(rgb) 2")
                ax.imshow(rgb)
                plt.xlabel('(R,G,B)=(T22, T33, T11)')
                decom_plotted = False

            drawing_poly = True
            vertices.append(data_to_pixel(event.xdata,
                                          event.ydata))  # add vertex

            if len(vertices) > 1:
                new_line, = ax.plot([vertices[-2][0], vertices[-1][0]],
                                    [vertices[-2][1], vertices[-1][1]], 'g-', lw=2)
                polygon_lines += [new_line]
            
            if len(vertices) > 1:  # update line if >1 points
                update_line(event)
            else:
                fig.canvas.draw()

def quit():
    # what to do at the end
    global opt
    global hv
    global pwr
    global sopt
    global aopt
    global popt
 
    for x in ['opt', 'hv', 'pwr', 'sopt', 'aopt', 'popt']:
        write_out(x)
    print('Cheerio')
    sys.exit(0)

def on_key(event):
    if event.key == 'escape':
        quit()
 
def on_close(event):
    quit()

fig.canvas.mpl_connect('button_press_event', on_press)
fig.canvas.mpl_connect('button_release_event', on_release)
fig.canvas.mpl_connect('motion_notify_event', update_line)

fig.canvas.mpl_connect('key_press_event', on_key)
fig.canvas.mpl_connect('close_event', on_close)
plt.show()
