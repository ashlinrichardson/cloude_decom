'''cloude_decom.py

Todo:
simple gui (point select)
polygon select ( shapefile or input to GUI )

e.g.

python3 cloude_decom.py  313 798
'''
from misc import read_config, read_binary, write_binary, write_hdr
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import cProfile
import pickle
import ctypes
import cmath
import copy
import math
import time
import sys
import os
args = sys.argv

special_rgb = '--special_rgb' in args
args_new = []

for arg in args:
    if arg[:2] == '--':
        pass
    else:
        args_new += [arg]
args = args_new


F = None
M_PI = math.pi
xp, yp = None, None  # target pixel, "graphics" convention (+x right, +y down)
nrow, ncol = None, None  # image dimensions
eps = np.finfo(np.float64).eps  # "machine epsilon" for 64-bit floating point number


t11_p, t22_p, t33_p, t12_r_p, t12_i_p, t13_r_p, t13_i_p, t23_r_p, t23_i_p =\
    None, None, None, None, None, None, None, None, None
t11c, t22c, t33c, t12c, t13c, t23c =\
    None, None, None, None, None, None
v1_v, v2_v, v3_v, e1_v, e2_v, e3_v =\
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


def lamcloude(a, b, c, z1, z2, z3):
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
    e1 = np.array([math.nan + 0j for i in range(npx)])
    e2 = copy.deepcopy(e1)
    e3 = copy.deepcopy(e1)
    v1 = copy.deepcopy(e1)
    v2 = copy.deepcopy(e1)
    v3 = copy.deepcopy(e1)

    for i in range(npx):
        X = lamcloude(a[i], b[i], c[i], z1[i], z2[i], z3[i])
        e1[i] = X[0]
        e2[i] = X[1]
        e3[i] = X[2]
        v1[i] = X[3]
        v2[i] = X[4]
        v3[i] = X[5]
    return [e1, e2, e3, v1, v2, v3]


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
    cmd =  ('write_binary(' + variable_name + '.tolist(), "' +
            variable_name + '.bin"); write_hdr("' +
            variable_name + '.hdr", ncol, nrow, 1, ["' +
            variable_name + '.bin"])'
    # print(cmd)
    exec(cmd)


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
    popt = np.angle(t13c * np.conjugate(t12c)) * 180. / M_PI   # cmath.phase(z2 * z1.conjugate()) * 180. / M_PI
    za = (t12c * np.conjugate(t12c) - t13c * np.conjugate(t13c)) + 1j * 2. * np.abs(t12c) * np.abs( t13c) #  abs(z1) * abs(z2)
    aopt = np.angle(za) * 90. / M_PI # cmath.phase(za) * 90. / M_PI
    ar = aopt * M_PI / 180.
    br = popt * M_PI / 180.
    
    print(o2d1)
    print(o3d1)
    # optimum weight vector
    w1 = np.conjugate(np.cos(ar) * o2d1 + np.sin(ar) * np.exp(1j * br) * o3d1) # .conjugate()
    w2 = np.conjugate(np.cos(ar) * o2d2 + np.sin(ar) * np.exp(1j * br) * o3d2) # .conjugate()
    w3 = np.conjugate(np.cos(ar) * o2d3 + np.sin(ar) * np.exp(1j * br) * o3d3) # .conjugate()
    
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


def nullspace_vectors(xp, yp):
    print('nullspace_vectors', xp, yp)
    n_use = 1;
    i = yp * ncol + xp
    t11 = t11_p[i]
    t22 = t22_p[i]
    t33 = t33_p[i]
    t12_r = t12_r_p[i]
    t12_i = t12_i_p[i]
    t13_r = t13_r_p[i]
    t13_i = t13_i_p[i]
    t23_r = t23_r_p[i]
    t23_i = t23_i_p[i]
    
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

xp, yp = int(args[3]) if len(args) > 3 else None,\
         int(args[2]) if len(args) > 2 else None    # cloude_decom T3 313 798 # col/ row for previous test target

# check specified col, row indices ( respectively ) are "in bounds"
if xp is not None and (xp < 0 or xp > ncol):
        err("x coord out of bounds")
if yp is not None and (yp < 0 or yp > nrow):
        err("y coord out of bounds")

print("complex arrays..")
N = len(t11_p)
t11c = [t11_p[i] + 0j for i in range(N)] # + 0j; # a = t11 + 0j
t22c = [t22_p[i] + 0j for i in range(N)] # + 0j; # b = t22 + 0j
t33c = [t33_p[i] + 0j for i in range(N)] # + 0j; # c = t33 + 0j

t12c = [t12_r_p[i] + t12_i_p[i] * 1j  for i in range(N)] # z1 = t12_r + t12_i * 1j
t13c = [t13_r_p[i] + t13_i_p[i] * 1j  for i in range(N)] # z2 = t13_r + t13_i * 1j
t23c = [t23_r_p[i] + t23_i_p[i] * 1j  for i in range(N)] # z3 = t23_r + t23_i * 1j

''' # aliases
t12c = z1
t13c = z2
t23c = z3
t11c = a
t22c = b
t33c = c
'''
t11c = np.array(t11c)
t22c = np.array(t22c) 
t33c = np.array(t33c)

# /* avoid 0 elements.. conditioning */
print("conditioning..")
eps2 = (t11c + t22c + t33c ) * (1.0e-9) + eps
t12c = t12c + eps2 * F 
t13c = t13c + eps2 * F 
t23c = t23c + eps2 * F
t11c = t11c + eps2 * F 
t22c = t22c + eps2 * F
t33c = t33c + eps2 * F 

pickle_filename = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + 'cloude_decom.pkl'
if os.path.exists(pickle_filename):
    print("unpickling..")
    [t11c, t22c, t33c, t12c, t13c, t23c, v1_v, v2_v, v3_v, e1_v, e2_v, e3_v] =\
        pickle.load(open(pickle_filename, 'rb'))
else:
    print("lamcloude..")
    [e1_v, e2_v, e3_V, v1_v, v2_v, v3_v] =\
        lamcloude_vectorised(t11c, t22c, t33c, t12c, t13c, t23c)
    
    print("rank1 t3..")
    [t11c, t12c, t13c, t22c, t23c, t33c] =\
        rank1_t3_vectorised(e1_v, v1_v, v2_v, v3_v)
    
    T11c = np.abs(t11c)
    T22c = np.abs(t22c)
    T33c = np.abs(t33c)

    for x in ['T11c', 'T22c', 'T33c']:
        write_out(x)

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
                 e1_v, e2_v, e3_v],
                open(pickle_filename, 'wb'))


print("null vectors..")
if xp is not None and yp is not None:
    [o2d1, o2d2, o2d3, 
     o3d1, o3d2, o3d3,
     o2d1c, o2d2c, o2d3c,
     o3d1c, o3d2c, o3d3c] = nullspace_vectors(xp, yp)

    [opt, hv, pwr, sopt, aopt, popt]  =\
        decom(o2d1, o2d2, o2d3,
              o3d1, o3d2, o3d3,
              o2d1c, o2d2c, o2d3c,
              o3d1c, o3d2c, o3d3c)

    for x in ['opt', 'hv', 'pwr', 'sopt', 'aopt', 'popt']:
        write_out(x)


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
rgb = np.zeros((nrow, ncol, 3))
if special_rgb:
    rgb[:, :, 0] = scale(np.array(dn)).reshape((nrow, ncol))
    rgb[:, :, 1] = scale(np.array(vn)).reshape((nrow, ncol))
    rgb[:, :, 2] = scale(np.array(sn)).reshape((nrow, ncol))
else:
    rgb[:, :, 0] = scale(np.array(t22_p)).reshape((nrow, ncol))
    rgb[:, :, 1] = scale(np.array(t33_p)).reshape((nrow, ncol))
    rgb[:, :, 2] = scale(np.array(t11_p)).reshape((nrow, ncol))
fig, ax = plt.subplots()
im = ax.imshow(rgb)
# plt.title('JAXA ALOS-1 data over SanFransisco, USA')
plt.xlabel('(R,G,B)=(T22, T33, T11)')
plt.tight_layout()


def on_press(event):  # called when point is clicked
    print("on_press(): now release mouse button over target location")
    ax.imshow(rgb)
    plt.xlabel('(R,G,B)=(T22, T33, T11)')
    plt.draw()  # Redraw the canvas


def on_release(event):
    print("on_release", event.xdata, event.ydata)
    x, y = event.xdata, event.ydata

    if x is not None and y is not None:  # ensure click within axes
        x = math.floor(x + 0.5)
        y = math.floor(y + 0.5)

        [o2d1, o2d2, o2d3,
         o3d1, o3d2, o3d3,
         o2d1c, o2d2c, o2d3c,
         o3d1c, o3d2c, o3d3c] = nullspace_vectors(x, y)

        [opt, hv, pwr, sopt, aopt, popt]  = decom(o2d1, o2d2, o2d3,
                                                  o3d1, o3d2, o3d3,
                                                  o2d1c, o2d2c, o2d3c,
                                                  o3d1c, o3d2c, o3d3c)
        ax.imshow(scale(opt).reshape((nrow, ncol)),
                  cmap='gray',
                  vmin=0,
                  vmax=1)  # Update the image
        plt.xlabel('opt.bin')
        plt.draw()


fig.canvas.mpl_connect('button_press_event', on_press)
fig.canvas.mpl_connect('button_release_event', on_release)
plt.show()
