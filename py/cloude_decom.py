from misc import read_config, read_binary, write_binary, write_hdr
import multiprocessing as mp
import numpy as np
import ctypes
import cmath
import math
import sys
import os

M_PI = math.pi
xp, yp = None, None
nrow, ncol = None, None
eps = np.finfo(np.float64).eps  # "machine epsilon" for 64-bit floating point number
o2d1, o2d2, o2de, o3d1, o3d2, o3d3 = None, None, None, None, None, None
out_r, out_g, out_b, out_e1, out_e2, out_e3, out_opt, out_sopt, out_v1 = None, None, None, None, None, None, None, None, None
t11_p, t22_p, t33_p, t12_r_p, t12_i_p, t13_r_p, t13_i_p, t23_r_p, t23_i_p = None, None, None, None, None, None, None, None, None

class vec3:
    def __init__(self, a, b, c):
        self.a, self.b, self.c = a, b, c

    def norm(self):
        return math.sqrt(abs(self.a)*abs(self.a) + abs(self.b)*abs(self.b) + abs(self.c)*abs(self.c));

    def normalize(self):
        norm = self.norm()  # give vector l2 length of 1 
        self.a, self.b, self.c = self.a / norm, self.b / norm, self.c / norm

    def __truediv__(self, a):  # define division by scalar ( X / a ) operator
        return vec3(self.a / a, self.b / a, self.c / a)
    
    def __sub__(self, b):
        return vec3(self.a - b.b, self.b - b.b, self.c - b.c)
    
    def __str__(self):  # tostring()
        return ' '.join([str(x) for x in [self.a, self.b, self.c]])


def solve_cubic(a, b, c, d):
    _2t13, _2t23, sqrt3 = 2. ** 0.3333333333333333, 2. ** 0.6666666666666666, math.sqrt(3.)
    
    t1, t2 = b*(-2.*b*b + 9.*a*c) - 27.*a*a*d, 3.*a*c -b*b
    t0 = (t1 +  (4.*(t2*t2*t2) + (t1*t1)) ** 0.5) #  ,0.5))
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
        _A, _B = -1. + 0j, a + d + f
        _C = (-(a*d) - a*f - d*f + b*b.conjugate() + c*c.conjugate() + e*e.conjugate())
        _D = d*(a*f - c*c.conjugate()) + e*(b*c.conjugate() - a*e.conjugate()) + b.conjugate()*(-(b*f) + c*e.conjugate())
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
    return vec3(1. + 0j, # cf(1.,0.),
                -((A.a)*(A.e)-_lambda*(A.e)-(A.c)* (A.b).conjugate() )/((A.b)*(A.e)-(A.d)*(A.c)+_lambda*(A.c)),
                (-(A.b)* (A.b).conjugate() -_lambda*(A.a)+(A.d)*(A.a)-(A.d)*_lambda+(_lambda*_lambda))/((A.b)*(A.e)-(A.d)*(A.c)+_lambda*(A.c)))


def eig(A): #  L, E1, E2, E3): # herm3<cf> &A , vec3<cf> &L, vec3<cf> &E1, vec3<cf> &E2, vec3<cf> &E3){
    lambdas = A.solve_characteristic()
    e1, e2, e3 = eigv(A, lambdas.a), eigv(A, lambdas.b), eigv(A, lambdas.c)
    l1, l2, l3 = lambdas.a, lambdas.b, lambdas.c  

    e1.normalize()
    e2.normalize()
    e3.normalize()

    X = [[abs(l1), e1], [abs(l2), e2], [abs(l3), e3]]
    X.sort(reverse=True)  # sort eigenvectors by eigenvalue ( decreasing order )

    L = [X[i][0] for i in range(len(X))]
    L = vec3(L[0], L[1], L[2])
    E1, E2, E3 = X[0][1], X[1][1], X[2][1]

    d1, d2, d3 = (A*E1)/(L.a) - E1, (A*E2)/(L.b) - E2, (A*E3)/(L.c) - E3;
    diff = d1.norm() + d2.norm() + d3.norm() 
    return [L, E1, E2, E3]


def read_T3(d):
    global t11_p, t22_p, t33_p, t12_r_p, t12_i_p, t13_r_p, t13_i_p, t23_r_p, t23_i_p 
    sep = os.path.sep
    t11_p = read_binary(d + sep + 'T11.bin')[3]
    t22_p = read_binary(d + sep + 'T22.bin')[3]
    t33_p = read_binary(d + sep + 'T33.bin')[3]
    t12_r_p = read_binary(d + sep + 'T12_real.bin')[3]
    t12_i_p = read_binary(d + sep + 'T12_imag.bin')[3]
    t13_r_p = read_binary(d + sep + 'T13_real.bin')[3]
    t13_i_p = read_binary(d + sep + 'T13_imag.bin')[3]
    t23_r_p = read_binary(d + sep + 'T23_real.bin')[3]
    t23_i_p = read_binary(d + sep + 'T23_imag.bin')[3]

def lamcloude(a, b, c, z1, z2, z3):
    p = 1./3.  # cf tra, z1p, z2p, z3p, fac0, fac1, fac2, fac3, s1, s2, deta, tr3;
    tra = (a + b + c) / 3.
    z1p, z2p, z3p = z1.conjugate(), z2.conjugate(), z3.conjugate() # conjugate(z1), conjugate(z2), conjugate(z3)
    fac0 = z1 * z1p + z2 * z2p + z3 * z3p

    s1 = a * b + a * c + b * c - fac0
    deta = a * b * c - c * z1 * z1p - b * z2 * z2p + z1 * z2p * z3 + z1p * z2 * z3p - a * z3 * z3p
    s2 = a * a - a * b + b * b - a * c - b * c + c * c + 3. * fac0
    fac1 = 27. * deta - 27. * s1 * tra + 54. * (tra ** 3.) # pow(tra, 3.)
    tr3 = fac1 + cmath.sqrt( (fac1 ** 2.)- 4. * (s2 ** 3.))  #pow(s2, 3.))
    fac2 = 1. + (1j * math.sqrt(3.))
    fac3 = 1. - (1j * math.sqrt(3.))

    e1 = (tra + pow(tr3, p) / (3. * pow(2., p)) + (s2 * pow(2., p) + eps) / (3. * pow(tr3, p) + eps)).real
    e2 = (tra - (fac2 * s2) / (3. * pow(tr3, p) * pow(2., 2. * p) + eps) - (fac3 * pow(tr3,p)) / (6.*pow(2., p) + eps)).real
    e3 = (tra - (fac3 * s2) / (3. * pow(2., 2.*p) * pow(tr3, p) + eps) - (fac2 * pow(tr3,p)) / (6.*pow(2., p) + eps)).real
    
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


def rank1_t3(e1, v1, v2, v3):  #  generate T3 rank 1
    t11c = e1 * v1 * (v1.conjugate())
    t12c = e1 * v1 * (v2.conjugate())
    t13c = e1 * v1 * (v3.conjugate())
    t22c = e1 * v2 * (v2.conjugate())
    t23c = e1 * v2 * (v3.conjugate())
    t33c = e1 * v3 * (v3.conjugate())
    return [t11c, t12c, t13c, t22c, t23c, t33c]


def decom(i):   # calculate decom for pixel at linear index "i"
    if i % 10000 == 0:
            print(i, 100. * (i + 1) / (nrow * ncol))
    debug = (i == ncol * yp + xp)  # check if we're on target pixel

    if True:
        t11 = t11_p[i]
        t22 = t22_p[i]
        t33 = t33_p[i]
        t12_r = t12_r_p[i]
        t12_i = t12_i_p[i]
        t13_r = t13_r_p[i]
        t13_i = t13_i_p[i]
        t23_r = t23_r_p[i]
        t23_i = t23_i_p[i]

        if debug:
            print("t11", t11)
            print("t22", t22)
            print("t33", t33)
    
        a = t11 + 0j
        b = t22 + 0j
        c = t33 + 0j
        z1 = t12_r + t12_i * 1j
        z2 = t13_r + t13_i * 1j
        z3 = t23_r + t23_i * 1j
    
        # aliases
        t12c = z1
        t13c = z2
        t23c = z3
        t11c = a
        t22c = b
        t33c = c
    
        # /* avoid 0 elements.. conditioning */
        eps2 = (a + b + c) * (1.0e-9) + eps;
        F = ((float(nrow) + float(ncol)) / 2.);
        z1 = z1 + eps2 * F
        z2 = z2 + eps2 * F 
        z3 = z3 + eps2 * F
        a = a + eps2 * F
        b = b + eps2 * F
        c = c + eps2 * F

        [e1, e2, e3, v1, v2, v3] = lamcloude(a, b, c, z1, z2, z3) #  e1, e2, e3, v1, v2, v3);
        if debug:
            print("lamcloude")
            print("e1", e1)
            print("e2", e2)
            print("e3", e3)
            print("v1", v1)
            print("v2", v2)
            print("v3", v3)

        # // rank 1 t3
        [t11c, t12c, t13c, t22c, t23c, t33c] = rank1_t3(e1, v1, v2, v3)

        # // generate alpha etc. eigenvector parameters
        alpha = math.acos(abs(v1));
        phi = cmath.phase(t12c);
        theta = cmath.phase((t22c - t33c) + 2. * 1j * t23c.real) / 4.
    
        # // generate RGB colour composite from multiple eigenvector angles
        dn = alpha * 2. / M_PI # // alpha angle in red channel
        theta2 = theta + (theta > M_PI / 4.) * (M_PI / 4. - theta)
        theta2 = theta2 + (theta2 < -M_PI / 4.) * (-M_PI / 4. - theta2)
        vn = (theta2 + M_PI / 4.) * 2. / M_PI   # // az slope is green
        sn = abs(phi) / M_PI  # // mag of Pauli phase is blue (180 is Bragg)
    
        out_r = dn
        out_g = vn
        out_b = sn
        out_e1 = e1
        out_e2 = e2
        out_e3 = e3
    
        #  // project data onto null channels // null_vecs=[o2d o3d];
        z1 = o2d1.conjugate()*v1 + o2d2.conjugate()*v2 + o2d3.conjugate()*v3  # // oconj=o2d';
        z2 = o3d1.conjugate()*v1 + o3d2.conjugate()*v2 + o3d3.conjugate()*v3  #// oconj=o3d';
    
        #  find optimum weights
        popt = cmath.phase(z2 * z1.conjugate()) * 180. / M_PI
        za = (z1*z1.conjugate() - z2*z2.conjugate()) + 1j * 2.*abs(z1)*abs(z2)
        aopt = cmath.phase(za) * 90. / M_PI
        ar = aopt * M_PI / 180.
        br = popt * M_PI / 180.
    
        # // optimum weight vector
        w1 = (math.cos(ar) * o2d1 + math.sin(ar) * cmath.exp(1j * br) * o3d1).conjugate()
        w2 = (math.cos(ar) * o2d2 + math.sin(ar) * cmath.exp(1j * br) * o3d2).conjugate()
        w3 = (math.cos(ar) * o2d3 + math.sin(ar) * cmath.exp(1j * br) * o3d3).conjugate()
    
        # find optimum subspace signal
        zopt = w1 * v1 + w2 * v2 + w3 * v3
        ip = abs(zopt * zopt.conjugate()) 
        ip_eps = ip + eps;
        sopt = 10. * math.log(ip_eps) / math.log(10.); # optimum normalised power
    
        sp = t11c + t22c + t33c  # span power
        abs_sp = abs(sp)
        pwr = 10. * math.log(abs(sp)) / math.log(10.)  # span channel
    
        sm = math.fabs(t33)
        hv = 10. * math.log(sm) / math.log(10.)
        sm2 = sopt + pwr
    
        opt = pow(10., sm2 / 10.); #// linear opt channel
    
        # out_opt[i] =  opt;
        # out_v1[i] = sopt;
        return [i, opt, sopt, out_r, out_g, out_b, out_e1, out_e2, out_e3]
        '''
        //out_hv[i] = hv;
        //out_sm[i] = sm;
    
        //out_pwr[i] = (float)pwr;
        //out_sopt[i] = (float)sopt;
        //out_abs_sp[i] = (float)abs_sp;
    
        //out_ar[i] = (float) abs(ar);
        //out_br[i] = (float) abs(br);
        '''
    else:
        pass


def worker(task_queue, result_array, job_count, chunk_size):
    """ Worker function to process tasks (one integer at a time) and store results in the shared array """
    while True:
        job = task_queue.get()  # get next job
        if job is None:  # sentinel value to stop worker
            break 
        start_idx = job * chunk_size  # start index for results chunk (this job)
        result_array[start_idx: start_idx + chunk_size] = decom(job)  # put the results at the appropriate location

def work_queue(job_count, num_workers, chunk_size):
    task_queue = mp.Queue()  # task queue object
    result_array = mp.Array('d', [math.nan] * (job_count * chunk_size))  # init array with null

    processes = []
    for _ in range(num_workers):  # start the workers
        p = mp.Process(target=worker, args=(task_queue, result_array, job_count, chunk_size))
        p.start()
        processes.append(p)

    for job in range(job_count):  # add jobs to task queue
        task_queue.put(job)

    for _ in range(num_workers):  # add sentinel values to stop the workers
        task_queue.put(None)

    for p in processes:  # wait to finish
        p.join()

    return result_array  # return list(result_array)  # convert to regular list 


x = read_config('../T3/config.txt')
nrow, ncol = x['nrow'], x['ncol']
npx = nrow * ncol
read_T3('../T3/')  # load the T3 matrix data

# initialize output variables
out_r, out_g, out_b, out_e1, out_e2, out_e3, out_opt, out_sopt, out_v1 = [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)], [math.nan for i in range(npx)]

xp, yp = int(sys.argv[2]), int(sys.argv[1])  # cloude_decom T3 313 798 # col/ row for test target

if xp < 0 or xp > ncol:
    err("x coord out of bounds")

if yp < 0 or yp > nrow: 
    err("y coord out of bounds")

n_use = 1;
i = yp * ncol + xp;
t11 = t11_p[i]
t22 = t22_p[i]
t33 = t33_p[i]
t12_r = t12_r_p[i]
t12_i = t12_i_p[i]
t13_r = t13_r_p[i]
t13_i = t13_i_p[i]
t23_r = t23_r_p[i]
t23_i = t23_i_p[i]

print("target info:")
print("T11", t11)
print("T22", t22)
print("T33", t33)


if len(sys.argv) > 3:
    pass
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
z1 = z1 + eps2 * F   # %randn(sx,sy);
z2 = z2 + eps2 * F   # %randn(sx,sy);
z3 = z3 + eps2 * F   # %randn(sx,sy);
a = a + eps2 * F   # %randn(sx,sy);
b = b + eps2 * F   # %randn(sx,sy);
c = c + eps2 * F   # %randn(sx,sy);

T = herm3(a, z1, z2, b, z3, c)
# herm3<cf> T(a, z1, z2, b, z3, c);
# cout << "T" << endl << T << endl;

# vec3<cf> L, E1, E2, E3;
# eig(T, L, E1, E2, E3);
L, E1, E2, E3 = eig(T)
print("L", L)  # cout << "L" << endl << L << endl;
print("E2", E2)
print("E3", E3)

# / dont forget to test eigs !!!!!!!
o2d1 = E2.a # [0] # at(E2, 0);
o2d2 = E2.b # [1] # at(E2, 1);
o2d3 = E2.c # [2] #at(E2, 2);

o3d1 = E3.a # [0] # at(E3, 0);
o3d2 = E3.b # [1] # at(E3, 1);
o3d3 = E3.c # [2] #  at(E3, 2);

# cout << "o2: " << o2d1 << o2d2 << o2d3 << endl;
# cout << "o3: " << o3d1 << o3d2 << o3d3 << endl;


# later generalize to i:

'''for i in range(npx):
    #if i % 1000 == 0:
    #    print(i)
    decom(i)
'''
# results = parfor(decom, range(nrow*ncol))
# print("results", results)


job_count = nrow * ncol  # 10  # Total number of jobs (more jobs than workers)
num_workers = 32  # Number of worker processes (threads)
chunk_size = 9  # number of elements returned by decom() function
results = work_queue(job_count, num_workers, chunk_size)

for i in range(nrow * ncol):
    job_i = int(results[i* chunk_size])
    out_opt[job_i] =  results[i * chunk_size + 1]
    out_sopt[job_i] = results[i + chunk_size + 2]
    out_r[job_i] = results[i + chunk_size + 3]
    out_g[job_i] = results[i + chunk_size + 4]
    out_b[job_i] = results[i + chunk_size + 5]


'''
[i, opt, sopt, out_r, out_g, out_b, out_e1, out_e2, out_e3]
'''
print("+w opt.bin")
write_binary(out_opt, "opt.bin"); write_hdr("opt.hdr", ncol, nrow, 1, ["opt.bin"])
write_binary(out_sopt, "sopt.bin"); write_hdr("sopt.hdr", ncol, nrow, 1, ["sopt.bin"])
write_binary(out_r, "r.bin"); write_hdr("r.hdr", ncol, nrow, 1, ["r.bin"])
write_binary(out_g, "g.bin"); write_hdr("g.hdr", ncol, nrow, 1, ["g.bin"])
write_binary(out_b, "b.bin"); write_hdr("b.hdr", ncol, nrow, 1, ["b.bin"])



# test case: 
# python3 cloude_decom.py  313 798
