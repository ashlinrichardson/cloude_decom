/* Poincare Orthogonality filter [1] based on Rank-1 Cloude decomposition [2] 
adapted from lamcloude.m and nullopt.m by Shane R Cloude

compile (terminal):
  g++ -O3 cloude.cpp -o cloude.exe -lm

run:
  ./cloude.exe T3 # for T3 matrix data in folder T3

output:
  OPT.bin

References
[1] Shane R Cloude, Ashlin Richardson, Generalized Poincare Orthogonality: A New Approach to POLSAR Data Analysis, proc APSAR 2021
[2] Shane R Cloude, Target Detection Using Rank-1 Polarimetric Processing, IEEE GRSL (2020)

C++ implementation 20210601 Ash Richardson, Senior Data Scientist, BC Wildfire Service */
#include<math.h>
#include<float.h>
#include<stdio.h>
#include<stdlib.h>
#include<memory.h>
#include<string.h>
#include"matrix2.h"
#include"matrix3.h"
#include<pthread.h>
#include<unistd.h>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  //define something for Windows (32-bit and 64-bit)
#elif __APPLE__
  // apple stuff
#else 
  // assume LINUX 
  #include<sys/sysinfo.h>
#endif

#define str string
#define eps DBL_EPSILON
#define for0(i,n) for(i = 0; i < n; i++) /* for loop shorthand */

#define N_IN 9 /* number of input files */
#define T11 0 /* T3 matrix input data indexing */
#define T12_re 1
#define T12_im 2
#define T13_re 3
#define T13_im 4
#define T22 5
#define T23_re 6
#define T23_im 7
#define T33 8

const char* T_fn[] = {
  "T11.bin",
  "T12_real.bin",
  "T12_imag.bin",
  "T13_real.bin",
  "T13_imag.bin",
  "T22.bin",
  "T23_real.bin",
  "T23_imag.bin",
  "T33.bin"};

float ** T; /* T3 matrix elements */
float ** T_f; /* filtered T3 matrix elements */
float ** out_d; /* output buffers */

#define N_OUT 10 /* number of output bands.. was 17 */
#define _e1_f 0 /* output data indexing */
#define _e2_f 1
#define _e3_f 2
#define _r_f 3
#define _g_f 4
#define _b_f 5
#define _opt_f 6
#define _v1_f 7
#define _v2_f 8
#define _v3_f 9
#define _hv_f 10
#define _sm_f 11
#define _pwr_f 12
#define _sopt_f 13
#define _abs_sp_f 14
#define _ar_f 15
#define _br_f 16

const char * out_fn[] = {
  "e1.bin",
  "e2.bin",
  "e3.bin",
  "r.bin",
  "g.bin",
  "b.bin",
  "OPT.bin",
  "v1.bin",
  "v2.bin",
  "v3.bin"};
/*,
  "hv.bin",
  "sm.bin",
  "pwr.bin",
  "sopt.bin",
  "abs_sp.bin",
  "ar.bin",
  "br.bin"};
*/

char sep(){
  #ifdef _WIN32
    return '\\'; /* windows path sep */
  #else
    return '/'; /* posix */
  #endif
}

#define MAX_ARRAYS 1024 /* more than how many arrays we need */
int n_arrays = 0; /* count arrays initialized */
void ** arrays; /* track mallocs, free at end */

void * alloc(size_t n){
  void * d = malloc(n); /* create array */
  if(!d) err("failed to allocate memory");
  memset(d, '\0', n); /* must touch memory on windows */
  arrays[n_arrays ++] = d; /* save pointer to free later */
  return d;
}

float * falloc(size_t n){
  return (float *)alloc(n * sizeof(float)); /* float32 array */
}

#define READ 1
#define WRITE 0

FILE * open(const char * fn, int mode){
  printf("+%s %s\n", mode?"r":"w", fn); /* open a file for read or write */
  FILE * f = fopen(fn, mode?"rb":"wb");
  if(!f){
     err("file access failed");
  }
  return f;
}

void read_config(char * file_name, int * nrow, int * ncol){
  /* based on PolSARPro by Eric POTTIER and Laurent FERRO-FAMIL */
  int i;
  size_t x;
  char tmp[4096];
  FILE * f = open(file_name, READ);
  for0(i, 2) x = fscanf(f, "%s\n", tmp);
  *nrow = atoi(tmp); // number of rows
  for0(i, 3) x = fscanf(f, "%s\n", tmp);
  *ncol = atoi(tmp); // number of cols
  fclose(f);
  printf("nrow %d ncol %d\n", *nrow, *ncol);
}

float * read(const char * file_name, size_t n_float){
  FILE * f = open(file_name, READ);
  float * d = falloc(n_float);
  size_t nr = fread(d, sizeof(float), n_float, f);
  if(nr != n_float){
    printf("Expected number of floats: %zu\n", n_float);
    printf("Number of floats read: %zu\n", nr);
    err("Unexpected float read count");
  }
  fclose(f);
  return d; /* return array of floats we read in */
}

#define STR_MAX 4096 /* string variable length */
void hwrite(char * bfn, size_t nrow, size_t ncol, size_t nband){
  size_t i;
  char hfn[STR_MAX];
  size_t L = strlen(bfn);
  strcpy(hfn, bfn); /* change ext to hdr*/
  hfn[L - 3] = 'h'; hfn[L - 2] = 'd'; hfn[L - 1] = 'r';

  FILE * f = open(hfn, WRITE);
  fprintf(f, "ENVI\n");
  fprintf(f, "samples = %zu\n", ncol);
  fprintf(f, "lines = %zu\n", nrow);
  fprintf(f,"bands = %zu\n", nband);
  fprintf(f, "header offset = 0\n");
  fprintf(f, "file type = ENVI Standard\n");
  fprintf(f, "data type = 4\n");
  fprintf(f, "interleave = bsq\n");
  fprintf(f, "byte order = 0\n");
  fprintf(f, "band names = {band 1");
  for0(i, nband - 1) fprintf(f, ",\nband %zu", i + 2);
  fprintf(f, "}\n");
  fclose(f);
}

void swp(double * a, double * b){
  double tmp = *a;
  *a = *b;
  *b = tmp;
}

#define mtx_lock pthread_mutex_lock
#define mtx_unlock pthread_mutex_unlock
pthread_mutex_t print_mtx;
pthread_attr_t pthread_attr; // specify threads joinable
pthread_mutex_t pthread_next_j_mtx; // work queue
size_t pthread_next_j; // next job to run
size_t pthread_start_j, pthread_end_j; // start and end indices for job
void (*pthread_eval)(size_t); // function pointer to execute in parallel, over range start_j:end_j inclusive

void init_mtx(){
  // mutex setup
  pthread_mutex_init(&print_mtx, NULL);
  pthread_mutex_init(&pthread_next_j_mtx, NULL);
}

void cprint(string s){
  pthread_mutex_lock(&print_mtx);
  cout << s << endl;
  pthread_mutex_unlock(&print_mtx);
}

void * worker_fun(void * arg){
  size_t k, my_next_j;
  k = (size_t)arg;
  // cprint(str("worker_fun(") + std::to_string(k) + str(")"));

  while(1){
    // try to pick up a job
    mtx_lock(&pthread_next_j_mtx);
      my_next_j = pthread_next_j ++; // index of data this thread should pick up if it can
    mtx_unlock(&pthread_next_j_mtx);

    if(my_next_j >= pthread_end_j)
      return(NULL);
    
    pthread_eval(my_next_j); // perform action segment
  }
}

void parfor(size_t start_j, size_t end_j, void(*eval)(size_t)){
  // ideally the worker fun would be an inline (inside of here)
  pthread_eval = eval; // set global function pointer
  pthread_start_j = start_j;
  pthread_end_j = end_j;

  pthread_next_j = start_j; // pthread_next_j_mtx is the lock on this variable
  size_t n_cores = sysconf(_SC_NPROCESSORS_ONLN);
  cout << "Number of cores: " << n_cores << endl;

  // allocate threads, make threads joinable
  pthread_attr_init(&pthread_attr);
  pthread_attr_setdetachstate(&pthread_attr, PTHREAD_CREATE_JOINABLE);
  pthread_t * my_pthread = new pthread_t[n_cores];
  size_t j;
  for0(j, n_cores) pthread_create(&my_pthread[j], &pthread_attr, worker_fun, (void *)j);

  // wait for threads to finish
  for0(j, n_cores) pthread_join(my_pthread[j], NULL);

  // delete my_pthread;
  // cprint(str("return parfor()"));
}

void lamcloude(cf a, cf b, cf c, cf z1, cf z2, cf z3, double & e1, double & e2, double & e3, cf & v1, cf & v2, cf & v3){
  double p = 1./3.;
  cf tra, z1p, z2p, z3p, fac0, fac1, fac2, fac3, s1, s2, deta, tr3;

  tra = (a + b + c) / 3.;
  z1p = conj(z1); z2p = conj(z2); z3p = conj(z3);

  fac0 = z1 * z1p + z2 * z2p + z3 * z3p;

  s1 = a * b + a * c + b * c - fac0;
  deta = a * b * c - c * z1 * z1p - b * z2 * z2p + z1 * z2p * z3 + z1p * z2 * z3p - a * z3 * z3p;
  s2 = a * a - a * b + b * b - a * c - b * c + c * c + 3. * fac0;
  fac1 = 27. * deta - 27. * s1 * tra + 54. * pow(tra, 3.);
  tr3 = fac1 + sqrt(pow(fac1, 2.) - 4. * pow(s2, 3.));
  fac2 = 1. + J*sqrt(3.);
  fac3 = 1. - J*sqrt(3.);

  e1 = real(tra + pow(tr3, p) / (3. * pow(2., p)) + (s2 * pow(2., p) + eps) / (3. * pow(tr3, p) + eps));
  e2 = real(tra - (fac2 * s2) / (3. * pow(tr3, p) * pow(2., 2. * p) + eps) - (fac3 * pow(tr3,p)) / (6.*pow(2., p) + eps));
  e3 = real(tra - (fac3 * s2) / (3. * pow(2., 2.*p) * pow(tr3, p) + eps) - (fac2 * pow(tr3,p)) / (6.*pow(2., p) + eps));

  if(e1 < e3) swp(&e1, &e3); // sort eigenvalues
  if(e1 < e2) swp(&e1, &e2);
  if(e2 < e3) swp(&e2, &e3);
  if(e1 < e3 || e1 < e2 || e2 < e3)
  printf("Warning: not sorted (%e, %e, %e)\n", e1, e2, e3);

  v2 = ((a - e1) * z3 - z1p * z2) / ((b - e1) * z2 - z3 * z1 + eps); // dominant eigenvector
  v3 = (e1 - a - z1 * v2) / (z2 + eps);
  v1 = 1.; // ones(size(v2));

  double av1 = abs(v1);
  double av2 = abs(v2);
  double av3 = abs(v3);
  double n = sqrt(av1 * av1 + av2 * av2 + av3 * av3) + eps;

  // normalised components as output
  v1 /= n; v2 /= n; v3 /= n;
}

void rank1_t3(double e1, cf v1, cf v2, cf v3, cf & t11c, cf & t12c, cf & t13c, cf & t22c, cf & t23c, cf & t33c){
  // generate T3 rank 1
  t11c = e1 * v1 * conj(v1); t12c = e1 * v1 * conj(v2); t13c = e1 * v1 * conj(v3);
  t22c = e1 * v2 * conj(v2); t23c = e1 * v2 * conj(v3);
  t33c = e1 * v3 * conj(v3);
}

int NROW, NCOL;

float * out_e1, * out_e2, * out_e3;
float * out_r, * out_g, * out_b;
float * out_opt;
float * out_v1, * out_v2, * out_v3; //, * out_hv, * out_sm, * out_pwr, * out_sopt, * out_abs_sp, * out_ar, * out_br;

float * t11_p;
float * t12_r_p;
float * t12_i_p;
float * t13_r_p;
float * t13_i_p;
float * t22_p;
float * t23_r_p;
float * t23_i_p;
float * t33_p;

cf o2d1;//( 0.069259936421070 + 0.000570219370078 * J);
cf o2d2;//( 0.996125064964563 + 0.000000000000000 * J);
cf o2d3;//( 0.054196590398518 + 0.000566211391969 * J);
cf o3d1;//( 0.005343171277222 - 0.000578560652077 * J);
cf o3d2;//(-0.054696348535595 + 0.000610840317812 * J);
cf o3d3;//( 0.998488383567507 + 0.000000000000000 * J);

void decom(size_t i){
  // intermediary variables
  double t11, t12_r, t12_i, t13_r, t13_i, t22, t23_r, t23_i, t33;
  double e1, e2, e3, p;
  cf a, b, c, z1, z2, z3;
  t11 = (double)t11_p[i]; t22 = (double)t22_p[i]; t33 = (double)t33_p[i];
  t12_r = (double)t12_r_p[i]; t12_i = (double)t12_i_p[i];
  t13_r = (double)t13_r_p[i]; t13_i = (double)t13_i_p[i];
  t23_r = (double)t23_r_p[i]; t23_i = (double)t23_i_p[i];

  a = t11; b = t22; c = t33;
  z1 = t12_r + t12_i * J;
  z2 = t13_r + t13_i * J;
  z3 = t23_r + t23_i * J;

  /* avoid 0 elements.. conditioning */
  cf eps2 = (a + b + c) * (1.0e-9) + eps;
  cf F(((double)NROW + (double)NCOL) / 2.);
  z1 = z1 + eps2 * F; // %randn(sx,sy);
  z2 = z2 + eps2 * F; // %randn(sx,sy);
  z3 = z3 + eps2 * F; // %randn(sx,sy);
  a = a + eps2 * F; // %randn(sx,sy);
  b = b + eps2 * F; // %randn(sx,sy);
  c = c + eps2 * F; // %randn(sx,sy);

  //run lamcloude
  cf v1, v2, v3;
  lamcloude(a, b, c, z1, z2, z3, e1, e2, e3, v1, v2, v3);

  // rank 1 t3
  cf t11c, t12c, t13c, t22c, t23c, t33c;
  rank1_t3(e1, v1, v2, v3, t11c, t12c, t13c, t22c, t23c, t33c);

  // generate alpha etc. eigenvector parameters
  double alpha = acos(abs(v1));
  double phi = angle(t12c);
  double theta = angle((t22c - t33c) + 2. * J * real(t23c)) / 4.;

  // generate RGB colour composite from multiple eigenvector angles
  double dn = alpha * 2. / M_PI; // alpha angle in red channel
  double theta2 = theta + (theta > M_PI / 4.) * (M_PI / 4. - theta);
  theta2 = theta2 + (theta2 < -M_PI / 4.) * (-M_PI / 4. - theta2);
  double vn = (theta2 + M_PI / 4.) * 2. / M_PI; // az slope is green
  double sn = abs(phi) / M_PI; // mag of Pauli phase is blue (180 is Bragg)

  out_r[i] = (float) dn;
  out_g[i] = (float) vn;
  out_b[i] = (float) sn;

  out_e1[i] = (float) e1;
  out_e2[i] = (float) e2;
  out_e3[i] = (float) e3;

  // project data onto null channels // null_vecs=[o2d o3d];
  z1 = conj(o2d1)*v1 + conj(o2d2)*v2 + conj(o2d3)*v3; // oconj=o2d';
  z2 = conj(o3d1)*v1 + conj(o3d2)*v2 + conj(o3d3)*v3; // oconj=o3d';

  // find optimum weights
  cf popt = angle(z2 * conj(z1)) * 180. / M_PI;
  cf za = (z1*conj(z1) - z2*conj(z2)) + J*2.*abs(z1)*abs(z2);
  cf aopt = angle(za) * 90. / M_PI;
  cf ar = aopt * M_PI / 180.;
  cf br = popt * M_PI / 180.;

  // optimum weight vector
  cf w1 = cos(ar)*o2d1 + sin(ar)*exp(J*br)*o3d1;
  w1 = conj(w1);
  cf w2 = cos(ar)*o2d2 + sin(ar)*exp(J*br)*o3d2;
  w2 = conj(w2);
  cf w3 = cos(ar)*o2d3 + sin(ar)*exp(J*br)*o3d3;
  w3 = conj(w3);

  // find optimum subspace signal
  cf zopt = w1 * v1 + w2 * v2 + w3 * v3;
  double ip = abs(zopt * conj(zopt));
  double ip_eps = ip + eps;
  double sopt = 10. * log(ip_eps) / log(10.); // optimum normalised power

  cf sp = t11c + t22c + t33c; //span power
  double abs_sp = abs(sp);
  double pwr = 10. * log(abs(sp)) / log(10.); //span channel

  double sm = fabs(t33);
  double hv = 10. * log10(sm);
  double sm2 = sopt + pwr;

  double opt = pow(10., sm2 / 10.); // linear opt channel

  out_opt[i] = (float) opt;
  out_v1[i] = (float) sopt;

  //out_hv[i] = (float)hv;
  //out_sm[i] = (float)sm;

  //out_pwr[i] = (float)pwr;
  //out_sopt[i] = (float)sopt;
  //out_abs_sp[i] = (float)abs_sp;

  //out_ar[i] = (float) abs(ar);
  //out_br[i] = (float) abs(br);
}

int main(int argc, char ** argv){
  printf("cloude_decom.cpp\n");
  if(argc < 2) err("cloude.exe [input T3 directory] # [null target row] [null target col] # null target window size");
  char * path = argv[1]; /* T3 matrix data path */
  int xp = -1; int yp = -1;
  int xp2 = -1; int yp2 = -1;
  int ws = -1; int dw = -1;

  if(argc > 3){
    xp = atoi(argv[3]); //
    yp = atoi(argv[2]);
  }

  printf("target pixel (x,y): %d, %d\n", xp, yp);
  if(argc > 4){
    ws = atoi(argv[4]);
    printf("window size: %d\n", ws);
    if(ws%2 != 1){
      printf("Error: window size not odd\n"); exit(1);
    }
    dw = (ws - 1) / 2;
  }
  printf("**** ws %d dw %d\n", ws, dw);

  int i, j, k, np, nrow, ncol, di, dj, ii, jj, x, ix, jx, nw;

  char fn[STR_MAX];
  strcpy(fn, path); // path cat operator should be fxn
  fn[strlen(path)] = sep();
  strcpy(fn + strlen(path) + 1, "config.txt"); /* path to config.txt */
  read_config(fn, &nrow, &ncol); /* read image dimensions */
  np = nrow * ncol; /* number of px */
  NROW = nrow;
  NCOL = ncol;

  arrays = (void **) malloc(sizeof(void *) * MAX_ARRAYS); /* ptrs to free later */
  memset(arrays, '\0', sizeof(void *) * MAX_ARRAYS); /* touch mem on winOS */
  n_arrays = 0; /* start from beginning */

  T = (float **) alloc(sizeof(float *) * N_IN); /* input file bufs */
  for0(k, N_IN){
    strcpy(fn, path);
    fn[strlen(path)] = sep();
    strcpy(fn + strlen(path) + 1, T_fn[k]); /* [path][sep][fn] eg T3/T11.bin */
    T[k] = read(fn, np); /* read each input data band */
  }

  out_d = (float **) alloc(sizeof(float *) * N_OUT); /* output bands buffers */
  for0(i, N_OUT) out_d[i] = falloc(np); /* allocate output space */
  //for0(k, N_IN) filter(T[k], T_f[k], nrow, ncol, wsi); /* filter */

  t11_p = T[T11];
  t12_r_p = T[T12_re];
  t12_i_p = T[T12_im];
  t13_r_p = T[T13_re];
  t13_i_p = T[T13_im];
  t22_p = T[T22];
  t23_r_p = T[T23_re];
  t23_i_p = T[T23_im];
  t33_p = T[T33];

  out_e1 = out_d[_e1_f]; out_e2 = out_d[_e2_f]; out_e3 = out_d[_e3_f];
  out_r = out_d[_r_f]; out_g = out_d[_g_f]; out_b = out_d[_b_f];
  out_opt = out_d[_opt_f];
  out_v1 = out_d[_v1_f]; out_v2 = out_d[_v2_f]; out_v3 = out_d[_v3_f];
  
  /*
  out_hv = out_d[_hv_f];
  out_sm = out_d[_sm_f];
  out_pwr = out_d[_pwr_f];
  out_sopt = out_d[_sopt_f];
  out_abs_sp = out_d[_abs_sp_f];
  out_ar = out_d[_ar_f];
  out_br = out_d[_br_f];
*/

  /*
  %select reference coherency matrix
  Tc=[t11c(yp,xp) t12c(yp,xp) t13c(yp,xp);
  conj(t12c(yp,xp)) t22c(yp,xp) t23c(yp,xp);
  conj(t13c(yp,xp)) conj(t23c(yp,xp)) t33c(yp,xp)];

  Tc=Tc/trace(Tc);%rank-1 T matrix for pixel
  [v,d]=eig(Tc,'vector');
  [md,ii]=sort(d);
  o2d=v(:,ii(1));o3d=v(:,ii(2));%null space vectors

  %check eigenvector selection is correct
  check=max([o2d'*Tc*o2d o3d'*Tc*o3d]);%should be zero
  %check2=abs([o2d'*o2d o3d'*o3d o2d'*o3d]); %should be [1 1 0]
  if check>0.01
  disp('Warning: eigenvector error');
  check
  end;
  disp(['ref pixel location = ' num2str([yp xp])]);
  null_vecs=[o2d o3d]` */

  /* null vectors to project onto.. default values */
  o2d1 = cf( 0.069259936421070 + 0.000570219370078 * J);
  o2d2 = cf( 0.996125064964563 + 0.000000000000000 * J);
  o2d3 = cf( 0.054196590398518 + 0.000566211391969 * J);
  o3d1 = cf( 0.005343171277222 - 0.000578560652077 * J);
  o3d2 = cf(-0.054696348535595 + 0.000610840317812 * J);
  o3d3 = cf( 0.998488383567507 + 0.000000000000000 * J);
  double t11, t22, t33, t12_r, t12_i, t13_r, t13_i, t23_r, t23_i;
  cf a, b, c, z1, z2, z3;
  float n_use = 1;
  if(argc > 3){
    i = yp * ncol + xp;
    t11 = (double)t11_p[i]; t22 = (double)t22_p[i];
    t33 = (double)t33_p[i];
    t12_r = (double)t12_r_p[i]; t12_i = (double)t12_i_p[i];
    t13_r = (double)t13_r_p[i]; t13_i = (double)t13_i_p[i];
    t23_r = (double)t23_r_p[i]; t23_i = (double)t23_i_p[i];

    if(ws > 1){
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
    printf("N_USE: %f\n", n_use);
    t11 /= n_use; t22 /= n_use; t33 /= n_use;
    t12_r /= n_use; t12_i /= n_use;
    t13_r /= n_use; t13_i /= n_use;
    t23_r /= n_use; t23_i /= n_use;

    a = t11; b = t22; c = t33;
    z1 = t12_r + t12_i * J;
    z2 = t13_r + t13_i * J;
    z3 = t23_r + t23_i * J;

    /* avoid 0 elements.. conditioning */
    cf eps2 = (a + b + c) * (1.0e-9) + eps;
    cf F(((double)nrow + (double)ncol) / 2.);
    z1 = z1 + eps2 * F; // %randn(sx,sy);
    z2 = z2 + eps2 * F; // %randn(sx,sy);
    z3 = z3 + eps2 * F; // %randn(sx,sy);
    a = a + eps2 * F; // %randn(sx,sy);
    b = b + eps2 * F; // %randn(sx,sy);
    c = c + eps2 * F; // %randn(sx,sy);

    herm3<cf> T(a, z1, z2, b, z3, c);
    cout << "T" << endl << T << endl;

    vec3<cf> L, E1, E2, E3;
    eig(T, L, E1, E2, E3);
    cout << "L" << endl << L << endl;

    // dont forget to test eigs !!!!!!!
    o2d1 = at(E2, 0);
    o2d2 = at(E2, 1);
    o2d3 = at(E2, 2);

    o3d1 = at(E3, 0);
    o3d2 = at(E3, 1);
    o3d3 = at(E3, 2);

    cout << "o2: " << o2d1 << o2d2 << o2d3 << endl;
    cout << "o3: " << o3d1 << o3d2 << o3d3 << endl;
  }

  // run the decom
  parfor(0, np, decom);

  FILE * out_f[N_OUT];
  for0(i, N_OUT){
    strcpy(fn, path);
    fn[strlen(path)] = sep();
    strcpy(fn + strlen(path) + 1, out_fn[i]);
    hwrite(fn, nrow, ncol, 1); /* write envi header */
    out_f[i] = open(fn, WRITE); /* open output file */
    nw = fwrite(out_d[i], sizeof(float), np, out_f[i]);
    if(nw != np) err("failed to write expected number of floats");
    fclose(out_f[i]);
  }

  for0(k, n_arrays)
  free(arrays[n_arrays]); /* free anything we malloc'ed */

  free(arrays);
  return 0;
}
