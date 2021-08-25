/* simple one-to one image viewer */
#include"util.h"
#include"image.h"
#include<fstream>
#include"newzpr.h"
#include<iostream>
#include<stdexcept>

extern int NWIN;
extern size_t IMG_NR;
extern size_t IMG_NC;
extern size_t IMG_NB;
extern string IMG_FN;
extern string IMG_HFN;

int main(int argc, char ** argv){
  IMG_FN = string("stack.bin"); // default image filename to load
  if(argc < 2) printf("imv.cpp: [infile] # [window size]\n");
  else IMG_FN = string(argv[1]);
  if(!exists(IMG_FN)) err("failed to open input file"); // check if input file exists

  NWIN = 1;
  if(argc > 2) NWIN = atoi(argv[2]); // analysis window size

  string hfn(getHeaderFileName(IMG_FN)); // this section: get image scale
  size_t nr, nc, nb, np;
  parseHeaderFile(hfn, nr, nc, nb);
  np = nr * nc;
  IMG_HFN = hfn;
  IMG_NR = nr;
  IMG_NC = nc;
  IMG_NB = nb;

  zprManager * myManager = zprManager::Instance(argc, argv); // window manager class
  SA<float> dat(np * nb); // the image data..
  SA<float> bb(np); // whole image, one-band buffer
  SCENE_DAT = &dat;

  FILE * f = fopen(IMG_FN.c_str(), "rb");
  if(!f) err("failed to open input file\n");
  fread(&dat[0], sizeof(float), np * nb, f);

  myImg a;
  a.initFrom(SCENE_DAT, nr, nc, nb);

  SCENE_MYIMG = &a;
  zprInstance * myZpr = myManager->newZprInstance(nr, nc, nb); // "Scene" / overview image WINDOW
  glImage * myImage = new glImage(myZpr, &a);
  SCENE_GLIMG = (void*)(glImage*)myImage;
  myZpr->setTitle(string("Scene "));

  myImage->rebuffer(true, false, true);

  glutMainLoop();
  return 0;
}
