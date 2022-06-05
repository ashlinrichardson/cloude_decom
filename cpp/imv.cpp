/* simple one-to one image viewer */
#include"util.h"
//#include"image.h"
#include<fstream>
#include"newzpr.h"
#include"image.h"
#include<iostream>
#include<stdexcept>

extern int NWIN;
extern size_t IMG_NR;
extern size_t IMG_NC;
extern size_t IMG_NB;
extern string IMG_FN;
extern string IMG_HFN;
extern SA<float> * SCENE_DAT;
extern myImg * SCENE_MYIMG;
extern void * SCENE_GLIMG;

void read_config(char * file_name, size_t * nrow, size_t * ncol){
  /* based on PolSARPro by Eric POTTIER and Laurent FERRO-FAMIL */
  int i;
  size_t x;
  char tmp[4096];
  FILE * f = fopen(file_name, "rb");
  for0(i, 2) x = fscanf(f, "%s\n", tmp);
  *nrow = atoi(tmp); // number of rows
  for0(i, 3) x = fscanf(f, "%s\n", tmp);
  *ncol = atoi(tmp); // number of cols
  fclose(f);
  printf("nrow %d ncol %d\n", *nrow, *ncol);
}

int main(int argc, char ** argv){
  size_t nr, nc, nb, np;
      
  IMG_FN = string("stack.bin"); // default image filename to load
  if(argc < 2) printf("cloude_view [infile] # [window size] # e.g.:\ncloude_view\ncloude_view stack.bin\ncloude_view stack.bin 3 # for window size 3");
  else IMG_FN = string(argv[1]);

  if(IMG_FN == string("stack.bin") && !exists(IMG_FN)){
    if(!exists("T11.bin") || !exists("T22.bin") || !exists("T33.bin")){
      err("T11.bin, T22.bin, T33.bin required for pauli visualization. Check a T3 dataset is in the present folder");
    }
    if(exists("T11.hdr")){
      parseHeaderFile("T11.hdr", nr, nc, nb);
      if(!exists("config.txt")){
        write_config(str("config.txt"), nr, nc);
      }
    }
    else if(exists("T11.bin.hdr")){
      parseHeaderFile("T11.bin.hdr", nr, nc, nb);
      if(!exists("config.txt")){
        write_config(str("config.txt"), nr, nc);
      }
    }
    else if(exists("config.txt")){
      read_config("config.txt", &nr, &nc);
    }
    else{
      err("expected header file at T11.hdr or T11.bin.hdr");
    }
    system("cat T22.bin T33.bin T11.bin > stack.bin");
    writeHeader("stack.hdr", nr, nc, 3);
  }
  else{
    if(exists("T11.hdr")){
      parseHeaderFile("T11.hdr", nr, nc, nb);
      if(!exists("config.txt")){
        write_config(str("config.txt"), nr, nc);
      }
    }
    else if(exists("T11.bin.hdr")){
      parseHeaderFile("T11.bin.hdr", nr, nc, nb);
      if(!exists("config.txt")){
        write_config(str("config.txt"), nr, nc);
      }
    }

  }

  if(!exists(IMG_FN)) err("failed to open input file"); // check if input file exists

  NWIN = 3;
  if(argc > 2) NWIN = atoi(argv[2]); // analysis window size
  printf("window size: %d\n", NWIN);
  string hfn(getHeaderFileName(IMG_FN)); // this section: get image scale
  parseHeaderFile(hfn, nr, nc, nb);
  np = nr * nc;
  IMG_HFN = hfn;
  IMG_NR = nr;
  IMG_NC = nc;
  IMG_NB = nb;

  zprManager * myManager = new zprManager(argc, argv); // window manager class
myZprManager = myManager;

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
  SCENE_GLIMG = (void*) myImage;
  myZpr->setTitle(string("Scene "));

  myImage->rebuffer(true, false, true);

  printf("mainloop\n");
  glutMainLoop();
  return 0;
}
