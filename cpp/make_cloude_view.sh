g++ 	-c -w  -O4  imv.cpp &
g++ 	-c -w  -O4  newzpr.cpp &
g++ 	-c -w  -O4  my_math.cpp &
g++ 	-c -w  -O4  image.cpp &
g++ 	-c -w  -O4  util.cpp &
wait
g++ 	imv.o newzpr.o my_math.o	image.o	util.o	-o	cloude_view  -lm -lpthread  -lGL -lGLU -lglut
