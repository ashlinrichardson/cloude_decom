echo "image"
g++ 	-c -w  -O3  image.cpp
echo "mymath"
g++ 	-c -w  -O3  my_math.cpp
echo "util"
g++ 	-c -w  -O3  util.cpp
echo "newzpr"
g++ 	-c -w  -O3  newzpr.cpp
echo "imv"
g++ 	-c -w  -O3  imv.cpp
echo "bin"
g++ 	imv.o newzpr.o my_math.o	image.o	util.o	-o	cloude_view  -lm -lpthread  -lGL -lGLU -lglut
