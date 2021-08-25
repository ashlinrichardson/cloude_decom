g++ 	-c -w  -O3  image.cpp
g++ 	-c -w  -O3  my_math.cpp
g++ 	-c -w  -O3  util.cpp
g++ 	-c -w  -O3  newzpr.cpp
g++ 	-c -w  -O3  imv.cpp

if [[ "$OSTYPE" == "darwin"* ]]; then
    g++ 	imv.o newzpr.o my_math.o	image.o	util.o	-o	cloude_view  -lm -lpthread  -framework OpenGL -framework GLUT
else
    g++   imv.o newzpr.o my_math.o  image.o util.o  -o  cloude_view  -lm -lpthread  -lGL -lGLU -lglut
fi
