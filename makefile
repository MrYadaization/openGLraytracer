all: rayT

rayT: main.o primitives.o material.o light.o colour.o ray.o
	g++ main.o primitives.o material.o light.o colour.o ray.o -lGLU -lglut -lGL -lGLEW -o rayT

main.o: raytracer.cpp
	g++ -c raytracer.cpp -o main.o

primitives.o: primitives.cpp primitives.h
	g++ -c primitives.cpp

material.o: material.cpp material.h
	g++ -c material.cpp

light.o: light.cpp light.h
	g++ -c light.cpp

colour.o: colour.cpp colour.h
	g++ -c colour.cpp

ray.o: ray.cpp ray.h
	g++ -c ray.cpp