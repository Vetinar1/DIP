# makefile for testing c++ integration into C

all : main.c libcool.so
	gcc -o main main.c -L./ -lcool

libcool.so : CoolCInterface.cpp
	g++ -c -fPIC CoolCInterface.cpp -I./
	ar cvr libcool.a CoolCInterface.o

