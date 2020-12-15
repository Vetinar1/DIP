# makefile for testing c++ integration into C

all : main.c libcoolc.so libcoolcpp.so
	gcc -o main main.c -L./ -lcoolc -lcoolcpp

libcoolc.so : libcoolcpp.so CoolCInterface.cpp
	g++ -fpic -shared CoolCInterface.cpp -L. -I. -lcoolcpp -o libcoolc.so

libcoolcpp.so : CoolSimplex.cpp CoolCool.cpp
	g++ -fpic -shared CoolSimplex.cpp CoolCool.cpp -o libcoolcpp.so -I.

#libcoolcpp.so : CoolTest.cpp
#	g++ -fpic -shared CoolTest.cpp -o libcoolcpp.so -I.

clean :
	rm *.so