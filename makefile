default:
	@echo "Please specify: c to compile example.c, cpp to compile example.cpp, so to compile libraries, clean to rm *.so *.o"

c : libraries example1.c
	gcc -o example1 example1.c -L./ -lcoolc -lcoolcpp
	@echo "Done. Don't forget to set your LD_LIBRARY_PATH when running the example!"
	
cpp : example2.cpp
	g++ -c example2.cpp Cool.cpp CoolManager.cpp Simplex.cpp PSI.cpp
	g++ -o example2 example2.o Cool.o CoolManager.o Simplex.o PSI.o
	
so : libraries

libraries : libcoolc.so libcoolcpp.so
	g++ -fpic -shared CInterface.cpp -L. -I. -lcoolcpp -o libcoolc.so
	g++ -fpic -shared Simplex.cpp Cool.cpp CoolManager.cpp -o libcoolcpp.so -I.

PSI :
	g++ -o psi PSImain.cpp PSI.cpp Simplex.cpp Cool.cpp CoolManager.cpp

clean :
	rm *.so *.o
