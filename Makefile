# size of the grid (in every direction) to map
SIZE:=4

all: sfc++

sfc++: sfc++.cpp sfc++.hpp Makefile
	g++ -DSIZE=$(SIZE) -I$$HOME/include -O3 -W -Wall -Wextra -pedantic sfc++.cpp -o sfc++

c: clean
clean:
	rm -f sfc++ 2d.dat 3d.dat

