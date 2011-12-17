# size of the grid (in every direction) to map
SIZE:=4

all: sfc++ parallel_example

sfc++: sfc++.cpp sfc++.hpp Makefile
	g++ -DSIZE=$(SIZE) -I$$HOME/include -O3 -W -Wall -Wextra -pedantic sfc++.cpp -o sfc++

parallel_example: parallel_example.cpp sfc++.hpp Makefile
	g++ -fopenmp -DSIZE=$(SIZE) -I$$HOME/include -O3 -W -Wall -Wextra -pedantic parallel_example.cpp -o parallel_example

c: clean
clean:
	rm -f sfc++ parallel_example 2d.dat 3d.dat serial.dat parallel.dat

