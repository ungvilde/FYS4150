# Project 3

build test with Mac Silicone: 
> CPATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib g++ test_FE_RK4_single_particle.cpp -std=c++11 -larmadillo src/*.cpp -I include -o test.exe

build test:
> g++ test_FE_RK4_single_particle.cpp -std=c++11 -larmadillo src/*.cpp -I include -o test.exe

run:
> ./test.exe

