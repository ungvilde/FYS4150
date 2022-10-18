# Project 3

build test with Mac Silicone: 
> CPATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib g++ test_FE_RK4_single_particle.cpp -std=c++11 -larmadillo src/*.cpp -I include -o test.exe

build test:
> g++ test_FE_RK4_single_particle.cpp -std=c++11 -larmadillo src/*.cpp -I include -o test.exe

run:
> ./test.exe

Idea:
- Make functions for all the different simulations/experiments, main.cpp is where we can run each experiment
- Save relevant values from each experiment
- Plot in python, where we define functions for different plot types, i.e plot_phase_space(position, velocity)

CPATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib g++ simulate_two_particles.cpp -std=c++11 -larmadillo src/*.cpp -I include -o simulate.exe