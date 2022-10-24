compile: 
	g++ -O2  simulate_one_two_particles.cpp -std=c++11 -larmadillo src/* -I include -o simulate.out
	g++ -O2  simulate_time_dependence.cpp -std=c++11 -larmadillo src/* -I include -o time_dependent.out

run: 
	./time_dependent.out
	./simulate.out

all: compile run