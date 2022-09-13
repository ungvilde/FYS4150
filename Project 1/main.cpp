#include <iostream> 
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <chrono>

#include "utils.hpp"
#include "algorithms.hpp"


int main()
{
    int N = 100; //number of steps along x-axis to compute exact values
    double h = 1./N; 
    std::vector<double> xvalues(N+1);
    std::vector<double> uvalues(N+1);
    xvalues[0] = 0.; //boundary conditions
    uvalues[0] = 0.;
    xvalues[N] = 1.;
    uvalues[N] = 0.;

    //then loop to assign values incremented by small step size to every element in vector
    for (int i = 1; i < N; i++)
    {
        xvalues[i] = xvalues[i-1] + h;
        uvalues[i] = u( xvalues[i] );
    }

    // save in file
    std::vector<std::vector<double>> datavalues;
    datavalues.push_back(xvalues);
    datavalues.push_back(uvalues);
    write_file("exact.txt", datavalues);

    // now we make approximations with varying number of steps along x-axis
    // using general algorithm
    int Nsteps = 10;
    while(Nsteps <= 10000000) //generated datasets for num steps up to 10**7
    {
    std::vector< std::vector<double> > approxvalues = general_algorithm(Nsteps);  
    std::string filename = "general_approx_N" +  std::to_string(Nsteps) + ".txt";
    write_file(filename, approxvalues);
    Nsteps=Nsteps*10;
    }
    
    // here we time the algorithms for varying number of steps
    Nsteps = 10; //initial num. steps
    std::ofstream ofile;
    ofile.open("timing.txt");
    while(Nsteps <= 1000000) //up to 10**6
    {
    double time_general = time_general_algorithm(Nsteps, 10);
    double time_special = time_special_algorithm(Nsteps, 10);

    ofile << log10(Nsteps) << " " << scientific_format(time_general, 10, 3) << " " << scientific_format(time_special, 10, 3) << std::endl;
    Nsteps=Nsteps*10;
    }
    ofile.close();
    return 0;
}

