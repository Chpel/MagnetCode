#include <iostream>
#include <cstdlib>
#include"mcmc.h"
#include <string>

int main(int argc, char* argv[]) {
    //std::cout << "Hello, World!" << std::endl;
    long int N = std::stoi(argv[1]); //std::atoi(argv[1]);
    //std::cout << N << " " << std::endl;
    double J = 0.001 * (double)std::stoi(argv[2]); //std::atof(argv[2]);
    std::cout << J << " " << std::endl;
    std::string l = argv[3];
    //int nSim = std::stoi(argv[3]);
    //double  h = 0.1*(double)std::stoi(argv[3]);
    time_t start, end;
    time(&start);
    
    //int N = 10;
    //double J = 0.5;
    //std::string l = "cubic";

    Protein p(N, l);

    p.MC(J, 0);
    //Protein p(20);
    //p.MC(0.8);


    time(&end);

    double dif = difftime(end, start);
    printf("Time = %lf \n", dif);

    /*for (int i = 4; i < 15 ; i++ )
    {
        Protein p(i);
        p.MC(0.5);
    }*/

    //p.MC(0.8);
    return 0;
}