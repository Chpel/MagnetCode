//
// Created by kamilla on 13.10.2019.
//
//#ifndef MC_CPP_MCMC_H
//#define MC_CPP_MCMC_H
#include <list>
#include <valarray>
#include <algorithm>
#include <random>
#include <tuple>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <math.h>
#include <set>
#include "observable.h"
#include "lattice.h"


class Protein {
public:

    Protein();
    Protein(long int length, std::string l);

    void Reconnect(int j);
    void Reconnect1(int j);
    bool IsEndInStuck();

    void MC(  double J=0, double h=0, int nSimulation = 0, long int steps_to_equilibrium = 40000000, long int mc_steps = 5000000000000, bool radius = false);
    //void MC(double J = 0, double h = 0, int nSimulation = 0, long int steps_to_equilibrium = 1000, long int mc_steps = 100000, bool radius = false);

    void save_calcs();

    void calc_bulk();

    long int radius();
    void radius_gyration();
    void radius_gyration1();

    void count_contacts();
    bool CheckAndFlipNodes(long int& coord, int& sign);

    void coord_form();
    void write_file(long int i);
public:

    Lattice* lattice;

    std::valarray<long int> ordered_coords;

    std::valarray<short> directions; //их n-1;
    //если в directions[10] стоит 0, значит, двигаемся вправо из координаты 10

    std::valarray<char> sequence_on_lattice;
    std::valarray<coord_check_t> next_monomers;
    std::valarray<coord_check_t> previous_monomers;
    coord_t end_conformation = 0, start_conformation = 0;

    std::queue<long int>  spins_in_cluster;

    mc_stats::ScalarObservable<long int> dists;
    mc_stats::ScalarObservable<double> gyration;

    mc_stats::ScalarObservable<long double> energy; //сохранение энергии
    mc_stats::ScalarObservable<long double> energy_sq;
    mc_stats::ScalarObservable<long double> energy_4;

    mc_stats::ScalarObservable<long double> magnetization; //сохранение намагниченности
    mc_stats::ScalarObservable<long double> magnetization_sq;
    mc_stats::ScalarObservable<long double> magnetization_4;

    std::vector<mc_stats::ScalarObservable<double>> bulk;
    std::vector<int> bulk_now;

    mc_stats::ScalarObservable<double> eigs1;
    mc_stats::ScalarObservable<double> eigs2;
    mc_stats::ScalarObservable<double> aratio;

    long int number_of_monomers = 0;
    long int E = 0; // -1* число топологических контактов текущей конформации
    long int current_H_counts = 0;

    std::map <long int, long long int> count_E;
    std::map <long int, long long int> count_M;

    std::vector<long double> sum_coord;
    double P_add = 0;

private:
    double J = 0;
    double h = 0;
    std::string l;
    int nSimulation = 0;

};


//#endif //MC_CPP_MCMC_H