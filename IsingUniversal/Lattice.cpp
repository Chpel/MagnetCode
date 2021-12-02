#include "lattice.h"

long long int int_pow(int lattice_side, int d) {
    long long int a = lattice_side;
    for (int i = 1; i < d; i++) {
        a *= lattice_side;
    }
    return a;
}