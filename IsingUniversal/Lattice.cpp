#include "lattice.h"

unsigned long long int int_pow(int lattice_side, int d) {
    unsigned long long int a = lattice_side;
    for (int i = 1; i < d; i++) {
        a *= lattice_side;
    }
    return a;
}