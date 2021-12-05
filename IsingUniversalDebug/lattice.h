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

unsigned long long int int_pow(int lattice_side, int d);

typedef unsigned long int coord_t;
typedef long long int coord_check_t;

class Lattice {
public:
	int lattice_side; //длина стороны решётки
	std::valarray<coord_t> map_of_contacts_int;
	std::vector<std::valarray<short>> coords;
    virtual int ndim() = 0; //число соседей узла
    virtual int d() = 0; //размерность пространства
	virtual void create_lattice(long int max_seq_size = 100) = 0;
};

class Square_Lattice_2D : public Lattice {
    int d() override { return 2; }
    int ndim() override { return 4; }
	void create_lattice(long int max_seq_size = 100) override {
        lattice_side = max_seq_size;
        map_of_contacts_int.resize(int_pow(lattice_side, d()) * ndim());
        coords.resize(d());
        for (int i = 0; i < d(); i++) {
            coords[i].resize(int_pow(lattice_side, d()));
        }
        long int x, y;
        div_t n;
        for (unsigned long int i = 0; i < int_pow(lattice_side, d()); i++) {
            map_of_contacts_int[ndim() * i] = i + 1;
            map_of_contacts_int[ndim() * i + 1] = i - 1;
            map_of_contacts_int[ndim() * i + 2] = i + lattice_side;
            map_of_contacts_int[ndim() * i + 3] = i - lattice_side;
            n = div(i, lattice_side);
            x = n.rem;
            y = n.quot;
            if (x == 0) {
                map_of_contacts_int[ndim() * i + 1] = i + lattice_side - 1;
            }
            if (x == (lattice_side - 1)) {
                map_of_contacts_int[ndim() * i] = i - (lattice_side - 1);
            }
            if (y == 0) {
                map_of_contacts_int[ndim() * i + 3] = lattice_side * (lattice_side - 1) + x;
            }
            if (y == (lattice_side - 1)) {
                map_of_contacts_int[ndim() * i + 2] = x;
            }
            coords[0][i] = x;
            coords[1][i] = y;
        }
	}
};

class Square_Lattice_3D : public Lattice {
    int d() override { return 3; }
    int ndim() override { return 6; }
    void create_lattice(long int max_seq_size = 100) override {
        lattice_side = max_seq_size;
        map_of_contacts_int.resize(int_pow(lattice_side, d()) * ndim());
        coords.resize(d());
        for (int i = 0; i < d(); i++) {
            coords[i].resize(int_pow(lattice_side, d()));
        }
        long int x, y, z;
        div_t n;
        for (unsigned long int i = 0; i < int_pow(lattice_side, d()); i++) {
            map_of_contacts_int[ndim() * i] = i + 1;
            map_of_contacts_int[ndim() * i + 1] = i - 1;
            map_of_contacts_int[ndim() * i + 2] = i + lattice_side;
            map_of_contacts_int[ndim() * i + 3] = i - lattice_side;
            map_of_contacts_int[ndim() * i + 4] = i + lattice_side * lattice_side;
            map_of_contacts_int[ndim() * i + 5] = i - lattice_side * lattice_side;
            n = div(i, lattice_side * lattice_side);
            z = n.quot;
            x = n.rem % lattice_side;
            y = n.rem / lattice_side;
            if (x == 0) {
                map_of_contacts_int[ndim() * i + 1] = i + lattice_side - 1;
            }
            if (x == (lattice_side - 1)) {
                map_of_contacts_int[ndim() * i] = i - (lattice_side - 1);
            }
            if (y == 0) {
                map_of_contacts_int[ndim() * i + 3] = lattice_side * (lattice_side - 1) + i;
            }
            if (y == (lattice_side - 1)) {
                map_of_contacts_int[ndim() * i + 2] = i - lattice_side * (lattice_side - 1);
            }
            if (z == 0) {
                map_of_contacts_int[ndim() * i + 5] = i + lattice_side * lattice_side * (lattice_side - 1);
            }
            if (z == lattice_side - 1) {
                map_of_contacts_int[ndim() * i + 4] = i - lattice_side * lattice_side * (lattice_side - 1);
            }
            coords[0][i] = x;
            coords[1][i] = y;
            coords[2][i] = z;
        }
    }
};

class Triangle_Lattice_2D : public Lattice {
    int d() override { return 2; }
    int ndim() override { return 6; }
    void create_lattice(long int max_seq_size = 100) override {
        lattice_side = max_seq_size;
        map_of_contacts_int.resize(int_pow(lattice_side, d()) * ndim());
        coords.resize(d());
        for (int i = 0; i < d(); i++) {
            coords[i].resize(int_pow(lattice_side, d()));
        }
        lattice_side = max_seq_size;
        //создается одномерный массив соседей на квадратной решетке
        long int x, y;
        div_t n;
        for (unsigned long int i = 0; i < int_pow(lattice_side, d()); i++) {
            map_of_contacts_int[ndim() * i] = i + 1;
            map_of_contacts_int[ndim() * i + 1] = i - 1;
            map_of_contacts_int[ndim() * i + 2] = i + lattice_side;
            map_of_contacts_int[ndim() * i + 3] = i - lattice_side;
            map_of_contacts_int[ndim() * i + 4] = i + lattice_side + 1;
            map_of_contacts_int[ndim() * i + 5] = i - lattice_side - 1;
            n = div(i, lattice_side);
            x = n.rem;
            y = n.quot;
            for (int j = 0; j < ndim(); j++) {
                if (x == 0) {
                    map_of_contacts_int[ndim() * i + 1] = i + lattice_side - 1;
                }
                if (x == (lattice_side - 1)) {
                    map_of_contacts_int[ndim() * i] = i - (lattice_side - 1);
                }
                if (y == 0) {
                    map_of_contacts_int[ndim() * i + 3] = lattice_side * (lattice_side - 1) + x;
                }
                if (y == (lattice_side - 1)) {
                    map_of_contacts_int[ndim() * i + 2] = x;
                }
                if ((y == (lattice_side - 1)) && (x == (lattice_side - 1))) {
                    map_of_contacts_int[ndim() * i + 4] = 0;
                }
                if ((y == (lattice_side - 1)) && (x < (lattice_side - 1))) {
                    map_of_contacts_int[ndim() * i + 4] = x + 1;
                }
                if ((y < (lattice_side - 1)) && (x == (lattice_side - 1))) {
                    map_of_contacts_int[ndim() * i + 4] = i + 1;
                }
                if ((y == 0) && (x == 0)) {
                    map_of_contacts_int[ndim() * i + 5] = lattice_side * lattice_side - 1;
                }
                if ((y == 0) && (x > 0)) {
                    map_of_contacts_int[ndim() * i + 5] = lattice_side * (lattice_side - 1) + x - 1;
                }
                if ((y > 0) && (x == 0)) {
                    map_of_contacts_int[ndim() * i + 5] = i - 1;
                }
                coords[0][i] = x;
                coords[1][i] = y;
            }
        }
    }
};

class Square_Lattice_4D : public Lattice {
    int d() override { return 4; }
    int ndim() override { return 8; }
    void create_lattice(long int max_seq_size = 100) override {
        lattice_side = max_seq_size;
        map_of_contacts_int.resize(int_pow(lattice_side, d()) * ndim());
        coords.resize(d());
        for (int i = 0; i < d(); i++) {
            coords[i].resize(int_pow(lattice_side, d()));
        }
        unsigned short x, y, z, t;
        div_t n;
        for (coord_t i = 0; i < int_pow(lattice_side, d()); i++) {
            map_of_contacts_int[ndim() * i] = i + 1;
            map_of_contacts_int[ndim() * i + 1] = i - 1;
            map_of_contacts_int[ndim() * i + 2] = i + lattice_side;
            map_of_contacts_int[ndim() * i + 3] = i - lattice_side;
            map_of_contacts_int[ndim() * i + 4] = i + lattice_side * lattice_side;
            map_of_contacts_int[ndim() * i + 5] = i - lattice_side * lattice_side;
            map_of_contacts_int[ndim() * i + 6] = i + lattice_side * lattice_side * lattice_side;
            map_of_contacts_int[ndim() * i + 7] = i - lattice_side * lattice_side * lattice_side;
            n = div(i, lattice_side * lattice_side * lattice_side);
            t = n.quot;
            n = div(n.rem, lattice_side * lattice_side);
            z = n.quot;
            n = div(n.rem, lattice_side);
            x = n.rem;
            y = n.quot;
            if (x == 0) {
                map_of_contacts_int[ndim() * i + 1] = i + lattice_side - 1;
            }
            if (x == (lattice_side - 1)) {
                map_of_contacts_int[ndim() * i] = i - (lattice_side - 1);
            }
            if (y == 0) {
                map_of_contacts_int[ndim() * i + 3] = lattice_side * (lattice_side - 1) + i;
            }
            if (y == (lattice_side - 1)) {
                map_of_contacts_int[ndim() * i + 2] = i - lattice_side * (lattice_side - 1);
            }
            if (z == 0) {
                map_of_contacts_int[ndim() * i + 5] = i + lattice_side * lattice_side * (lattice_side - 1);
            }
            if (z == lattice_side - 1) {
                map_of_contacts_int[ndim() * i + 4] = i - lattice_side * lattice_side * (lattice_side - 1);
            }
            if (t == 0) {
                map_of_contacts_int[ndim() * i + 7] = i + lattice_side * lattice_side * lattice_side * (lattice_side - 1);
            }
            if (t == lattice_side - 1) {
                map_of_contacts_int[ndim() * i + 6] = i - lattice_side * lattice_side * lattice_side * (lattice_side - 1);
            }
            coords[0][i] = x;
            coords[1][i] = y;
            coords[2][i] = z;
            coords[3][i] = t;
        }
    }
};