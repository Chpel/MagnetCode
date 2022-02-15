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

size_t int_pow(int lattice_side, int d);

typedef unsigned long long int coord_t;
typedef long long int coord_check_t;

class Lattice {
public:
	int lattice_side; //длина стороны решётки
	std::valarray<size_t> map_of_contacts_int;
	std::vector<std::valarray<short>> coords;
    virtual int ndim() = 0; //число соседей узла
    virtual int d() = 0; //размерность пространства
	virtual void create_lattice(long int max_seq_size = 100) = 0;
};

class Square_Lattice_2D : public Lattice {
	int ndim() override;
	int d() override;
	void create_lattice(long int max_seq_size = 100) override;
};

class Square_Lattice_3D : public Lattice {
	int ndim() override;
	int d() override;
	void create_lattice(long int max_seq_size = 100) override;
};

class Triangle_Lattice_2D : public Lattice {
	int ndim() override;
	int d() override;
	void create_lattice(long int max_seq_size = 100) override;
};

class Square_Lattice_4D : public Lattice {
	int ndim() override;
	int d() override;
	void create_lattice(long int max_seq_size = 100) override;
};