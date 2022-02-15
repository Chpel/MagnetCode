#include "mcmc.h"
#include <iostream>
#include <fstream>
#include <cassert>

//Lattice::Lattice() {};
std::random_device generator;
std::uniform_real_distribution<double> distribution(0.0, 1.0);
//std::mt19937 generator(1234);
std::random_device generators1;
//std::mt19937 generators1(123);
std::random_device generators2;
std::uniform_int_distribution<int> distribution2(0, 1);
//std::mt19937 generators3(12377);
std::random_device generators3;

Protein::Protein() {}

Protein::Protein(long int n, std::string l, double k_u = 1) {
	double metric_exp = 1;
    /** type | previous | next   **/
    if (l == "square") {
        lattice = new Square_Lattice_2D;
    }
    else if (l == "triangle") {
        lattice = new Triangle_Lattice_2D;
    }
    else if (l == "cubic") {
        lattice = new Square_Lattice_3D;
    }
    else if (l == "hypercubic") {
        lattice = new Square_Lattice_4D;
		double metric_exp = 0.5;
    }
    this->l = l;
	int side = (int)(k_u * pow(n, metric_exp)) + 5;
	std::cout << "Creating " << l << " lattice with side = "<< side << std::endl;
    lattice->create_lattice(side); //создание решетки, на 5 больше, чем длина цепочки
    std::cout << l << " lattice ready! Len = " << lattice->map_of_contacts_int.size() << ". Metric exponent = " << metric_exp <<std::endl;
    //массив числа соседей
    bulk.resize(lattice->ndim() - 1);
    bulk_now.resize(lattice->ndim() - 1);
    // 0 - отстуттвие элемента на решетке
    
    size_t len = int_pow(lattice->lattice_side, lattice->d());
    std::cout << "Additional arrays resizing. Len = " << len << std::endl;
    sequence_on_lattice.resize(len, 0); //последовательность мономеров
    next_monomers.resize(len, -1); //номер-ссылка на следующий узел
    previous_monomers.resize(len, -1); //номер-ссылка на предыдующий узел
    directions.resize(len, -1); //направления из {0,1,2,3}
    //ordered_coords.resize(len, -1); //пока что для последовательной нумерации

    number_of_monomers = n;
    start_conformation = 0;
	coord_check_t current_pos = 0;
	int step = 0;
	//end_conformation = n - 1;

    sum_coord.resize(lattice->d());
	short step_on_map[] =  {1, -1, lattice->lattice_side};
										
	for (int i = 1; i < n - 1; i++)
    {
		next_monomers[current_pos] = current_pos + step_on_map[step];
		directions[current_pos] = step;
		current_pos += step_on_map[step];
		sequence_on_lattice[current_pos] = 1;
		previous_monomers[current_pos] = current_pos - step_on_map[step];
		
		sum_coord[0] += lattice->coords[0][current_pos];
		sum_coord[1] += lattice->coords[1][current_pos];
		
		if ((step == 0) && (lattice->coords[0][current_pos] == lattice->lattice_side - 2)) {step = 2;}
		else if ((step == 2) && (lattice->coords[1][current_pos] % 2 == 0)) {
			if (lattice->coords[0][current_pos] == 0) {step = 0;}
			else {step = 1;}
		}
		else if ((step == 1) && (lattice->coords[0][current_pos] == 0)) {step = 2;}
		//ordered_coords[i] = ;
    }		
	next_monomers[current_pos] = current_pos + step_on_map[step];
	directions[current_pos] = step;
	current_pos += step_on_map[step];
	sequence_on_lattice[current_pos] = 1;
	previous_monomers[current_pos] = current_pos - step_on_map[step];
	
	sum_coord[0] += lattice->coords[0][current_pos];
	sum_coord[1] += lattice->coords[1][current_pos];
	
    //ordered_coords[0] = 0;
    //ordered_coords[n - 1] = n - 1;

    sequence_on_lattice[0] = 1;
	end_conformation = current_pos;

    current_H_counts = n;
    E = -(n - 1);


    //сначала все направления - движение вправо
    

    bulk_now[0] = number_of_monomers - 2;

	std::cout << l << " lattice is ready!" << std::endl;
	
	save_walks();
	exit(1);
	
}

void Protein::count_contacts()
{
    long int hh = 0;
    coord_check_t current_position = start_conformation;
    coord_t step;
    long int mag = 0;
    for (int i = 0; i < number_of_monomers; i++) {
        for (int j = 0; j < lattice->ndim(); j++) {
            step = lattice->map_of_contacts_int[lattice->ndim() * current_position + j];
            if (sequence_on_lattice[step] != 0)
            {
                hh = hh + sequence_on_lattice[current_position] * sequence_on_lattice[step];
            }
        }
        mag = mag + sequence_on_lattice[current_position];
        current_position = next_monomers[current_position];
    }

    E = -(hh / 2);
    current_H_counts = mag;
}


bool Protein::IsEndInStuck()
{
    int hh = 0;
    coord_t step;
    for (int j = 0; j < lattice->ndim(); j++) {
        step = lattice->map_of_contacts_int[lattice->ndim() * end_conformation + j];
        if (sequence_on_lattice[step] != 0) {
            hh = hh + 1;
        }
    }
    return hh == lattice->ndim();
}

/*void Protein::Reconnect1(int j) {
    int inverse_steps[4] = { 1, 0, 3, 2 };
    int reflect_directions[4][4] =
    { {3, 2, 0, 1}, //90
     {1, 0, 3, 2}, //180
     {2, 3, 1, 0}, //270
     {0, 1, 2, 3}
    };

    //int j = directions[previous_monomers[end_conformation]]; //направление последнего ребра

    long int step = lattice->map_of_contacts_int[lattice->ndim() * end_conformation + j];
    long int new_end = next_monomers[step];

    next_monomers[step] = end_conformation;

    directions[step] = inverse_steps[j];
    directions[new_end] = -1;

    long int c = end_conformation;
    long int new_c;
    while (c != new_end)
    {
        new_c = previous_monomers[c];

        next_monomers[c] = previous_monomers[c];
        //previous_monomers[new_c]=c;
        directions[c] = inverse_steps[directions[new_c]];
        c = new_c;
    }
    long int temp_prev_next = next_monomers[new_end];

    previous_monomers[end_conformation] = step;
    c = end_conformation;

    while (c != new_end)
    {
        new_c = next_monomers[c];
        previous_monomers[new_c] = c;
        c = new_c;
    }

    end_conformation = new_end;


    previous_monomers[new_end] = temp_prev_next;
    next_monomers[new_end] = -1;

}*/

void Protein::calc_bulk()
{
    for (int dm = 2; dm <= lattice->ndim(); dm++) {
        bulk_now[dm - 2] = 0;
    }

    coord_check_t current = start_conformation;
    coord_t step;
    int k = 0;
    for (int e = 0; e < number_of_monomers; e++)
    {
        k = 0;
        for (int j = 0; j < lattice->ndim(); j++) {
            step = lattice->map_of_contacts_int[lattice->ndim() * current + j];
            if (sequence_on_lattice[step] != 0) {
                k += 1;
            }
        }

        if (k > 1) {
            bulk_now[k - 2] += 1;
        }

        

        current = next_monomers[current];

    }


    //std::cout << current << " ";

}





void Protein::Reconnect(int j) {
    std::vector<int> inverse_steps;
    inverse_steps.resize(lattice->ndim());
    for (int dm = 0; dm < lattice->ndim(); dm++) {
        if (dm % 2 == 0) {
            inverse_steps[dm] = dm + 1;
        }
        else {
            inverse_steps[dm] = dm - 1;
        }
    }
    //    int reflect_directions[4][4] =
    //            {{3, 2, 0, 1}, //90
    //             {1, 0, 3, 2}, //180
    //             {2, 3, 1, 0}, //270
    //             {0, 1, 2, 3}
    //            };
    coord_check_t c;

    coord_check_t step = lattice->map_of_contacts_int[lattice->ndim() * end_conformation + j];
    coord_check_t new_end = next_monomers[step];

    next_monomers[step] = end_conformation;

    directions[step] = inverse_steps[j];

    c = end_conformation;

    //std::cout << " end_conformation " <<  end_conformation << std::endl;
    if (end_conformation > int_pow(lattice->lattice_side, lattice->d()))
    {
        std::cout << " HZ why " << std::endl;
        return;

    }
    coord_check_t new_c;
    while (c != new_end)
    {
        // std::cout << "c  " << c << std::endl;
        assert(c >= 0 && c < int_pow(lattice->lattice_side, lattice->d()));
        if (c > int_pow(lattice->lattice_side, lattice->d()))
        {
            std::cout << " HZ why " << std::endl;
            return;
        }

        if (previous_monomers[c] > int_pow(lattice->lattice_side, lattice->d()))
        {
            std::cout << "nor  prev " << std::endl;
        }

        new_c = previous_monomers[c];

        if (new_c < 0 || c > int_pow(lattice->lattice_side, lattice->d()))
        {
            std::cout << " we have problems " << std::endl;
        }

        next_monomers[c] = previous_monomers[c];


        directions[c] = inverse_steps[directions[new_c]];
        c = new_c;
    }
    coord_check_t temp_prev_next = next_monomers[new_end];


    previous_monomers[end_conformation] = step;
    c = end_conformation;


    while (c != new_end)
    {
        new_c = next_monomers[c];
        previous_monomers[new_c] = c;
        c = new_c;
    }

    end_conformation = new_end;


    //std::cout <<"mmmm" << " " << step << std::endl;

    previous_monomers[new_end] = temp_prev_next;
    next_monomers[new_end] = -1;
    directions[new_end] = -1;
}


void Protein::MC(double J_in, double h_in, int Simulation, long int steps_to_equilibrium, long int mc_steps, bool bradius)
{
	
    std::uniform_int_distribution<int> distribution1(0, lattice->ndim() - 1);
    nSimulation = Simulation;
    J = J_in;
    h = h_in;
	
	k_too_big = 0;
    //[i,j]: i-поворот, j-направление (против часовой)
    int reflect_directions[4][4] =
    { {3, 2, 0, 1 }, //90
     {1,0,3,2 }, //180
     {2,3, 1, 0 }, //270
    { 0, 1, 2, 3}
    };

    std::vector<int> inverse_steps;
    inverse_steps.resize(lattice->ndim());
    for (int dm = 0; dm < lattice->ndim(); dm++) {
        if (dm % 2 == 0) {
            inverse_steps[dm] = dm + 1;
        }
        else {
            inverse_steps[dm] = dm - 1;
        }
    }

    //double step_rd; //Для выбора апдейта: обычный или реконнект
    double q_rd, p1, p_metropolis; //Для вероятности принятия шага
    int rand_path; // = distribution1(generators1); //выбирается направление: 0 - переставляем начало в конец
    double typeOfUpdate; //0 - простой; 1 - реконнект
    coord_t step;
    int step_on_lattice;//выбор одного из соседей
    coord_t new_point; //номер новой точки в цепочке
    long int new_E, new_H;
    int hh;
    coord_check_t temp, del;
	char oldspin;


    //std::uniform_int_distribution<long int> distribution_spin(0, number_of_monomers-1);
    std::uniform_int_distribution<coord_t> distribution_spin(0, int_pow(lattice->lattice_side, lattice->d()) - 1);
    //std::mt19937 generator_spin(123);
    std::random_device generator_spin;


    //вероятность добавить спин в кластер в кластерном апдейте
    P_add = 1 - exp(-2 * J); //пока так для h=0

    double p_for_local_update = 0.6;
    double p_for_reconnect = 1.0; //p_for_local_update - p_for_reconnect = вероятность реконнекта

    //spins_in_cluster.resize(number_of_monomers, false);

    long int all_steps = steps_to_equilibrium + mc_steps;

    for (long int i = 0; i < all_steps + 2; i++) {
        //std::cout << "STEP : " << i << std::endl;
        typeOfUpdate = distribution(generator);
        if (typeOfUpdate < p_for_local_update) {
            hh = 0;
            rand_path = distribution2(generators3);

            if (rand_path == 0) {//переставляем начало в конец

                step_on_lattice = distribution1(generators1);
                new_point = lattice->map_of_contacts_int[lattice->ndim() * end_conformation + step_on_lattice];
                oldspin = sequence_on_lattice[start_conformation];

                if (sequence_on_lattice[new_point] == 0) { //проверка, что в узле нет мономеров

                    //делаем апдейт

                    //добавляем в конец
                    next_monomers[end_conformation] = new_point;
                    sequence_on_lattice[new_point] = 2 * distribution2(generators3) - 1; //выбор спина
                    previous_monomers[new_point] = end_conformation;
                    end_conformation = new_point;

                    //удаляем начало
                    temp = start_conformation;
                    start_conformation = next_monomers[start_conformation];
                    next_monomers[temp] = -1;
                    previous_monomers[start_conformation] = -1;

                    //смотрим потери
                    for (int j = 0; j < lattice->ndim(); j++) {
                        step = lattice->map_of_contacts_int[lattice->ndim() * temp + j];
                        if (sequence_on_lattice[step] != 0) {
                            hh = hh - sequence_on_lattice[temp] * sequence_on_lattice[step];
                        }
                    }

                    //смотрим выигрыш
                    for (int j = 0; j < lattice->ndim(); j++) {
                        step = lattice->map_of_contacts_int[lattice->ndim() * end_conformation + j];
                        if (sequence_on_lattice[step] != 0) {
                            hh = hh + sequence_on_lattice[end_conformation] * sequence_on_lattice[step];
                        }
                    }

                    new_E = E + hh;
                    //new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[start_conformation];
                    new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[temp];

                    p1 = exp(-(-(new_E - E) * J - (new_H - current_H_counts) * h));
                    p_metropolis = std::min(1.0, p1);
                    q_rd = distribution(generator);
                    if (q_rd < p_metropolis) {
                        E = new_E;
                        current_H_counts = new_H;
                        sequence_on_lattice[temp] = 0; //делаю здесь, так как проще считать энергию(!!!)
                        for (int dm = 0; dm < lattice->d(); dm++) {
                            sum_coord[dm] = sum_coord[dm] + lattice->coords[dm][end_conformation] - lattice->coords[dm][temp];
                        }
                        //корректируем информацию о направлениях
                        directions[temp] = -1;
                        directions[previous_monomers[end_conformation]] = step_on_lattice;

                        //перенумерация
                        //ordered_coords[temp] = -1;
                        //ordered_coords[new_point] = end_conformation;

                    }
                    else {//отменяем изменения
                        //удаляем конец
                        del = end_conformation;
                        end_conformation = previous_monomers[end_conformation];
                        next_monomers[end_conformation] = -1;
                        previous_monomers[del] = -1;
                        sequence_on_lattice[del] = 0;

                        //возвращаем начало
                        previous_monomers[start_conformation] = temp;
                        next_monomers[temp] = start_conformation;
                        start_conformation = temp;
                        sequence_on_lattice[start_conformation] = oldspin;

                    }

                }
                else {//места нет, выходим из шага
                    //continue;
                }

            }
            else {//переставляем конец в начало

                step_on_lattice = distribution1(generators1);
                new_point = lattice->map_of_contacts_int[lattice->ndim() * start_conformation + step_on_lattice];
                oldspin = sequence_on_lattice[end_conformation];

                if (sequence_on_lattice[new_point] == 0) { //проверка, что в узле нет мономеров

                    //делаем апдейт
                    //добавляем в начало
                    previous_monomers[start_conformation] = new_point;
                    sequence_on_lattice[new_point] = 2 * distribution2(generators3) - 1; //выбор спина
                    next_monomers[new_point] = start_conformation;
                    start_conformation = new_point;

                    //удаляем конец
                    temp = end_conformation;
                    end_conformation = previous_monomers[end_conformation];
                    if (previous_monomers[end_conformation] < 0) {
                        std::cout << "problem update " << std::endl;
                    }
                    previous_monomers[temp] = -1;
                    next_monomers[end_conformation] = -1;

                    //смотрим потери
                    for (int j = 0; j < lattice->ndim(); j++) {
                        step = lattice->map_of_contacts_int[lattice->ndim() * temp + j];
                        if (sequence_on_lattice[step] != 0) {
                            hh = hh - sequence_on_lattice[temp] * sequence_on_lattice[step];
                        }
                    }

                    //смотрим выигрыш
                    for (int j = 0; j < lattice->ndim(); j++) {
                        step = lattice->map_of_contacts_int[lattice->ndim() * start_conformation + j];
                        if (sequence_on_lattice[step] != 0) {
                            hh = hh + sequence_on_lattice[start_conformation] * sequence_on_lattice[step];
                        }
                    }

                    new_E = E + hh;
                    new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[temp];

                    //p1 = exp(-(new_E - E) * J - (new_H - current_H_counts) * h);

                    p1 = exp(-(-(new_E - E) * J - (new_H - current_H_counts) * h));
                    p_metropolis = std::min(1.0, p1);
                    q_rd = distribution(generator);

                    if (q_rd < p_metropolis) {
                        E = new_E;
                        current_H_counts = new_H;
                        sequence_on_lattice[temp] = 0; //делаю здесь, так как проще считать энергию(!!!)
                        for (int dm = 0; dm < lattice->d(); dm++) {
                            sum_coord[dm] = sum_coord[dm] + lattice->coords[dm][end_conformation] - lattice->coords[dm][temp];
                        }
                        //sum_X = sum_X + lattice->x_coords[start_conformation];
                        //sum_Y = sum_Y + lattice->y_coords[start_conformation];

                        //корректируем информацию о направлениях
                        directions[end_conformation] = -1;
                        directions[start_conformation] = inverse_steps[step_on_lattice];

                    }
                    else {//отменяем изменения
                     //удаляем начало
                        del = start_conformation;
                        start_conformation = next_monomers[start_conformation];
                        previous_monomers[start_conformation] = -1;
                        next_monomers[del] = -1;
                        sequence_on_lattice[del] = 0;

                        //возвращаем конец
                        next_monomers[end_conformation] = temp;
                        previous_monomers[temp] = end_conformation;
                        end_conformation = temp;
                        sequence_on_lattice[end_conformation] = oldspin;

                        if (previous_monomers[temp] < 0) {
                            std::cout << "problem return " << std::endl;
                        }
                        if (temp < 0) {
                            std::cout << "problem return temp" << std::endl;
                        }
                    }
                }
                else {
                    //некуда идти
                }

            }
        }
        else {

            coord_t coord = distribution_spin(generator_spin);

            if (sequence_on_lattice[coord] != 0)
            { //вероятность такого события 1/n, делаем кластерный апдейт
                int sign = sequence_on_lattice[coord];

                std::valarray<bool> used_coords;
                used_coords.resize(int_pow(lattice->lattice_side,lattice->d()), false);

                std::queue<coord_t> Cluster;

                Cluster.push(coord);
                used_coords[coord] = true;

                while (!Cluster.empty()) {
                    temp = Cluster.front();
                    Cluster.pop();

                    for (int j = 0; j < lattice->ndim(); j++)
                    {
                        step = lattice->map_of_contacts_int[lattice->ndim() * temp + j];
                        double p = distribution(generator);
                        //???
                        if (sequence_on_lattice[step] == sign && p < P_add &&
                            !used_coords[step]) {
                            Cluster.push(step);
                            used_coords[step] = true;
                            sequence_on_lattice[step] *= -1;
                        }
                    }
                }
                sequence_on_lattice[coord] *= -1;

                count_contacts();

                //std::cout << "cluster type done " << std::endl;
            }
            else
            { //обычно будем здесь, здесь попытка сделать реконнект
                step_on_lattice = distribution1(generators1);
                new_point = lattice->map_of_contacts_int[lattice->ndim() * end_conformation + step_on_lattice];

                //проверка, что проверенный узел занят спином
                if (sequence_on_lattice[new_point] != 0 && next_monomers[new_point] != -1 &&
                    new_point != previous_monomers[end_conformation]) {
                    Reconnect(step_on_lattice);
                }
            }


        }

        //ИЛЬЯ: подсчёт значений (пока без coord_form)
        if (i > steps_to_equilibrium && i % 100000 == 0)
        {
            save_calcs();
            //radius();


            calc_bulk();
            for (int dm = 2; dm <= lattice->ndim(); dm++) {
                bulk[dm - 2] << 1.0 * bulk_now[dm - 2] / number_of_monomers;
            }

            //coord_form();

        }

		if (i > steps_to_equilibrium) {
			k_too_big += isTooBig();
		}

        
        if (i > steps_to_equilibrium && i % 100000000 == 0)
        {

            write_file(i);

			//save_walks();
        }


    }

}

int Protein::isTooBig() {
	std::vector <long long int> coords;
	std::vector <long long int> min_coords;
	std::vector <long long int> max_coords;
    coords.resize(lattice->d(),0);
    min_coords.resize(lattice->d(), 0);
    max_coords.resize(lattice->d(), 0);

    int direction;

    long int current = start_conformation;


    std::vector <std::vector<int>> steps = { {1,0,0,0}, {-1,0,0,0}, {0,1,0,0}, {0,-1,0,0},
											 {0,0,1,0}, {0,0,-1,0}, {0,0,0,0}, {0,0,0,-1}};
	long long int x;
    for (long long int i = 1; i < number_of_monomers; i++)
    {
        direction = directions[current];
		for (int j = 0; j < lattice->d(); j++) {
			x = coords[j] + steps[direction][j];
			if (x > max_coords[j]) {max_coords[j] = x;}
			if (x < min_coords[j]) {min_coords[j] = x;}
			coords[j] = x;
        }

        current = next_monomers[current];
    }
	
	for (int j = 0; j < lattice->d(); j++) {
		if ((max_coords[j] - min_coords[j]) > lattice->lattice_side) return 1;
	}
	return 0;
	
}

/*
bool Protein::CheckAndFlipNodes (long int& coord, int& sign )
{
    for (int j = 0; j < lattice->ndim2(); j++) {
     long int step = lattice->map_of_contacts_int[lattice->ndim2() * coord + j];
       // std::cout << coord << "  start "<< step <<  std::endl;
        if (sequence_on_lattice[step] != 0) {
            if (sequence_on_lattice[step] ==sign )
            {
                double p = distribution(generator);
                if (p < P_add)
                {
                    CheckAndFlipNodes (step, sign );
                    spins_in_cluster.push(step);
                    //Можно идти в глубину?
                    //spins_in_cluster[]
                }
            }
        }
    }
   // std::cout << coord << " finish function  " <<  std::endl;
return false;
}*/

/*void Protein::radius()
{
    long int point1x = end_conformation % lattice->lattice_side;
    long int point1y = end_conformation / lattice->lattice_side;
    long int point1xs = start_conformation % lattice->lattice_side;
    long int point1ys = start_conformation / lattice->lattice_side;
  //расстояние на торе
    long int xdiff = abs(point1x- point1xs);
    if (xdiff > (lattice->lattice_side  / 2))
    xdiff = lattice->lattice_side - xdiff;
    long int ydiff = abs(point1y- point1ys);
    if (ydiff > (lattice->lattice_side / 2))
    ydiff = lattice->lattice_side - ydiff;
    long int r = xdiff *xdiff  + ydiff*ydiff;
    dists << r;
}*/

long int Protein::radius()
{
    long int point1_c;
    long int point1s_c;
    long int diff;
    long int r = 0;
    //расстояние на торе
    for (int dm = 0; dm < lattice->d(); dm++) {
        point1_c = lattice->coords[dm][end_conformation];
        point1s_c = lattice->coords[dm][start_conformation];

        diff = abs(point1_c - point1s_c);
        if (diff > (lattice->lattice_side / 2))
            diff = lattice->lattice_side - diff;
        r += diff * diff;
    }
    
    dists << r;
    return r;
}



void Protein::radius_gyration()
{


    long double r_g = 0;
    coord_check_t current = start_conformation;
    long double y = 0, x = 0;
    long double point1x = 0, point1y = 0;
    std::vector<long int> point1;
    point1.resize(lattice->d(), 0);
    //long double point1x = 1.0*sum_X/number_of_monomers;
    //long double point1y = 1.0*sum_Y/number_of_monomers;
    long double xdiff, ydiff;
    std::vector<long double> diff;
    diff.resize(lattice->d(), 0);
    //point1x = start_conformation % lattice->lattice_side;
    //point1y = start_conformation / lattice->lattice_side;

    current = start_conformation;
    long double r;
    long int point1xs, point1ys;
    std::vector<long int> point1s;
    point1s.resize(lattice->d(), 0);
    coord_check_t second_current;
    for (int e = 0; e < number_of_monomers; e++)
    {
        second_current = start_conformation;

        for (int dm = 0; dm < lattice->d(); dm++) {
            point1[dm] = lattice->coords[dm][current];
        }

        for (int e1 = 0; e1 < number_of_monomers; e1++) {

            for (int dm = 0; dm < lattice->d(); dm++) {
                point1s[dm] = lattice->coords[dm][second_current];
            }
            
            //расстояние на торе
            r = 0;
            for (int dm = 0; dm < lattice->d(); dm++) {
                diff[dm] = abs(point1[dm] - point1s[dm]);
                if (diff[dm] > (lattice->lattice_side / 2))
                    diff[dm] = lattice->lattice_side - diff[dm];

                r += diff[dm] * diff[dm];
            }
            r_g = r_g + r;

            second_current = next_monomers[second_current];

        }
        current = next_monomers[current];

        //std::cout << current << " ";

    }
    //std::cout << std::endl;
    gyration << 0.5 * r_g / number_of_monomers / number_of_monomers;

}

void Protein::save_walks() {
	
	std::string filename;
    std::ofstream out_result;

    filename = "Ising_Walk_" + l + " " + std::to_string(J) + "_" + std::to_string(h) + "_" + std::to_string(number_of_monomers) + ".txt";
    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
	
	size_t current = start_conformation;
	for (size_t i = 0; i < number_of_monomers; i++) {
		for (size_t k = 0; k < lattice->d(); k++) {
			out_result << lattice->coords[k][current] << " ";
		}
		out_result << std::endl;
		current = next_monomers[current];
	}
	
	out_result.close();
}

//ИЛЬЯ: пока я не решил проблему для собственных значений при разной размерности. пока лежит так

/*void Protein::coord_form() {

    std::vector <long long int> xs, ys;
    xs.push_back(0);
    ys.push_back(0);

    int direction;

    long int current = start_conformation;


    std::vector <std::vector<int>> steps = { {1,0}, {-1,0}, {0,1}, {0,-1} };

    for (long long int i = 1; i < number_of_monomers; i++)
    {
        direction = directions[current];

        long long int x = xs.back() + steps[direction][0];
        long long int y = ys.back() + steps[direction][1];

        xs.push_back(x);
        ys.push_back(y);

        current = next_monomers[current];


    }

    long double r_g = 0;
    long double r, r_notdig;
    long double y = 0, x = 0;
    long double point1x = 0, point1y = 0;
    //long double point1x = 1.0*sum_X/number_of_monomers;
    //long double point1y = 1.0*sum_Y/number_of_monomers;
    long double xdiff, ydiff;
    double  A_element = 0, D_element = 0, BC_element2 = 0;


    long long int second_current;
    for (int e = 0; e < number_of_monomers; e++)
    {
        second_current = start_conformation;



        for (int e1 = 0; e1 < number_of_monomers; e1++) {


            //расстояние на торе

            xdiff = xs[e1] - xs[e];
            ydiff = ys[e1] - ys[e];

            r = xdiff * xdiff + ydiff * ydiff;
            r_g = r_g + r;



            A_element += xdiff * xdiff;
            D_element += ydiff * ydiff;
            BC_element2 += xdiff * ydiff;
            //r_notdig = xdiff*ydiff;
            //NoDiagnalEleents = NoDiagnalEleents + r_notdig;



            second_current = next_monomers[second_current];

        }
        current = next_monomers[current];


        //std:: cout << eigs1.mean() << " " << eigs2.mean() << " " << aratio.mean()  << " " << gyration.mean() << std::endl;
        //std::cout << current << " ";

    }

    gyration << 0.5 * r_g / number_of_monomers / number_of_monomers;
    A_element = 1.0 * A_element / number_of_monomers / number_of_monomers / 2.0;
    D_element = 1.0 * D_element / number_of_monomers / number_of_monomers / 2.0;

    BC_element2 = 1.0 * BC_element2 / number_of_monomers / number_of_monomers / 2.0;

    double D = (A_element + D_element) * (A_element + D_element) - 4 * (A_element * D_element - BC_element2 * BC_element2);

    double eigvals1 = ((A_element + D_element) + sqrt(D)) * 0.5;

    double eigvals2 = ((A_element + D_element) - sqrt(D)) * 0.5;

    eigs1 << eigvals1;
    eigs2 << eigvals2;

    aratio << 1.0 * (eigvals1 - eigvals2) * (eigvals1 - eigvals2) / ((eigvals1 + eigvals2) * (eigvals1 + eigvals2));



    //std::cout << " Compare from coords " << std:: endl;
    //std::cout << A_element+D_element <<  " - " << 0.5*r_g/number_of_monomers/number_of_monomers << std::endl ;
    //std:: cout << eigs1.mean() << " " << eigs2.mean() << " " << aratio.mean()  << " " << gyration.mean() << std::endl;
}*/

void Protein::write_file(long int i) {

    std::string filename;
    std::ofstream out_result;

    filename = "Geometry_Ising_" + l + " " + std::to_string(J) + "_" + std::to_string(h) + "_" + std::to_string(number_of_monomers) + ".txt";
    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";
    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq ";
    out_result << "lambda1 err_lambda1 lambda2 err_lambda2 acperical err_aspherical ";
    for (int dm = 2; dm <= lattice->ndim(); dm++) {
        out_result << "bulk" << dm << " err_bulk" << dm << " ";
    }
    out_result << " steps" << " k_too_big"<< std::endl;

    out_result << number_of_monomers << " " << J << " " << h << " ";
    out_result << dists.mean() << " " << dists.errorbar() << " " << gyration.mean() << " " << gyration.errorbar() << " ";

    out_result << eigs1.mean() << " " << eigs1.errorbar() << " ";
    out_result << eigs2.mean() << " " << eigs2.errorbar() << " ";
    out_result << aratio.mean() << " " << aratio.errorbar() << " ";

    for (int dm = 2; dm <= lattice->ndim(); dm++) {
        out_result << bulk[dm - 2].mean() << " " << bulk[dm - 2].errorbar() << " ";
    }

    out_result << i << " " << k_too_big;

    out_result << std::endl;

    out_result.close();


    out_result.close();



    filename = "BC_" + l + " " + std::to_string(J) + "_" + std::to_string(h) + "_" + std::to_string(number_of_monomers) + "_" + std::to_string(nSimulation) + ".txt";
    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";
    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq ";
    out_result << "mean_e err_mean_e mean_e_sq err_mean_e_sq mean_e_fourth err_mean_e_fourth ";
    out_result << "mean_m err_mean_m mean_m_sq err_mean_m_sq mean_m_fourth err_mean_m_fourth " << std::endl;

    out_result << number_of_monomers << " " << J << " " << h << " ";
    out_result << dists.mean() << " " << dists.errorbar() << " " << gyration.mean() << " " << gyration.errorbar() << " ";

    out_result << energy.mean() << " " << energy.errorbar() << " ";
    out_result << energy_sq.mean() << " " << energy_sq.errorbar() << " ";
    out_result << energy_4.mean() << " " << energy_4.errorbar() << " ";

    out_result << magnetization.mean() << " " << magnetization.errorbar() << " ";
    out_result << magnetization_sq.mean() << " " << magnetization_sq.errorbar() << " ";
    out_result << magnetization_4.mean() << " " << magnetization_4.errorbar() << " ";
    out_result << i << " ";

    out_result << std::endl;

    out_result.close();


    out_result.close();

}

//ИЛЬЯ: также пока лежит

/*void Protein::radius_gyration1()
{
    long double r_g = 0;
    long int current = start_conformation;
    long double y = 0, x = 0;
    long double point1x = 0, point1y = 0;
    std::vector<long double> point1;
    point1.resize(lattice->d(), 0);
    //long double point1x = 1.0*sum_X/number_of_monomers;
    //long double point1y = 1.0*sum_Y/number_of_monomers;
    long double xdiff, ydiff;
    point1x = start_conformation % lattice->lattice_side;
    point1y = start_conformation / lattice->lattice_side;

    for (int e = 0; e < number_of_monomers; e++)
    {

        long int point1xs = lattice->x_coords[current];
        long int point1ys = lattice->y_coords[current];

        //расстояние на торе
        xdiff = abs(point1x - point1xs);
        if (xdiff > (lattice->lattice_side / 2))
            xdiff = lattice->lattice_side - xdiff;

        ydiff = abs(point1y - point1ys);
        if (ydiff > (lattice->lattice_side / 2))
            ydiff = lattice->lattice_side - ydiff;

        x = x + xdiff;
        y = y + ydiff;

        //r = xdiff *xdiff  + ydiff*ydiff;

        current = next_monomers[current];

        //std::cout << current << " ";

    }

    point1x = x / number_of_monomers;
    point1y = y / number_of_monomers;

    current = start_conformation;
    long double r;
    long int point1xs, point1ys;
    for (int e = 0; e < number_of_monomers; e++)
    {

        point1xs = lattice->x_coords[current];
        point1ys = lattice->y_coords[current];

        //расстояние на торе
        xdiff = abs(point1x - point1xs);
        if (xdiff > (lattice->lattice_side / 2))
            xdiff = lattice->lattice_side - xdiff;

        ydiff = abs(point1y - point1ys);
        if (ydiff > (lattice->lattice_side / 2))
            ydiff = lattice->lattice_side - ydiff;

        r = xdiff * xdiff + ydiff * ydiff;

        r_g = r_g + r;


        current = next_monomers[current];

        //std::cout << current << " ";

    }
    //std::cout << std::endl;
    gyration << 1.0 * r_g / number_of_monomers;

}*/

void Protein::save_calcs()
{

    energy << 1.0 * (E) / number_of_monomers;
    energy_sq << 1.0 * (E) / number_of_monomers * 1.0 * (E) / number_of_monomers;
    energy_4 << 1.0 * (E) / number_of_monomers * 1.0 * (E) / number_of_monomers * 1.0 * (E) / number_of_monomers * 1.0 * (E) / number_of_monomers;

    magnetization << 1.0 * abs(current_H_counts) / number_of_monomers;
    magnetization_sq << 1.0 * current_H_counts / number_of_monomers * 1.0 * current_H_counts / number_of_monomers;
    magnetization_4 << 1.0 * current_H_counts / number_of_monomers * 1.0 * current_H_counts / number_of_monomers * 1.0 * current_H_counts / number_of_monomers * 1.0 * current_H_counts / number_of_monomers;


    count_E[E] = count_E[E] + 1;
    count_M[current_H_counts] = count_M[current_H_counts] + 1;

    radius();


}