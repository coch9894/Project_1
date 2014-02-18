#pragma once

struct Individual
{
	double fitness;
	double* vector;
};

class Population
{
public:
	Population(char fit);
	~Population(void);

	double min_range;
	double max_range;

	char func;

	int index;

	Individual* individual[1000];

	void Print_Ind_Vector(Individual* ind);
	void Print_Vector();
	void Print_Fitness();
	
	void Fitness();
	void Single_Fitness(Individual* ind);
	void Sort();

	void Tournament_Selection();
	void Roulette_Selection();

	void Mutate(Individual* ind[]);
	void Crossover();
};

