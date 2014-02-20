#include "StdAfx.h"
#include "Population.h"
#include <time.h>

int partition( Individual* input[], int left, int right )
{
	double pivot = abs(input[right]->fitness);

	while( left < right )
	{
		while( abs(input[left]->fitness) < pivot )
			left++;
		while( abs(input[right]->fitness) > pivot )
			right--;
		if( input[left]->fitness == input[right]->fitness )
			left++;
		else
		{
			Individual* temp = input[left];
			input[left] = input[right];
			input[right] = temp;
		}
	}
	return right;
}

void quicksort( Individual* input[], int left, int right )
{
    if ( left < right )
    {
        int j = partition(input, left, right);      
        quicksort(input, left, j-1);
        quicksort(input, j+1, right);
    }
}

Population::Population(char fit)
{
	func = fit;
	// set range
	// Sphere || Rastrigin
	if( func == 's' || func == 'i' )
	{
		min_range = -5.12;
		max_range = 5.12;
	}
	// Rosenbrock
	if( func == 'r' )
	{
		min_range = -2.048;
		max_range = 2.048;
	}
	//Schwefel
	if( func == 'c' )
	{
		min_range = -512.03;
		max_range = 511.97;
	}
	//Ackely
	if( func == 'a' )
	{
		min_range = -30;
		max_range = 30;
	}
	//Griewangk
	if( func == 'g' )
	{
		min_range = -600;
		max_range = 600;
	}
	ratio = 0.01;
	// generate 1000 random individuals
	srand(time(NULL));
	for( int i = 0; i < POP_SIZE; i++ )
	{
		individual[i] = new Individual;
		individual[i]->vector = new double[30];
		for( int j = 0; j < 30; j++ )
		{
			
			double random = ((double)rand()) / (double)RAND_MAX;
			double x = random * (max_range - min_range);
			individual[i]->vector[j] = min_range + x;
			//cout << individual[i]->vector[j] << endl;
			//individual[i]->vector[j] = -420;
		}
	}
	Fitness();
}


Population::~Population(void)
{
}

void Population::Print_Ind_Vector( Individual* ind )
{
	//print individual ind
	for( int i = 0; i < 30; i++ )
	{
		cout << ind->vector[i] << " ";
	}
	cout << endl;
}

void Population::Print_Vector()
{
	//print all of the individuals
	for(int i = 0; i < POP_SIZE; i++ )
	{
		cout << "Individual " << i+1 << endl;
		for( int j = 0; j < 30; j++ )
		{
			cout << individual[i]->vector[j] << " ";
		}
		cout << endl << endl;
	}
}

void Population::Print_Fitness()
{
	//print the fitnesses
	Fitness();
	for( int i = 0; i < POP_SIZE; i++ )
	{
		cout << "Individual " << i+1 << " Fitness: " << individual[i]->fitness << endl;
	}
}

void Population::Fitness()
{
	//calculate all the fitness values
	//cout << "Calculating Fitness Values" << endl;

	char fit = func;
	
	// Sphere
	if( fit == 's' )
	{
		double y = 0;
		for( int j = 0; j < POP_SIZE; j++ )
		{
			y = 0;
			for(int i = 0; i < 30; i++)
			{
				y += (individual[j]->vector[i] * individual[j]->vector[i]);
			}
			individual[j]->fitness = y;
		}
	}

	// Rosenbrock
	else if( fit == 'r' )
	{
		double y = 0;
		for( int j = 0; j < POP_SIZE; j++ )
		{
			y = 0;
			for( int i = 0; i < 29; i++ )
			{
				y = y + ( 100 * ( individual[j]->vector[i+1] - (individual[j]->vector[i] * individual[j]->vector[i] ) ) + ((individual[j]->vector[i]-1)*(individual[j]->vector[i]-1)));
			}
			individual[j]->fitness = y;
		}
	}

	// Rastrigin
	else if( fit == 'i' )
	{
		double y = 0;
		for( int j = 0; j < POP_SIZE; j++ )
		{
			y = 300;
			for( int i = 0; i < 29; i++ )
			{
				y = y + ( (individual[j]->vector[i] * individual[j]->vector[i]) - 10*( cos(2*PI*individual[j]->vector[i])));
			}
			individual[j]->fitness = y;
		}
	}

	// Schwefel
	else if( fit == 'c' )
	{
		double y = 0;
		for( int j = 0; j < POP_SIZE; j++ )
		{
			y = 418.9829 * 30;
			for( int i = 0; i < 30; i++ )
			{
				y = y + ( individual[j]->vector[i] * sin( sqrt( abs(individual[j]->vector[i]) ) ) );
			}
			individual[j]->fitness = y;
		}
	}

	// Ackley
	else if( fit == 'a' )
	{
		double x = 0;
		double y = 0;
		double z = 0;
		for( int j = 0; j < POP_SIZE; j++ )
		{
			x = 20 + exp(1.0);
			y = 0;
			z = 0;
			for( int i = 0; i < 30; i++ )
			{
				z = z + ( individual[j]->vector[i] * individual[j]->vector[i] );
				y = y + ( cos( 2* PI * individual[j]->vector[i] ) );
			}
			z = -0.2*sqrt((1/30) * z);
			y = exp( (1/30) * y );
			x = x - 20*exp( z - y );
			individual[j]->fitness = x;
		}
	}

	// Griewangk
	else if( fit == 'g' )
	{
		double x = 0;
		double y = 0;
		double z = 0;
		for( int j = 0; j < POP_SIZE; j++ )
		{
			x = 1;
			y = 0;
			z = 0;
			for( int i = 0; i < 30; i++ )
			{
				z = z + (( individual[j]->vector[i] * individual[j]->vector[i] ) / 4000 );
				y = y * ( cos( individual[j]->vector[i] / sqrt((double)i+1) ) );
			}
			x = 1;
			x = x + z - y;
			individual[j]->fitness = x;
		}
	}

	// Sort
	Sort();
	//cout << "End Calculating Fitness" << endl;
}

void Population::Single_Fitness(Individual* ind)
{
	//single fitness
	char fit = func;
	
	// sphere
	if( fit == 's' )
	{
		double y = 0;
		for( int j = 0; j < 1; j++ )
		{
			y = 0;
			for(int i = 0; i < 30; i++)
			{
				y += (ind->vector[i] * ind->vector[i]);
			}
			ind->fitness = y;
		}
	}

	// Rosenbrock
	else if( fit == 'r' )
	{
		double y = 0;
		for( int j = 0; j < 1; j++ )
		{
			y = 0;
			for( int i = 0; i < 29; i++ )
			{
				y = y + ( 100 * ( ind->vector[i+1] - (ind->vector[i] * ind->vector[i] ) ) + ((ind->vector[i]-1)*(ind->vector[i]-1)));
			}
			ind->fitness = y;
		}
	}

	// Rastrigin
	else if( fit == 'i' )
	{
		double y = 0;
		for( int j = 0; j < 1; j++ )
		{
			y = 300;
			for( int i = 0; i < 29; i++ )
			{
				y = y + ( (ind->vector[i] * ind->vector[i]) - 10*( cos(2*PI*ind->vector[i])));
			}
			ind->fitness = y;
		}
	}

	// Schwefel
	else if( fit == 'c' )
	{
		double y = 0;
		for( int j = 0; j < 1; j++ )
		{
			y = 418.9829 * 30;
			for( int i = 0; i < 29; i++ )
			{
				y = y + ( ind->vector[i] * sin( sqrt( abs(ind->vector[i]) ) ) );
			}
			ind->fitness = y;
		}
	}

	// Ackley
	else if( fit == 'a' )
	{
		double x = 0;
		double y = 0;
		double z = 0;
		for( int j = 0; j < 1; j++ )
		{
			x = 20 + exp(1.0);
			y = 0;
			z = 0;
			for( int i = 0; i < 30; i++ )
			{
				z = z + ( ind->vector[i] * ind->vector[i] );
				y = y + ( cos( 2* PI * ind->vector[i] ) );
			}
			z = -0.2*sqrt((1/30) * z);
			y = exp( (1/30) * y );
			x = x - 20*exp( z - y );
			ind->fitness = x;
		}
	}

	// Griewangk
	else if( fit == 'g' )
	{
		double x = 0;
		double y = 0;
		double z = 0;
		for( int j = 0; j < 1; j++ )
		{
			x = 1;
			y = 0;
			z = 0;
			for( int i = 0; i < 30; i++ )
			{
				z = z + (( ind->vector[i] * ind->vector[i] ) / 4000 );
				y = y * ( cos( ind->vector[i] / sqrt((double)i+1) ) );
			}
			x = x + z - y;
			ind->fitness = x;
		}
	}
}

void Population::Sort()
{
	//sort the individuals
	quicksort(individual, 0, POP_SIZE - 1);
}

void Population::Tournament_Selection()
{
	cout << "Tournament Selection:" << endl;

	double low = 10000;
	int count = 1;
	int repeat = 0;
	while( low > 0.01 || low < -0.01 )
	{
		if( repeat == 2000 || count == 1000 )
		{
			//Print_Fitness();
			break;
		}
		//tournament selection
		int k = 0;
		int p = 0;
		int N = 5;
		Individual* IND[POP_SIZE];

		//IND[0] = individual[0];
		//IND[1] = individual[0];

		while( k < POP_SIZE )
		{
			Individual* ind = new Individual;
			ind->vector = new double[30];
			ind->fitness = 10000;
			for(int j = 0; j < 30; j++ )
			{
				ind->vector[j] = max_range;
			}
	
			for( int i = 0; i < N; i++ )
			{
				int random = ((int)rand()) / (int)RAND_MAX;
				int diff = POP_SIZE - 1 - 0;
				int r = random * diff;
				int index = 0 + r;
	
				if( abs(individual[index]->fitness) < abs(ind->fitness) )
				{
					ind->fitness = individual[index]->fitness;
					for(int j = 0; j < 30; j++ )
					{
						ind->vector[j] = individual[index]->vector[j];
					}
				}
			}
			IND[k] = ind;
			k++;
		}

		// crossover
		Crossover(IND);
	
		// mutate
		//Print_Ind_Vector(IND[0]);
		Mutate(IND);
		//Print_Ind_Vector(IND[0]);


		for( int i = 0; i < POP_SIZE; i++ )
		{
			individual[i]->fitness = IND[i]->fitness;
			
			for( int j = 0; j < 30; j++ )
			{
				individual[i]->vector[j] = IND[i]->vector[j];
			}
			//delete IND[i]->vector;
			IND[i] = NULL;
			delete IND[i];
			IND[i] = NULL;

		}
	
		Fitness();

		double saved = low;
	
		while( p != POP_SIZE )
		{
			if( abs(individual[p]->fitness) < abs(low) )
			{
				low = individual[p]->fitness;
				index = p;
			}
			p++;
		}

		//Print_Ind_Vector(individual[index]);

		if( low == saved )
			repeat++;
		else
			repeat = 0;
		
		//cout << individual[index]->fitness << endl;

		count++;
	}
	
	cout << "Final fitness of best solution after " << count << " rounds: " << low << endl;
}

void Population::Roulette_Selection()
{
}

void Population::Mutate( Individual* ind[] )
{

	for( int i = 0; i < POP_SIZE; i++ )
	{

		for( int j = 0; j < 30; j++ )
		{
			double random = ((double)rand()) / (double)RAND_MAX;
			double diff = 1.0 - 0.0;
			double r = random * diff;
			//double z  = abs(ind[i]->fitness);

			if( r > 0.50 && r < 0.75 )
			{
				ind[i]->vector[j] -= ratio;

				//Single_Fitness(ind[i]);

				//if( abs(ind[i]->fitness) > z )
					//ind[i]->vector[j] += ratio;
			}
			else if( r < 0.50 && r > 0.25 )
			{
				ind[i]->vector[j] += ratio;

				//Single_Fitness(ind[i]);

				//if( abs(ind[i]->fitness) > z )
					//ind[i]->vector[j] -= ratio;
			}

			if( ind[i]->vector[j] > max_range )
				ind[i]->vector[j] = max_range;
			if( ind[i]->vector[j] < min_range )
				ind[i]->vector[j] = min_range;

			Single_Fitness(ind[i]);

		}
	}
}

void Population::Crossover(Individual* ind[])
{
	for( int j = 0; j < POP_SIZE; j+=2 )
	{
		for( int i = 0; i < 30; i++ )
		{
			double random = ((double)rand()) / (double)RAND_MAX;
			double diff = 1.0 - 0.0;
			double r = random * diff;

			if( r > 0.50 )
			{
				double temp = ind[j]->vector[i];
				ind[j]->vector[i] = ind[j+1]->vector[i];
				ind[j+1]->vector[i] = temp;
			}
		}
	}
}