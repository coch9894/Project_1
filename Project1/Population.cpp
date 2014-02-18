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
	if( func == 's' || func == 'g' )
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
	// generate 1000 random individuals
	srand(time(NULL));
	for( int i = 0; i < 1000; i++ )
	{
		individual[i] = new Individual;
		individual[i]->vector = new double[30];
		for( int j = 0; j < 30; j++ )
		{
			double random = ((double)rand()) / (double)RAND_MAX;
			double diff = max_range - min_range;
			double r = random * diff;
			individual[i]->vector[j] = min_range + r;
		}
		//individual[i]->vector should now = foo;
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
	for(int i = 0; i < 1000; i++ )
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
	for( int i = 0; i < 1000; i++ )
	{
		cout << "Individual " << i+1 << " Fitness: " << individual[i]->fitness << endl;
	}
}

void Population::Fitness()
{
	//calculate all the fitness values
	//cout << "Calculating Fitness Values" << endl;

	char fit = func;
	
	// sphere
	if( fit == 's' )
	{
		double y = 0;
		for( int j = 0; j < 1000; j++ )
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
	if( fit == 'r' )
	{
		double y = 0;
		for( int j = 0; j < 1000; j++ )
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
	if( fit == 'g' )
	{
		double y = 0;
		for( int j = 0; j < 1000; j++ )
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
	if( fit == 'c' )
	{
		double y = 0;
		for( int j = 0; j < 1000; j++ )
		{
			y = 418.9829 * 30;
			for( int i = 0; i < 30; i++ )
			{
				y = y + ( individual[j]->vector[i] * sin( sqrt( individual[j]->vector[i] ) ) );
			}
			individual[j]->fitness = y;
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
	if( fit == 'r' )
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
	if( fit == 'g' )
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
	if( fit == 'c' )
	{
		double y = 0;
		for( int j = 0; j < 1; j++ )
		{
			y = 418.9829 * 30;
			for( int i = 0; i < 29; i++ )
			{
				y = y + ( ind->vector[i] * sin( sqrt( ind->vector[i] ) ) );
			}
			ind->fitness = y;
		}
	}
}

void Population::Sort()
{
	//sort the individuals
	quicksort(individual, 0, 999);
}

void Population::Tournament_Selection()
{
	cout << "Tournament Selection:" << endl;

	double low = 10000;
	int count = 1;
	int repeat = 0;
	while( low > 0.01 || low < -0.01 )
	{
		if( repeat == 1000000 || count == 1000 )
		{
			//do nothing
			Print_Fitness();
			break;
		}
		//tournament selection
		int k = 0;
		int p = 0;
		int N = 5;
		Individual* IND[1000];
	
		while( k < 1000 )
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
				int diff = 999 - 0;
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
	
		// mutate
		Mutate(IND);
	
		for( int i = 0; i < 1000; i++ )
		{
			individual[i]->fitness = IND[i]->fitness;
			
			for( int j = 0; j < 30; j++ )
			{
				individual[i]->vector[j] = IND[i]->vector[j];
			}
			//delete IND[i]->vector;
			delete IND[i];
		}

		Crossover();
	
		Fitness();

		double saved = low;
	
		while( p != 1000 )
		{
			if( abs(individual[p]->fitness) < abs(low) )
			{
				low = individual[p]->fitness;
				index = p;
			}
			p++;
		}

		if( low == saved )
			repeat++;
		else
			repeat = 0;
		
		//cout << low << " ";

		count++;
	}
	
	cout << "Final fitness of best solution after " << count << " rounds: " << low << endl;
}

void Population::Roulette_Selection()
{
}

void Population::Mutate( Individual* ind[] )
{

	for( int i = 0; i < 1000; i++ )
	{

		for( int j = 0; j < 30; j++ )
		{
			double random = ((double)rand()) / (double)RAND_MAX;
			double diff = 1.0 - 0.0;
			double r = random * diff;
			double z  = abs(ind[i]->fitness);

			if( r > 0.50 )
			{
				ind[i]->vector[j] -= 0.01;

				Single_Fitness(ind[i]);

				if( abs(ind[i]->fitness) > z )
					ind[i]->vector[j] += 0.01;

				if( ind[i]->vector[j] > max_range )
					ind[i]->vector[j] = max_range;
				if( ind[i]->vector[j] < min_range )
					ind[i]->vector[j] = min_range;
			}
			else
			{
				ind[i]->vector[j] += 0.01;

				Single_Fitness(ind[i]);

				if( abs(ind[i]->fitness) > z )
					ind[i]->vector[j] -= 0.01;

				if( ind[i]->vector[j] > max_range )
					ind[i]->vector[j] = max_range;
				if( ind[i]->vector[j] < min_range )
					ind[i]->vector[j] = min_range;
			}

		}
	}
}

void Population::Crossover()
{
	double temp;
	for( int j = 0; j < 1000; j+=2 )
	{
		for(int i = 0; i < 15; i++ )
		{
			temp = individual[j]->vector[i];
			individual[j]->vector[i] = individual[j+1]->vector[i+14];
			individual[j+1]->vector[i+14] = temp;
		}
	}
}