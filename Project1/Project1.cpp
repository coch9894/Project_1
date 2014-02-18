// Project1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


int _tmain(int argc, _TCHAR* argv[])
{
	char input = 'r';
	
	while( input == 'r' || input == 'p' || input == 'i' )
	{
		// options
		// s = sphere
		// c = schwefel
		// r = rosenbrock
		// g = rastrigin
		Population p1('g');
		//p1.Print_Fitness();
		
		if( input != 'p' && input != 'i' )
		{
			p1.Tournament_Selection();
		}

		cout << "r will re-run, p to print fitnesses, i to print best individual, and other characters will quit:";
		cin >> input;
		if( input == 'p' )
		{
			p1.Print_Fitness();
		}
		if( input == 'i' )
		{
			p1.Print_Ind_Vector(p1.individual[p1.index]);
		}
	}

	return 0;
}

