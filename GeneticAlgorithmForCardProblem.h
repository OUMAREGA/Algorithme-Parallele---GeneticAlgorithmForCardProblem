#ifndef GENETICALGORITHMFORCARDPROBLEM_H
#define GENETICALGORITHMFORCARDPROBLEM_H

#include <iostream >
#include <conio.h>
#include <Windows.h>
#include <cmath>
#include <mpi.h>
#include <cstdlib >
#include <algorithm>
#include <ctime>
#include <cstdlib>

class geneticAlgorithmForCardProblem
{
public:
	geneticAlgorithmForCardProblem(int maxiteration, int popsize, int dim, double pcross, double pmut);
	~geneticAlgorithmForCardProblem();
	void run();
	double random();
	int maxIterations() const;
	int d_popSize() const;
	int dimension() const;
	double crossoverProbability() const;
	double mutationProbability const;
	double SUMTARG();
	double PRODTARG();
	int bestIndice();
	int worstIndice();
	void initialPopulation();
	void display(int interaction, int individual);
	int population(int i, int j) const;
private:
	int d_maxIterations = 1000;
	int d_popSize = 30;
	int d_dimension = 10;
	double d_crossoverProbability = 0.5;
	double d_mutationProbability = 0.1;
	double d_SUMTARG = 36;
	double d_PRODTARG = 360;
	int ** d_population;
	double rnd;
	double evaluate(int n);
	
};
#endif // !GENETICALGORITHMFORCARDPROBLEM_H

