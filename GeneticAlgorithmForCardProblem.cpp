#include "GeneticAlgorithmForCardProblem.h"

using namespace std;

double geneticAlgorithmForCardProblem::random()
{
	time_t t;
	srand((unsigned)time(&t));
	rnd = rand() % d_popSize;
	return rnd;
}
geneticAlgorithmForCardProblem::geneticAlgorithmForCardProblem(int maxiteration, int popsize, int dim, double pcross, double pmut) : d_maxIterations{ maxiteration }, d_popSize{ popsize }, d_dimension{ dim }, d_crossoverProbability{ pcross }, d_mutationProbability{ pmut }
{
	d_population = new int *[popsize];
	for (int i = 0; i < popsize; ++i)
	{
		d_population[i] = new int[dim];
	}
}
geneticAlgorithmForCardProblem::~geneticAlgorithmForCardProblem()
{
	for (int i = 0; i < d_popSize; i++)
		delete[] d_population[i];
	delete[] d_population;
}
int geneticAlgorithmForCardProblem::maxIterations() const
{
	return d_maxIterations;
}
/*int geneticAlgorithmForCardProblem::evaluate(int a, int b)
{
	if (evaluate(a) > evaluate(b))
		return b;
	else
		return a;
}*/
void geneticAlgorithmForCardProblem::initialPopulation()
{
	//for entire population
	for (int i = 0; i < d_popSize; i++)
	{
		//for all genes
		for (int j = 0; j < d_dimension; j++)
		{
			//randomly create gene values
			if (random() < 0.5)
				d_population[i][j] = 0;
			else
				d_population[i][j] = 1;
		}
	}
}
void geneticAlgorithmForCardProblem::display(int iteration, int individual) //display(int tab)
{
	cout << "\r\n==============================\r\n" << endl;
	int somme = 0;
	int produit = 1;
		cout << "After " << iteration << " Iterations, Solution sum pile (should be 36) cards are : " << endl;
	for (int i = 0; i < d_dimension; i++) {
		if (d_population[individual][i] == 0) {
			cout << i + 1 << endl;
			somme += i + 1;
		}
	}
	cout << "La somme vaut " << somme;
	cout << "\r\nAnd Product pile (should be 360)  cards are : " << endl;
	for (int i = 0; i < d_dimension; i++) {
		if (d_population[individual][i] == 1) {
			cout << i + 1 << endl;
			produit *= i + 1;
		}
	}
	cout << "Le produit vaut " << produit;
}
double geneticAlgorithmForCardProblem::evaluate(int n)
{
	//initialise field values
	int sum = 0, prod = 1;
	double scaled_sum_error, scaled_prod_error, combined_error;
	//loop though all genes for this population member
	for (int i = 0; i < d_dimension; i++)
	{
		//if the gene value is 0, then put it in the sum (pile 0), and calculate sum
		if (d_population[n][i] == 0)
			sum += (1 + i);
		//if the gene value is 1, then put it in the product (pile 1), and calculate sum
		else
			prod *= (1 + i);
	}
	//work out how food this population member is, based on an overall error
	//for the problem domain
	//NOTE : The fitness function will change for every problem domain.
	scaled_sum_error = (sum - d_SUMTARG) / d_SUMTARG;
	scaled_prod_error = (prod - d_PRODTARG) / d_PRODTARG;
	combined_error = abs(scaled_sum_error) + abs(scaled_prod_error);

	return combined_error;
}
int geneticAlgorithmForCardProblem::bestIndice()
{
	double fitness;
	double best = evaluate(0);
	int indice = 0;


	for (int i = 1; i < d_dimension; i++)
	{
		//if the gene value is 0, then put it in the sum (pile 0), and calculate sum

		fitness = evaluate(i);
		if (fitness < best)
		{
			best = fitness;
			indice = i;
		}
		return indice;
	}
}
int geneticAlgorithmForCardProblem::worstIndice()
{
	double fitness;
	double worst = evaluate(0);
	int indice = 0;


	for (int i = 1; i < d_dimension; i++)
	{
		//if the gene value is 0, then put it in the sum (pile 0), and calculate sum

		fitness = evaluate(i);
		if (fitness > worst)
		{
			worst = fitness;
			indice = i;
		}
		return indice;
	}
}
int geneticAlgorithmForCardProblem::dimension() const
{
	return d_dimension;
}
int geneticAlgorithmForCardProblem::population(int i, int j) const
{
	return d_population[i][j];
}
//Runs the Microbial GA to solve the problem domain
//Where the problem domain is specified as follows

//You have 10 cards numbered 1 to 10.
//You have to divide them into 2 piles so that:

//The sum of the first pile is as close as possible to 36
//And the product of all in second pile is as close as poss to 360
void geneticAlgorithmForCardProblem::run()
{

	//declare pop member a,b, winner and loser
	int a, b, Winner, Loser;
	//initialize the population (randomly)
	initialPopulation();
	int BestTab[10];
	int nextProc, prevProc;
	MPI_Status staus;
	//start a tournament
	for (int iter = 0; iter < d_maxIterations; iter++)
	{
		cout << "iter : " << iter << endl;
		//pull 2 population members at random
		//a = rand() % PopSize ;
		double a = random();
		double b = random();
		//have a fight, see who has best genes
		if (evaluate(a) < evaluate(b))
		{
			Winner = a;
			Loser = b;
		}
		else
		{
			Winner = b;
			Loser = a;
		}

		//Possibly do some gene jiggling, on all genes of loser
		//again depends on randomness (simulating the natural selection
		//process of evolutionary selection)
		for (int i = 0; i < d_dimension; i++)
		{
			//Crossover operator
			if (random() < d_crossoverProbability)
				d_population[Loser][i] = d_population[Winner][i];

			//Mutation operator
			if (random() < d_mutationProbability)
				d_population[Loser][i] = 1 - d_population[Loser][i];
			
			if(iter%10==0)//(iter == 10 || iter == 20 || iter == 30 || iter == 40 || iter == 50 || iter == 60 || iter == 70 || iter == 80 || iter == 90)
			{


			int bestindice = bestIndice();
			for (int i = 0; i < d_dimension; i++)
			BestTab[i] = d_population[bestindice][i];
			MPI_Send(&BestTab, d_dimension, MPI_INT, nextProc, 0, MPI_COMM_WORLD);
			MPI_Recv(&BestTab, d_dimension, MPI_INT, prevProc, 0, MPI_COMM_WORLD, &status);
			int indicew = worstIndice();
			for (int j = 0; j < d_dimension; j++)
			d_population[indicew][j] = BestTab[j];
			}
			
			//then test to see if the new population member is a winner
			if (evaluate(Loser) == 0.0)
			{
				display(iter, Loser);
				iter = d_maxIterations;
				int c = a;
				for (int i = 0; i < d_dimension; i++)
					BestTab[i] = d_population[c][i];
				MPI_Send(BestTab, d_dimension, MPI_INT, nextProc, 0, MPI_COMM_WORLD);
				break;
			}
		}
	}

}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int rank, tag = 0, valeurAenvoye, valeurRecu, nbElement, nbProcessus, destination, source, currentProc, prevProc, nextProc,indiceBest,indiceWorst;
	MPI_Status status;
	geneticAlgorithmForCardProblem GA{ 1000,30, 10, 0.5, 0.1 };
	const int a = GA.dimension();
	int tableauBest[30];
	int tableauRecu[30];
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nbProcessus);
	prevProc = (nbProcessus + rank - 1) % nbProcessus;
	nextProc = (rank + 1) % nbProcessus;
	GA.initialPopulation();
	if (rank == 0)
	{
		
		int best = GA.bestIndice();
		for (int i = 0; i < GA.dimension(); ++i)
		{
			tableauBest[best] = GA.population[best][i];
			MPI_Send(&tableauBest, GA.dimension(), MPI_INT, nextProc, tag, MPI_COMM_WORLD);
			MPI_Recv(&tableauRecu, GA.dimension(), MPI_INT, prevProc, tag, MPI_COMM_WORLD, &status);
			cout << "==================================================================" << endl;
			cout << "moi processus  : " << rank << " j ai recu la  solution " <<tableauRecu[i] << " du processus : " << prevProc<< endl;
			GA.display(1000,30);
		}
	}
	else
	{
		MPI_Recv(&tableauRecu, GA.dimension(), MPI_INT, prevProc, tag, MPI_COMM_WORLD, &status);
		indiceWorst = GA.worstIndice();
		for (int j = 0; j < GA.dimension(); j++)
			GA.population[indiceWorst][j] = tableauRecu[j];
		int bestindice = GA.bestIndice();
		for (int i = 0; i < GA.dimension(); i++)
			tableauBest[i] = GA.population[bestindice][i];
		MPI_Send(&tableauBest, GA.dimension(), MPI_INT, nextProc, tag, MPI_COMM_WORLD);
	}
	system("pause");
	MPI_Finalize();
}