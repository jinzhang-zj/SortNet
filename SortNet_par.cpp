// Generate sorting network using Genetic Algorithm
// Created by Jin Zhang
// Copyright (c) 2014 Jin Zhang. All rights reserved.
// Implements the algorithm suggested by W. Daniel Hillis(1990,1992)

#include <stdio.h>	// printf, NULL
#include <stdlib.h>	// rand, strand
#include <algorithm>	// sort
#include <time.h>
#include <utility>	// pair
#include <iostream>	// cout
#include <time.h>	// time
#include <omp.h>

#define chromeL 50

using namespace std;

class SortNet{

	// parasits: the input arrays
	// hosts: solutions
	int** parasites;
	pair<int, int>** hosts;
	int arraySize;
	int testCases;
	int solSize;
	int optimal;
	double optFit;
	public:
	
	SortNet(int,int,int);
	void life(int);
	double fitness(int, int**, double*);
	void output();
	void solCrossover(int*, int);
	void inputCrossover(int*, int);
	void optSol();
	~SortNet();
};

// generate array to be sorted
// arraySize: size of array
// testCases: number of input arrays to be sorted
// solSize: number of solutions maintained
SortNet::SortNet(int arraySize, int testCases, int solSize){
	
	this->testCases = testCases;
	this->arraySize = arraySize;
	this->solSize = solSize;

	// construct parasites
	this->parasites = new int*[this->testCases];
	for (int i=0; i< testCases; i++){
		this->parasites[i] = new int[arraySize];
		for (int j=0; j< arraySize; j++){
			this->parasites[i][j] = rand() % 100;
		}
	}

	// construct hosts
	this->hosts = new pair<int, int>*[solSize*2];
	for (int i=0; i< solSize*2; i++){
		// the lenght of each chromosome is chromeL
		this->hosts[i] = new pair<int, int>[chromeL];
		for (int j=0; j< chromeL; j++){
			int left = rand()%arraySize;
			int right = rand()%arraySize;
			if (left < right)	this->hosts[i][j] = make_pair(left,right);
			else if (left > right)	this->hosts[i][j] = make_pair(right, left);
			else if (left > 0)	this->hosts[i][j] = make_pair(left--,right);
			else			this->hosts[i][j] = make_pair(left, right++);
		}
	}
};


// swap a and b
void swap(int* a, int* b){
	int tmp = *a;
	*a = *b;
	*b = tmp;
}

// helper function for quick sort
// result: array to be sorted
void quickSortHelper(int* result, int start, int end){
	if (start < end){
		int pivotVal = result[start];
		int point = start;

		swap(result[start], result[end]);

		for (int i=start; i< end; i++){
			if (result[i] <= pivotVal){
				swap(result[i], result[point]);
				point++;
			}
		}
		swap(result[point], result[end]);

		if (point > 0)		quickSortHelper(result, start, point-1);
		if (point < end)	quickSortHelper(result, point+1, end);
	}
};

// quick sort for int array
// array: original array
// size: size of the array
// result: sorted array
void quickSort(int* array, int size, int* result){
	copy(array, array+size, result);
	int start = 0;
	int end = size - 1;
	quickSortHelper(result, start, end);
};

// comparator for sort fitness and return index
bool compare(const pair<double, int>& pr1, const pair<double, int>& pr2){
	return pr1.first > pr2.first;
}

// Evolve both the input/solution arrays
// generations: maximal number iterations allowed
void SortNet::life(int generations){
	cout << "life cycle begins:" << endl;
	for (int generation=0; generation< generations; generation++){
		//cout << "generation " << generation << endl;
		//this->output();
		// 1. evaluate fitness
		// fitMat: store fitness for each solution with respect to each input
		double** fitMat = new double*[this->solSize]; 
	
		// correctly sorted array
		int** sortedArray = new int* [this->testCases];
		#pragma omp parallel for
		for (int i=0; i< this->testCases; i++){
			sortedArray[i] = new int[this->arraySize];
			quickSort(this->parasites[i], this->arraySize, sortedArray[i]);
		}

		// averaged & normalized fitness for each solution/question
		double solFit[this->solSize];
		pair<double, int> queFit[this->testCases];
		double sum = 0;

		#pragma omp parallel for
		for (int i=0; i< this->solSize; i++){
			fitMat[i] = new double[this->testCases];
			solFit[i] = this->fitness(i, sortedArray, fitMat[i]);
		}

		this->optFit = 0;
		// sum, and find the optimal solution with maximal fitness
		for (int i=0; i< this->solSize; i++){
			sum += solFit[i];
			if (solFit[i] > this->optFit){
				this->optFit = solFit[i];
				this->optimal = i;
			}
		}

		// normalization
		#pragma omp parallel for
		for (int i=0; i< this->solSize; i++){
			solFit[i] /= sum;
		}

		sum = 0;

		// first is fitness, second is index
		// input array fitness
		for (int i=0; i< this->testCases; i++){
			queFit[i] = make_pair(0,i);
			for (int j=0; j< this->solSize; j++)
				queFit[i].first += 1 - fitMat[j][i];
			queFit[i].first /= this->solSize;
			sum += queFit[i].first;
		}

		// normalization
		for (int i=0; i< this->testCases; i++){
			queFit[i].first /= sum;
		}

		//cout << "array fitness: ";
		//for (int i=0; i< this->testCases; i++){
		//	cout << queFit[i] << " ";
		//}
		//cout << endl;

		// 2. select m out of total solutions with replacement, cross over, generate n new solutions
		// opitional: top k solutions from previous iteration can be kept
		// sols: rank of selected solutions
		int m = this->solSize/2;
		int i = 0;
		int sols[m];
		while(m){
			int idx = rand() % this->solSize;
			double prob = (double)rand() / RAND_MAX;
			if (prob < solFit[idx]){	
				sols[i++] = idx;
				m--;
			}
		}
		// solution evolution,  cross over
		this->solCrossover(sols,i);
	
		// select m out of total arrays
		sort(queFit, queFit + this->testCases,  compare);

		//cout << "queFit: ";
		//for (i=0; i< this->testCases; i++){
		//	cout << "(" << queFit[i].first << "," << queFit[i].second << ") ";
		//}
		//cout << endl;

		m = this->testCases/4;
		int ques[m];
		for (int i=0; i<m; i++){
			ques[i] = queFit[i].second;
		}
		// input evolution, best fit survives!
		this->inputCrossover(ques, m);

		//while(n){
		//	int idx = rand() % this->testCases;
		//	double prob = (double) rand() / RAND_MAX;
		//	if (prob < queFit[i]){
		//		ques[i++] = idx;
		//		n--;
		//	}
		//}

		// 3. mutation in hosts
		// mutation frequency: 0.001
		// #pragma omp parallel for
		for (int i=0; i< this->solSize*2; i++){
			for (int j=0; j< chromeL; j++){
				// mutation with probability 0.001
				if ( rand()%1000 < 1){
					int left = rand() % this->arraySize;
					int right = rand() % this->arraySize;
					if (left < right)	this->hosts[i][j] = make_pair(left,right);
					else if (left > right)	this->hosts[i][j] = make_pair(right,left);
					else if (left > 0)	this->hosts[i][j] = make_pair(left--,right);
					else			this->hosts[i][j] = make_pair(left,right++);
				}
			}
		}

		if (! (generation%500)){
			cout << "generation: " << generation << endl;
			this->optSol();
		}
		// cleaning up
		delete [] fitMat;
		delete [] sortedArray;
	}
	cout << "life cycle ends." << endl;
};

// solution cross over
void SortNet::solCrossover(int* bestSol, int size){
	pair<int, int> newHosts[this->solSize*2][chromeL];
	//#pragma omp parallel for
	for (int i=0; i< this->solSize; i++){
		int far = bestSol[rand() % size];
		int mot = bestSol[rand() % size];
		
		// generate gamete and crossover
		// gamete from father
		int breakpoint = rand() % (chromeL - 5) + 2;
		for (int j=0; j< chromeL; j++){
			if (j <= breakpoint)	newHosts[2*i][j] = this->hosts[2*far][j];
			else			newHosts[2*i][j] = this->hosts[2*far+1][j];
		}
		// gamete from mother
		breakpoint = rand() % (chromeL - 5) + 2;
		for (int j=0; j< chromeL; j++){
			if (j <= breakpoint)	newHosts[2*i+1][j] = this->hosts[2*mot][j];
			else			newHosts[2*i+1][j] = this->hosts[2*mot+1][j];
		}
	}
	// update the solution
	// #pragma omp parallel for
	for (int i=0; i< this->solSize*2; i++){
		for (int j=0; j< chromeL; j++)
			this->hosts[i][j] = newHosts[i][j];
	}
};

// input cross over
void SortNet::inputCrossover(int* bestIn, int size ){
	int newParasite[this->testCases][this->arraySize];
	//#pragma omp parallel for
	for (int i=0; i< this->testCases; i++){
		// top size input will enter next generation
		if (i < size){
			for (int j=0; j< this->arraySize; j++)
				newParasite[i][j] = this->parasites[bestIn[i]][j];
		}else{
			for (int j=0; j< this->arraySize; j++)
				newParasite[i][j] = rand() % 100;
		}
		// cross over
		/*
		int far = bestIn[rand() % size];
		int mot = bestIn[rand() % size];

		int breakpoint = rand() % this->arraySize;

		// directly cross over father and mother
		for (int j=0; j< this->arraySize; j++){
			if (j <= breakpoint)	newParasite[i][j] = this->parasites[far][j];
			else			newParasite[i][j] = this->parasites[mot][j];
		}*/
	}
	// update the input
	// #pragma omp parallel for
	for (int i=0; i< this->testCases; i++){
		for (int j=0; j< this->arraySize; j++)
			this->parasites[i][j] = newParasite[i][j];
	}
};


// helper function: compare and swap the values
void swapcomp(pair<int, int>& pr, int* result){
	if(result[pr.first] > result[pr.second])	swap(result[pr.first], result[pr.second]);
};

// Calculate the fitness of each solution with respect to all input array
// i: the ith solution
// answer: correct sorted array
// fitMat: fitness matrix to be filled in
// return averaged fitness for given solution
double SortNet::fitness(int i, int** answer, double * fitMat){
	double avgfit = 0;
	for (int j=0; j< this->testCases; j++){
		pair<int,int>* chrome1 = this->hosts[i*2];
		pair<int,int>* chrome2 = this->hosts[i*2+1];
		int result[this->arraySize];

		copy(this->parasites[j], this->parasites[j]+this->arraySize, result);

		for (int k=0; k< chromeL; k++){
			swapcomp(chrome1[k], result);	
			swapcomp(chrome2[k], result);	
		}

		// swapped results
		//cout << "swapped results:" << endl;
		//for (int k=0; k< this->arraySize; k++){
		//	cout << result[k] << " ";
		//}
		//cout << endl;

		double fit=0;
		for (int k=0; k< this->arraySize; k++){
			if (result[k] == answer[j][k])	fit++;
		}

		fitMat[j] = fit/this->arraySize;
		avgfit += fitMat[j];
	}
	//cout << avgfit << endl;
	return avgfit/this->testCases;
};

// display all the current information
// for debug purpose
void SortNet::output(){
	cout << "displaying sorting network." << endl;
	cout << "parasites: " << endl;
	for (int i=0; i<this->testCases; i++){
		for (int j=0; j<this->arraySize; j++){
			cout << this->parasites[i][j] << " ";
		}
		cout << endl;
	}
	/*cout << "hosts: " << endl;
	for (int i=0; i<this->solSize*2; i++){
		for (int j=0; j<chromeL; j++){
			cout << "(" << this->hosts[i][j].first << "," << this->hosts[i][j].second << ")" << " ";
		}
		cout << endl;
	}*/
};

// print out the solution length, and fitness of the best solution
void SortNet::optSol(){
	pair<int,int> optSol[chromeL*2];
	int l = 0;
	
	for (int i=0; i< chromeL; i++){
		if (this->hosts[2*this->optimal][i] == this->hosts[2*this->optimal + 1][i]){
			optSol[l++] = this->hosts[2*this->optimal][i];
		}
		else{
			optSol[l++] = this->hosts[2*this->optimal][i];
			optSol[l++] = this->hosts[2*this->optimal+1][i];
		}
	}
	cout << "optimal solution has length " << l;
	cout << " and fitness " << this->optFit << endl;
	for (int i=0; i<l; i++){
		cout << "(" << optSol[i].first << "," << optSol[i].second << ") ";
	}
	cout << endl;
}

// destructor
SortNet::~SortNet(){
	if (this->parasites != NULL)
		delete [] this->parasites;
	if (this->hosts != NULL)
		delete [] this->hosts;
};

int main(int argc, char* argv[])
{
	srand(time(NULL));
	if (argc != 5){
		cerr <<  argv[0] << "  [input_num] [solution_num]  [generations] [num_of_threads]" << endl;
		return 0;
	}

	//cout << atoi(argv[1]) << atoi(argv[2]) << atoi(argv[3]) << endl;
	omp_set_num_threads(atoi(argv[4]))
	double start_time = omp_get_wtime();	
	// network size, number of inputs, number of solutions
	SortNet* sn = new SortNet(16,atoi(argv[1]), atoi(argv[2]));
	cout << "constructing sorting network done." << endl;
	//testSort();
	sn->life(atoi(argv[3]));
	sn->optSol();
	cout << "Total time: " << omp_get_wtime() - start_time << endl;
	//sn->output();	
	return 0;
}
