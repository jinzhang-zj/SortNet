// Generate sorting network using Genetic Algorithm
// Created by Jin Zhang
// Copyright (c) 2014 Jin Zhang. All rights reserved.
// Implements the algorithm suggested by W. Daniel Hillis(1990,1992)

#include <stdio.h>	// printf, NULL
#include <stdlib.h>	// rand, strand
#include <time.h>
#include <utility>	// pair
#include <iostream>	// cout
#include <time.h>	// time
#define chromeL 8

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
	public:
	
	SortNet(int,int,int);
	void evolve(int);
	double fitness(int, int**, double*);
	void output();
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
			this->parasites[i][j] = rand() % 50;
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

// Evolve both the input/solution arrays
// generations: maximal number iterations allowed
void SortNet::evolve(int generations){
	for (int generation=0; generation< generations; generation++){
		// 1. evaluate fitness
		// fitMat: store fitness for each solution with respect to each input
		double** fitMat = new double*[this->solSize]; 
	
		// correctly sorted array
		int** sortedArray = new int* [this->testCases];
		for (int i=0; i< this->testCases; i++){
			sortedArray[i] = new int[this->arraySize];
			quickSort(this->parasites[i], this->arraySize, sortedArray[i]);
		}

		// output sorted array
		for (int i=0; i< this->testCases; i++){
			for (int j=0; j< this->arraySize; j++)
				cout << sortedArray[i][j] << " ";
			cout << endl;
		}


		// averaged & normalized fitness for each solution/question
		double solFit[this->solSize];
		double queFit[this->testCases];
		double sum = 0;
		double max = -10;

		//cout << "get initial fitness start" << endl;

		//#pragma omp parallel for
		for (int i=0; i< this->solSize; i++){
			fitMat[i] = new double[this->testCases];
			solFit[i] = this->fitness(i, sortedArray, fitMat[i]);
		}

		//cout << "get initial fitness" << endl;

		// sum, and find the optimal solution with maximal fitness
		for (int i=0; i< this->solSize; i++){
			sum += solFit[i];
			if (solFit[i] > max){
				max = solFit[i];
				this->optimal = i;
			}
		}

		// normalization
		for (int i=0; i< this->solSize; i++){
			solFit[i] /= sum;
		}

		for (int i=0; i< this->solSize; i++){
			cout << solFit[i] << " ";
		}
		cout << endl;

		sum = 0;
		// input array fitness
		for (int i=0; i< this->testCases; i++){
			queFit[i] = 0;
			for (int j=0; j< this->solSize; j++)
				queFit[i] += 1 - fitMat[j][i];
			queFit[i] /= this->solSize;
			sum += queFit[i];
		}

		// normalization
		for (int i=0; i< this->testCases; i++){
			queFit[i] /= sum;
		}

		//cout << "evluation of fitness done" << endl;
		
		// 2. select m out of n solutions with replacement, cross over, generate n new solutions
		// opitional: top k solutions from previous iteration can be kept
		// sols: rank of selected solutions
		vector<int> sols;
		srand(time(NULL));
		int m = 100;
		while(m){
			int idx = rand() % this->solSize;
			double prob = (double)rand() / RAND_MAX;
			if (prob < solFit[idx]){	
				sols.push_back(idx);
				m--;
			}
		}



		// 3. mutation with some probability


		// cleaning up
		delete [] fitMat;
		delete [] sortedArray;
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
	cout << "hosts: " << endl;
	for (int i=0; i<this->solSize*2; i++){
		for (int j=0; j<chromeL; j++){
			cout << "(" << this->hosts[i][j].first << "," << this->hosts[i][j].second << ")" << " ";
		}
		cout << endl;
	}
};

// destructor
SortNet::~SortNet(){
	if (this->parasites != NULL)
		delete [] this->parasites;
	if (this->hosts != NULL)
		delete [] this->hosts;
};

// void sort
void testsort(){
	int a[7] = {1,1,2,2,1,1,2};
	int* sa = new int[7];
	quickSort(a, 7, sa);

	int d[9] = {3,5,4,2,1,7,8,4,9};
	int* sd = new int[9];
	quickSort(d, 9, sd);

	for (int i=0; i<7; i++)
		cout << sa[i] << " ";
	cout << endl;

	for (int i=0; i<9; i++)
		cout << sd[i] << " ";
	cout << endl;

	delete [] sa, sd;
};

int main(int argc, char* argv[])
{
	
	// network size, number of inputs, number of solutions
	SortNet* sn = new SortNet(8,1,2);
	cout << "constructing sorting network done." << endl;
	sn->output();	
	//testSort();
	sn->evolve(1);
	cout << "evolving done." << endl;
	return 0;
}
