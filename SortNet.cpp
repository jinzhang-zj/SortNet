// Generate sorting network using Genetic Algorithm
// Created by Jin Zhang
// Copyright (c) 2014 Jin Zhang. All rights reserved.
// Implements the algorithm suggested by W. Daniel Hillis(1990,1992)

#include <stdio.h>	// printf, NULL
#include <stdlib.h>	// rand, strand
#include <time.h>
#include <utility>	// pair
#include <iostream>	// cout

using namespace std;

class SortNet{

	// parasits: the input arrays
	// hosts: solutions
	int** parasites;
	pair<int, int>** hosts;
	int arraySize;
	int testCases;
	int solSize;
	public:
	
	SortNet(int,int,int);
	void evolve(int);
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
		// the lenght of each chromosome is 50
		this->hosts[i] = new pair<int, int>[50];
		for (int j=0; j< 50; j++){
			int left = rand()%arraySize;
			int right = rand()%arraySize;
			if (left < right)	this->hosts[i][j] = make_pair(left,right);
			else if (left > right)	this->hosts[i][j] = make_pair(right, left);
			else if (left > 0)	this->hosts[i][j] = make_pair(left--,right);
			else			this->hosts[i][j] = make_pair(left, right++);
		}
	}
};


// helper function to get sum of the array
double sumArray(double *array, int size){
	double sum = 0;
	for (int i=0; i< size; i++)
		sum += array[i];
	return sum;
}

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
// steps: maximal number iterations allowed
void SortNet::evolve(int steps){
	for (int step=0; step< steps; step++){
		// 1. evaluate fitness
		double hostFit[this->solSize];
		double parasiteFit[this->testCases];
	
		// correctly sorted array
		int** sortedArray = new int* [this->testCases];
		for (int i=0; i< this->testCases; i++){
			sortedArray[i] = new int[this->arraySize];
			quickSort(this->parasites[i], this->arraySize, sortedArray[i]);
		}

		for (int i=0; i<this->testCases; i++){
			for (int j=0; j<this->arraySize; j++)
				cout << sortedArray[i][j] << " ";
			cout << endl; 
		}

		// 2. select m out of n solutions, cross over, generate n new solutions
		// opitional: top k solutions from previous iteration can be kept
		// 3. mutation with some probability
	}
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
		for (int j=0; j<50; j++){
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
void testSort(){
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
	SortNet* sn = new SortNet(16,2,2);
	cout << "constructing sorting network done." << endl;
	sn->output();	
	//testSort();
	sn->evolve(1);
	return 0;
}
