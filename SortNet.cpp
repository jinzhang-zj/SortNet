// Generate sorting network using Genetic Algorithm
// Created by Jin Zhang
// Copyright (c) 2014 Jin Zhang. All rights reserved.
//
// Implements the algorithm suggested by W. Daniel Hillis(1990,1992)

#include <stdio.h>
#include <stdlib.h>
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
	void output();
	~SortNet();
};

// generate array to be sorted
// arraySize: size of array
// testCases: number of input arrays to be sorted
// solSize: number of solutions maintained
SortNet::SortNet(int arraySize, int testCases, int solSize){
	cout << "constructing sorting network" << endl;
	
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


// display all the current information
// for debug purpose
void SortNet::output(){
	cout << "displaying sorting network" << endl;

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
}

int main(int argc, char* argv[])
{
	SortNet* sn = new SortNet(16,15,100);
	sn->output();	

	return 0;
}
