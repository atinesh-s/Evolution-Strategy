#include <stdio.h>
#include "../include/ESSolver.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <iostream>
#include <iomanip>

using namespace std;

int maxFE;
int minVal;
int maxVal;
int selStrategy;
int recOpO;
int recOpS;
int mutOp;
double epsilon;

void es_thread(int ne, int function_eval, double& result, int, int, int, int, int);
double standard_deviation(vector<double> data);
double my_avg(vector<double> x);


int main(int argc, char *argv[]){
	
	// ./main #function #dimension #runs
    
    time_t timer = time(0); 

    ofstream times ("results/time.txt", ofstream::app);
    ofstream best ("results/best.txt", ofstream::app);
    ofstream measures ("results/measures.txt", ofstream::app);

	// ----- Parameters ----- //
	
	// selStrategy = 0, Selection Strategy Type 0 - Comma Selection
    // selStrategy = 1, Selection Strategy Type 1 - Plus Selection
	selStrategy = 1;
	
	// Object Variable
	// recOpO = 0, Recombination Type 0 - Discrete recombination/Dominant recombination
    // recOpO = 1, Recombination Type 1 - Intermediate recombination
	recOpO = 1;
	
	// Strategy Parameter
	// recOpS = 0, Recombination Type 0 - Discrete recombination/Dominant recombination
    // recOpS = 1, Recombination Type 1 - Intermediate recombination
	recOpS = 1;
	
	// mutOp = 0, Mutation Type 0 - Uncorrelated Mutation with 1 Step Size
    // mutOp = 1, Mutation Type 1 - Uncorrelated Mutation with n Step Size
	mutOp = 1;
	epsilon = 0.001; // Minimum allowed value for strategy parameter

	int function_eval = atoi(argv[1]);
	int nVar 		  = atoi(argv[2]);
	int NE 		      = atoi(argv[3]);
	
	int mu   = 30;
    int lbd  = 200;
    int rhoO = 2;
    int rhoS = mu;
	maxFE    = 5000*nVar;	// default 5000*D
	
	//ackley
	if( function_eval == 1){
		minVal = -32; maxVal = 32;
	}
	//rastrigin
	else if( function_eval == 2){
		minVal = -5; maxVal = 5;
	}    
	//rosenbrock
	else if( function_eval == 3){
	    minVal = -10; maxVal = 10;
	}    
	//schwefel
	else if( function_eval == 4){
	    minVal = -100; maxVal = 100;
	}    
	//sphere
	else if( function_eval == 5){
		minVal = -100; maxVal = 100;
	}

	vector<double> best_val(NE);

	//vector of threads
    vector<thread> ths;

	//NE Experiments
	for (int i = 0; i < NE; ++i){
		ths.emplace_back(es_thread, i+1, function_eval, ref(best_val[i]), nVar, mu, lbd, rhoO, rhoS);
		// es_thread(function_eval, best_val[i]);
		// printf("%f\n", best_val[i]);
	}

	//Joining threads
    for ( auto t = ths.begin(); t != ths.end(); ++t){
    	t->join();// ths[i].join();
    }

	best << "es_" << function_eval << "_" << nVar << "= [" ;
    // Print best_val vector
	for (int i = 0; i < NE; ++i){
		best << setw(10) << best_val[i];
	}
	best << "];"<< endl;

    measures << my_avg(best_val) << " ± " << standard_deviation(best_val) << endl;

	time_t timer2 = time(0);
	times << "Total time: " <<  difftime(timer2, timer) << endl;

	return 0;
}


void es_thread(int ne, int function_eval, double& result, int nVar, int mu, int lbd, int rhoO, int rhoS){
	// ----- Running Evolution Strategy ----- //
	ESSolver es(ne, nVar, mu, lbd, rhoO, rhoS);
	es.function_eval = function_eval;
	es.Init();
	es.Optimize();
	result = es.best_val;
}


double standard_deviation(vector<double> data){

    double mean=0.0, sum_deviation=0.0;
    
    int i;
    for(i=0; i < data.size(); ++i){
        mean += data[i];
    }
    mean=mean/data.size();
    
    for(i=0; i < data.size(); ++i){
        sum_deviation += (data[i]-mean)*(data[i]-mean);
    }    
    
    return sqrtf(sum_deviation/data.size());
}    


double my_avg(vector<double> x ){
    double avg_best = 0;
    //Getting the average
    for( int i = 0 ; i < x.size() ; ++i ){
        avg_best += x[i];
    }
    avg_best = avg_best / x.size();
    return avg_best;
}