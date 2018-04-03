#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <queue>
#include <fstream>
#include <string>
#include "../include/ESSolver.h"
#include <random>

using namespace std;

const double PI = atan(1.0)*4;
// default_random_engine rng( random_device{}() );


ESSolver::ESSolver(int _ne, int _nVar, int _mu, int _lbd, int _rhoO, int _rhoS):ne(_ne), nVar(_nVar), mu(_mu), lbd(_lbd), rhoO(_rhoO), rhoS(_rhoS)
{	
    muPop.resize(mu);
    muSigma.resize(mu);
    muSigma2.resize(mu);
    
    for (int i = 0; i < mu; i++){
        muPop[i].resize(nVar);
        muSigma2[i].resize(nVar);
    }
    
    lambdaPop.resize(lbd);
    lambdaSigma.resize(lbd);
    lambdaSigma2.resize(lbd);
    
    for (int i = 0; i < lbd; i++){
        lambdaPop[i].resize(nVar);
        lambdaSigma2[i].resize(nVar);
    }
    
    mulambdaPop.resize(mu+lbd);
    mulambdaSigma.resize(mu+lbd);
    mulambdaSigma2.resize(mu+lbd);
    
    for (int i = 0; i < mu+lbd; i++){
        mulambdaPop[i].resize(nVar);
        mulambdaSigma2[i].resize(nVar);
    }
    
    rhoPop.resize(rhoO);
    rhoSigma.resize(rhoS);
    rhoSigma2.resize(rhoS);
    
    for (int i = 0; i < rhoO; i++){
        rhoPop[i].resize(nVar);
    }
    
    for (int i = 0; i < rhoS; i++){
        rhoSigma2[i].resize(nVar);
    }
    
	return;
}


ESSolver::~ESSolver(void)
{
    return;
}


void ESSolver::Init()
{
	// Initializing object variable

	for (int i = 0; i < mu; i++)
	{
		for (int j = 0; j < nVar; j++)
		{
			muPop[i][j] = rand_r(minVal, maxVal);
		}
	}
	
	// Initializing strategy parameter
	if (mutOp == 0)	// Mutation Type 0
	{
	    for (int i = 0; i < mu; i++)
	    {
            muSigma[i] = rand01_r();
	    }
	}
	else			// Mutation Type 1
	{
		for (int i = 0; i < mu; i++)
	    {
		    for (int j = 0; j < nVar; j++)
		    {
			    muSigma2[i][j] = rand01_r();
		    }
	    }
    }
}


void ESSolver::Optimize()
{
	ofstream file;
	string filename;	
	// function_dimension_NE
	filename = "plot/" + to_string(function_eval) + "_" + to_string(nVar) + "_" + to_string(ne) + ".txt";
	file.open(filename);
	
	int usedFE = 0;
	int fe     = 0; 		// Fitness utilized during initialization
	while (fe < maxFE)
	{
		usedFE = 0;
        // recombination
        Recombination(recOpO, recOpS, rhoO, rhoS, mutOp);
        
        // mutation
        Mutation(mutOp);
        
        // selection
		Selection(selStrategy, mutOp);
        
        //cout << "\nFitness Evaluation : " << fe << ", Best Value : " << getgBestVal();
        
        // Writing data to File
		file << fe << " " << getgBestVal() << "\n"; 
        
        if (selStrategy == 0)
        {
			usedFE = usedFE + lbd;
		}
		else
		{
			usedFE = usedFE + (mu+lbd);
		}
		
        fe = fe + usedFE;
	}
	
	cout << endl;
	
	file.close();
	
	best_val = getgBestVal();
	
	cout << "Global Best Value : " << getgBestVal() << endl;
}


void ESSolver::sampleRhoPop(int _rhoO)
{
	// Randomly sample rho individuals from muPop
	for (int i = 0; i < _rhoO; i++)
	{   
		for (int j = 0; j < nVar; j++)
		{
            rhoPop[i][j] = muPop[rand_r(1, mu) - 1][j];
		}
	}
}


void ESSolver::sampleRhoSigma(int _rhoS, int _mutOp)
{
	// Randomly sample rho individuals from muSigma
	if (_mutOp == 0) // Mutation Type 0
	{
		for (int i = 0; i < _rhoS; i++)
	    {   
            rhoSigma[i] = muSigma[rand_r(1, mu) - 1];
	    }
    }
    else			 // Mutation Type 1
    {
    	for (int i = 0; i < _rhoS; i++)
	    {   
		    for (int j = 0; j < nVar; j++)
		    {
                rhoSigma2[i][j] = muSigma2[rand_r(1, mu) - 1][j];
		    }
	    }
    }
}


double ESSolver::avg(double *arr, int size, int isSigma)
{
	double sum = 0;
	double avg;
	
	if (isSigma == 0)
	{
		for (int i = 0; i < size; i++)
		{
			sum = sum + arr[i];
		}
		avg = sum/size;
	}
    else
    {
		for (int i = 0; i < size; i++)
		{
			sum = sum + arr[i]*arr[i];
		}
		avg = sqrt(sum/size);
	}
    
    delete []arr;
    return avg;
}


void ESSolver::Recombination(int _recOpO, int _recOpS, int _rhoO, int _rhoS, int _mutOp)
{
	// Recombination of Object Variables
	if (_recOpO == 0) // Discrete/Dominant Recombination
    {
    	for (int i = 0; i < lbd; i++)
	     {
		     sampleRhoPop(_rhoO);
		
		     // Recombination of Object Variables
		     for (int j = 0; j < nVar; j++)
		     {
			     lambdaPop[i][j] = rhoPop[ rand_r(1, _rhoO) - 1 ][j];
		     }
	     }
    }
    else              // Intermediate Recombination
    {
    	for (int i = 0; i < lbd; i++)
	     {
		     sampleRhoPop(_rhoO);
		       
		     // Recombination of Object Variables
		     for (int j = 0; j < nVar; j++)
		     {
			     lambdaPop[i][j] = avg( copyDim(j, &rhoPop, _rhoO), _rhoO, 0 );
		     }
		}
    }
	
    // Recombination of strategy parameters
	if (_recOpS == 0) // Discrete/Dominant Recombination
    {
    	for (int i = 0; i < lbd; i++)
	     {
		     sampleRhoSigma(_rhoS, _mutOp);
		     
		     // Recombination of strategy parameters
		     if (_mutOp == 0) // Mutation Type 0
		     {
		         lambdaSigma[i] = rhoSigma[ rand_r(1, _rhoS) - 1];
		     }
		     else			  // Mutation Type 1
		     {
			     for (int j = 0; j < nVar; j++)
		         {
			         lambdaSigma2[i][j] = rhoSigma2[ rand_r(1, _rhoS) - 1 ][j];
		         }
			 }
	     }
    }
    else              // Intermediate Recombination
    {
    	for (int i = 0; i < lbd; i++)
	     {
		     sampleRhoSigma(_rhoS, _mutOp);
		    
		     // Recombination of strategy parameters
		     if (_mutOp == 0) // Mutation Type 0
		     { 
				 double *trialVec = new double[rhoS];

				 for (int j = 0; j < rhoS; j++)
				 {
					 trialVec[j] = rhoSigma[j];
				 }
				 
		         lambdaSigma[i] = avg(trialVec, _rhoS, 1);
		     }
		     else			  // Mutation Type 1
		     {
			     for (int j = 0; j < nVar; j++)
		         {
			         lambdaSigma2[i][j] = avg( copyDim(j, &rhoSigma2, _rhoS), _rhoS, 1 );
		         }
			 }
	     }
    }
}


void ESSolver::Mutation(int _mutOp)
{
	double taup  = 1/(sqrt(2*nVar));
    double tau   = 1/(sqrt(2*sqrt(nVar)));
	
	if (_mutOp == 0) // Mutation Type 0 - Uncorrelated Mutation with 1 Step Size
    {
		for (int i = 0; i < lbd; i++)
	     {
		     // Mutation of strategy parameters
             lambdaSigma[i] = lambdaSigma[i] * exp(tau * rand01_r());
             
             if (lambdaSigma[i] < epsilon) {
                lambdaSigma[i] = epsilon;
             }
		     
		     // Mutation of Object Variables
		     for (int j = 0; j < nVar; j++) {
				 
			     lambdaPop[i][j] = lambdaPop[i][j] + lambdaSigma[i] * rand01_r();
			     
			     // Checking Boundary
			     if (lambdaPop[i][j] < minVal) {
				     lambdaPop[i][j] = minVal;
				 }
				 if (lambdaPop[i][j] > maxVal) {
				     lambdaPop[i][j] = maxVal;
				 }
		     }
	     }
	}
	else 			 // Mutation Type 1 - Uncorrelated Mutation with n Step Size
	{
		 for (int i = 0; i < lbd; i++)
	     {
		     double randF = rand01_r();
			 
		     // Mutation of Object Variables
		     for (int j = 0; j < nVar; j++) {
				 
			     lambdaSigma2[i][j]  = lambdaSigma2[i][j] * exp(taup * randF + tau * rand01_r() );
				 
			     if (lambdaSigma2[i][j] < epsilon) {
                    lambdaSigma2[i][j] = epsilon;
                 }
			
			     lambdaPop[i][j] = lambdaPop[i][j] + lambdaSigma2[i][j] * rand01_r();
			
			     // Checking Boundary
			     if (lambdaPop[i][j] < minVal) {
				     lambdaPop[i][j] = minVal;
				 }
				 if (lambdaPop[i][j] > maxVal) {
				     lambdaPop[i][j] = maxVal;
				 }
		     }
	     }
	}
}


void ESSolver::Selection(int _selStrategy, int _mutOp)
{   
    if (_selStrategy == 0) // Selection Strategy Type 0 - Comma Selection
    {
		double *val = new double[lbd];
		int *idx; // contains sorted index

		for (int i = 0; i < lbd; i++)
		{
			val[i] = objFun( copyVec(i, &lambdaPop) );
		}
		
		idx = sortArray(val, lbd);
		    
		updatePopSigma(idx, _selStrategy, _mutOp);
		
		// Updating Best Value
		gBestVal = val[ idx[0] ];
		
		delete []val;
		delete idx;
	}
	else 				// Selection Strategy Type 1 - Plus Selection
	{
		double *val = new double[mu+lbd];
		int *idx; // contains sorted index
		
		mergeMuLbd(val, _mutOp);
		
		idx = sortArray(val, mu+lbd);
		
		updatePopSigma(idx, _selStrategy, _mutOp);
		
		// Updating Best Value
		gBestVal = val[ idx[0] ];
		
		/*
		double * trial = copyVec(idx[0], &mulambdaPop);
		
		for (int i = 0; i < nVar; i++){
			cout << trial[i] << "   ";
		}
		cout << "\n\n";*/
		
		delete []val;
		delete idx;
	}
}


double ESSolver::getgBestVal()
{
	return gBestVal;
}


double ESSolver::objFun (const double *_pos)
{
	return evaluation(_pos, function_eval); //Function Evaluation Selector
}


double ESSolver::evaluation(const double *_pos, int function_eval)
{
	double F, bias, shifted;
	if( function_eval == 1 ){ // Ackley 0
        bias = 0.0 ;    
        double a = 20, b = 0.2;
        double exp_1 = 0, exp_2 = 0;
        for(int j = 0 ; j < nVar ; ++j){
            exp_1 += ( _pos[j] - bias ) * ( _pos[j] - bias );
            exp_2 += cos( 2.0*PI*( _pos[j] - bias ) );
        }    
        F = -a * exp(-b*sqrtf(1.0/nVar * exp_1)) - exp(1.0/nVar*exp_2) + a + exp(1);
        // printf("%f\n", F);
    }    
    else if( function_eval == 2 ){ // Rastrigin 0
		bias = ( shifted == 1 ) ? -5.0 : 0.0 ;
		double sum = 0;
		double exp_1 = 0, exp_2 = 0;
		for(int j = 0 ; j < nVar ; ++j){
			exp_1 = ( _pos[j] - bias )*( _pos[j] - bias );
			exp_2 = cos( 2.0*PI*( _pos[j] - bias ) );
			sum += exp_1 - 10*exp_2;
		}
		F = 10*nVar + sum;
	}
	else if( function_eval == 3 ){ // Rosenbrock 1
		bias = ( shifted == 1 ) ? -8.0 : 0.0 ;
		double sum = 0;
		double exp_1 = 0, exp_2 = 0;
		for(int j = 0 ; j < nVar - 1; ++j){
			exp_1 = ( _pos[j+1] - bias ) - powf(( _pos[j] - bias ),2);
			exp_2 = powf(( _pos[j] - bias )-1.0,2);
			sum += 100*powf(exp_1,2) + exp_2;
			// cout << sum << endl;
		}
		F = sum;
	}
	else if( function_eval == 4 ){ // Schwefel 420.9687
		bias = ( shifted == 1 ) ? -50.0 : 0.0 ;
		double sum = 0;
		double exp_1 = 0;
			for(int j = 0 ; j < nVar ; ++j){
			exp_1 = _pos[j] * sin(sqrtf(abs(_pos[j])));
			sum += exp_1;
		}
		F = 418.9829*nVar - sum;
	}
	else if( function_eval == 5 ){ // Sphere 0
		bias = ( shifted == 1 ) ? -50.0 : 0.0 ;
		double sum = 0;
		double exp_1 = 0;
		for(int j = 0 ; j < nVar ; ++j){
			exp_1 = powf(( _pos[j] - bias ),2);
			sum += exp_1;
		}
		F = sum ;
	}
	
	delete []_pos;
	return F;
}


int ESSolver::rand_r(int minVal, int maxVal)
{
    default_random_engine rng( std::random_device{}() );
    // mt19937 rng(time(0));
    uniform_int_distribution<int> dist( minVal, maxVal);
    return dist(rng);
}


double ESSolver::rand01_r()
{
    default_random_engine rng( std::random_device{}() );
    // mt19937 rng(time(0));
	uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}


int * ESSolver::sortArray(double *_val, int size)
{   
    priority_queue < pair<double, int>, vector< pair<double, int> >, greater< pair<double, int> > > pQueue;
  
    for (int i = 0; i < size; i++) 
    {
		pQueue.push(pair<double, int>(_val[i], i));
    }
  
    int *idx = new int[mu];
  
	for (int i = 0; i < mu; i++) 
	{
		idx[i] = pQueue.top().second;
		pQueue.pop();
	}
    
    return idx;
}


void ESSolver::updatePopSigma(int *_idx, int _selSt, int _mutOp)
{
    if (_selSt == 0) // Selection Strategy Type 0 - Comma Selection
    {
		for (int i = 0; i < mu; i++)
		{
			// Updating object variable
			for (int j = 0; j < nVar; j++)
			{
				muPop[i][j] = lambdaPop[_idx[i]][j];
			}
			
			// Updating strategy parameter
			if (_mutOp == 0) 	// Mutation Type 0 - Uncorrelated Mutation with 1 Step Size
			{
			    muSigma[i] = lambdaSigma[_idx[i]];
			}
			else				// Mutation Type 1 - Uncorrelated Mutation with n Step Size
			{
				for (int j = 0; j < nVar; j++)
			    {
				    muSigma2[i][j] = lambdaSigma2[_idx[i]][j];
			    }
			}
		}
	}
	else 		  // Selection Strategy Type 1 - Plus Selection
	{
		for (int i = 0; i < mu; i++)
		{
			// Updating object variable
			for (int j = 0; j < nVar; j++)
			{
				muPop[i][j] = mulambdaPop[_idx[i]][j];
			}
			
			// Updating strategy parameter
			if (_mutOp == 0)	// Mutation Type 0 - Uncorrelated Mutation with 1 Step Size
			{
			    muSigma[i] = mulambdaSigma[_idx[i]];
			}
			else				// Mutation Type 1 - Uncorrelated Mutation with n Step Size
			{
				for (int j = 0; j < nVar; j++)
			    {
				    muSigma2[i][j] = mulambdaSigma2[_idx[i]][j];
			    }
			}
		}
	}
}


double * ESSolver::copyVec(int row, vector<vector<double>> *lPop)
{
    double *trialVec = new double[nVar];

    for (int j = 0; j < nVar; j++)
    {
        trialVec[j] = (*lPop)[row][j];
    }
    
    return trialVec;
}


double * ESSolver::copyDim(int col, vector<vector<double>> *lPop, int size)
{
    double *trialVec = new double[size];

    for (int j = 0; j < size; j++)
    {
        trialVec[j] = (*lPop)[j][col];
    }
    
    return trialVec;
}


void ESSolver::mergeMuLbd(double *_val, int _mutOp)
{
    for (int i = 0; i < mu+lbd; i++)
    {
    	if (i < mu) 	
		{
    	    for (int j = 0; j < nVar; j++) {
    	        mulambdaPop[i][j] = muPop[i][j]; 
    	    }
    
            _val[i] = objFun( copyVec(i, &muPop) );
    	    
            if (_mutOp == 0)	// Mutation Type 0 - Uncorrelated Mutation with 1 Step Size
            {
                mulambdaSigma[i] = muSigma[i];
            }    
            else				// Mutation Type 1 - Uncorrelated Mutation with n Step Size
            {
				//cout << "\n i = " << i;
                for (int j = 0; j < nVar; j++) {
    	            mulambdaSigma2[i][j] = muSigma2[i][j]; 
    	        } 	
            }
        }
        else 			
		{
        	for (int j = 0; j < nVar; j++) {
    	        mulambdaPop[i][j] = lambdaPop[i - mu][j]; 
    	    }
    	        
            _val[i] = objFun( copyVec(i - mu, &lambdaPop) );
				
            if (_mutOp == 0)	// Mutation Type 0 - Uncorrelated Mutation with 1 Step Size
            {
                 mulambdaSigma[i] = lambdaSigma[i - mu];
            } 
            else				// Mutation Type 1 - Uncorrelated Mutation with n Step Size
            {
            	for (int j = 0; j < nVar; j++) {
    	            mulambdaSigma2[i][j] = lambdaSigma2[i - mu][j]; 
    	        }
            }  
        }
    }
}