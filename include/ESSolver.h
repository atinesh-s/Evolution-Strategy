#include <vector>

using namespace std;

extern int maxFE;
extern int minVal;
extern int maxVal;
extern int selStrategy;
extern int recOpO;
extern int recOpS;
extern int mutOp;
extern double epsilon;

class ESSolver
{
public:

	ESSolver(int, int, int, int, int, int);
	~ESSolver(void);
	
    void Init();										// Initialization
    void Optimize();
    void Recombination(int, int, int, int, int);
    void Mutation(int);
    void Selection(int, int);
    
    double *getgBestPos();
    double getgBestVal();
    
    double objFun (const double *);
   
    int rand_r(int, int);								 // Generates random number between min and max
    double rand01_r();									 // Generates random number between 0 and 1
    
    int *sortArray(double *, int);						 // Sorts array
    void updatePopSigma(int *, int, int);
    double *copyVec(int, vector<vector<double>> *);		 // Copies a row
    double *copyDim(int, vector<vector<double>> *, int); // Copies dimension value accross whole population
    void mergeMuLbd(double *, int);						 // Merge Populations
    void sampleRhoPop(int);							 	 // Sample rho individuals out of mu
    void sampleRhoSigma(int, int);
	double avg(double *, int, int);						 // Computes Average

    int function_eval;
    double evaluation(const double *_pos, int function_eval);	//Function Evaluation Selector
    double best_val;

private:
    int nVar;
    int ne;
    int mu;
    int lbd;
    int rhoO;
    int rhoS;
    
	double gBestVal;   			 			// Global Best Value
	
	vector<vector<double>> muPop; 			// Object Variable
	vector<double> muSigma;					// Strategy Parameter - Mutation Type 0
	vector<vector<double>> muSigma2;		// Strategy Parameter - Mutation Type 1

	vector<vector<double>> lambdaPop; 		// Object Variable
	vector<double> lambdaSigma; 			// Strategy Parameter - Mutation Type 0
	vector<vector<double>> lambdaSigma2; 	// Strategy Parameter - Mutation Type 1

	vector<vector<double>> mulambdaPop; 	// Object Variable
	vector<double> mulambdaSigma;			// Strategy Parameter - Mutation Type 0
	vector<vector<double>> mulambdaSigma2;	// Strategy Parameter - Mutation Type 1

	vector<vector<double>> rhoPop; 			// Object Variable
	vector<double> rhoSigma;				// Strategy Parameter - Mutation Type 0
	vector<vector<double>> rhoSigma2;		// Strategy Parameter - Mutation Type 1
};
