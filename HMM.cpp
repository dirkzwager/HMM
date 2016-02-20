#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Eisner.h"
#include "Matrix.h"

using namespace std;

//no separation of defintion & implementation because class is templated..
template <class T>
class HMM {
public:
	//HMM() { S=NULL; A=NULL; B=NULL, M=NULL; p==NULL;	 }
	HMM() { S=NULL; O=NULL; p=NULL; Ns=0; No=0; }
        ~HMM() {;}
	HMM(int *_S, int _Ns, int *_O, int _No, Matrix<T> _A, Matrix<T> _B, T *_p);
	HMM(int *_S, int _Ns, int *_O, int _No, Matrix<T> _A, Matrix<T> _B);//uniform initial distribution
	double observationProbabilitiy   (int *obs, int t);
	void   printOptimalStateSequence (int *obs, int t);	
	double sequenceProbability       (int *observations);
	double	  wrForward  (int *obs, int t);
	double	  wrBackward (int *obs, int t);
	void print();

	double	  Forward    (int *obs, int t);
	double	  rForward   (int *obs, int n, int t);
	double	  rBackward  (int *obs, int n, int t, int tO);
	double	  Backward   (int *obs, int t);
	int*	  Viterbi    (int *obs, int t);
	void 	  BaumWelch  (int *obs, int t);
private:
	T		  emmision   (int state, int obs) 		 { return B.cellVal(obs,state); }//y,x
	T		  transition (int from, int to) 			 { return A.cellVal(to, from); }
	T		  model      (int i, int j) 					 { return M.cellVal(i, j); }
	int		  argmax     (T *arg, int n);
	
	int       *S;		 //array with labels of possible states.
	int       Ns;		 //Number of states
	int		  *O;		 //array with labels of possible Observations
	int 	  No;		 //number of possible observations 
	Matrix<T>  A;		 //Matrix containing state transition probabilities
	Matrix<T>  B;		 //Matrix containing emission probabilities
	Matrix<T>  M;		 //Matrix containing the model with states St
	T		  *p;		 //initial distribution of states.	
	
};

template <class T>
HMM<T>::HMM(int *_S, int _Ns, int *_O, int _No, Matrix<T> _A, Matrix<T> _B, T *_p) {
	S  = _S; 
	Ns = _Ns;
	O  = _O; 
	No = _No; 
	A  = _A;
	B  	= _B;
	p = _p;
}

template <class T>
HMM<T>::HMM(int *_S, int _Ns, int *_O, int _No, Matrix<T> _A, Matrix<T> _B) {
	S  = _S;
	Ns = _Ns;
	A  = _A;
	O  = _O;
	No = _No;
	B  = _B;
	p  = new T[Ns];
	
	for (int i=0; i<Ns; i++)
		p[i] = (T)(1.0/(double)Ns);
}

template <class T> 
double HMM<T>::Forward(int *obs, int t) {
	double *fp = new double[Ns];//array containing forward probabilities for each time step
	for (int i=0; i<Ns; i++) fp[i] = p[i] * emmision(i, obs[0]);//init p*emmision
	
	//for each observation occuring at time tn
	for (int n=1; n<t; n++) {
		double f[Ns]; 
		for (int i=0; i<Ns; i++) f[i] = 0.0;
		for (int s=0; s<Ns; s++) {
			for (int k=0; k<Ns; k++) { 
				double tmp = (fp[k] * transition(k,s)) * emmision(s,obs[n]);
				f[s] += tmp;			 	
			 }
		}
		for (int i=0; i<Ns; i++)
			fp[i] = f[i];
	}
	double sum = 0.0;
	for (int i=0; i<Ns; i++) sum += fp[i];
    return sum;	
}

template <class T> 
double HMM<T>::wrForward(int *obs, int t) {
	double sum = 0.0;
	
	for (int n=0; n<Ns; n++) {
		double f = rForward(obs, n, t);
		sum += f;
		cout << "f:" << f << " sum " << sum << endl;
	}
	return sum;
}

template <class T> 
double HMM<T>::rForward(int *obs, int n, int t) {
	double sum = 0.0;
	
	if (t == 1) 
		return p[n] * emmision(n, obs[0]);		
	for (int k=0; k<Ns; k++) {
		double f = rForward(obs,k,t-1) * transition(k,n);
		sum += f;
	//	cout << "f:" << f << " sum " << sum << endl;
	}
	return sum * emmision(n,obs[t-1]);
} 

template <class T> 
double HMM<T>::wrBackward(int *obs, int t) {
	double sum = 0.0;
	
	for (int n=0; n<Ns; n++) 
		sum += p[n] * rBackward(obs, n, t-1, 0); 
	return sum;
}

template <class T> 
double HMM<T>::rBackward(int *obs, int n, int t, int tO) {
	cout << tO << endl;
	double sum = 0.0;
	
	if (tO == t)
		return 1.0;
	for (int k=0; k<Ns; k++) 
		sum += rBackward(obs,k,t,tO+1) * transition(n,k) * emmision(k,obs[tO+1]);
	return sum;
} 

template <class T> 
double HMM<T>::Backward(int *obs, int t) {
	double *fp = new double[Ns]; 
	for (int i=0; i<Ns; i++) fp[i] = 1.0;
	
	for (int n=t-1; n>0; n--) {
		double f[Ns]; for (int i=0; i<Ns; i++) f[i] = 0.0;
		for (int s=0; s<Ns; s++) {
			for (int k=0; k<Ns; k++) {
				double tmp = (fp[k] * transition(s,k)) * emmision(k,obs[n]);
				f[s] += tmp;
 			 	//cout << endl << "fp:" << tmp << " = " << "fp-1: " << fp[k] << "\t  *\t" << " A("; 
			 	//printStates(s); cout << ","; printStates(k); cout << "): " <<  transition(s,k) << "\t  *\tB(";
			 	//printStates(k); cout << ","; printObservations(obs[n]); cout << "): " << emmision(k, obs[n]);
			}
		}
		for (int i=0; i<Ns; i++) 
			fp[i] = f[i];
	}
	double sum = 0.0;
	for (int i=0; i<Ns; i++) sum += fp[i]*p[i]*emmision(i, obs[0]);
    return sum;	
}

template <class T>
int* HMM<T>::Viterbi(int *obs, int t) {
	int *path = new int[t]; 
	T *fp = new T[Ns]; for (int i=0; i<Ns; i++) fp[i] = p[i] * emmision(i, obs[0]);
	
	path[0] = argmax(fp,Ns);//get state with largest begin value
	
	for (int n=1; n<t; n++) {	
		double pt[Ns]; for (int i=0; i<Ns; i++) pt[i] = 0.0;
		double largest;
		int idx = -1;
		for (int s=0; s<Ns; s++) {
			for (int k=0; k<Ns; k++) {
				double f = (fp[k] * transition(k,s)) * emmision(s, obs[n]);
				if (idx == -1) { 
					largest = f;
					idx = 0;
				}		
				else if (f>largest) {
					largest = f;
					idx = s;
				}
				pt[s] += f;
			}
		}
		for (int i=0; i<Ns; i++) fp[i] = pt[i];			
		path[n] = idx;
	}

	return path;
}

template <class T>
void HMM<T>::BaumWelch(int *obs, int t) {
//	Matrix M(col, row);
	
	for (int n=0; n<t; n++) {
	
	}
}

//returns array index pointing to largest value in array.
template <class T> 
int HMM<T>::argmax(T *arg, int n) {
	int idx = 0;
	for (int i=1; i<n; i++) 
		if (arg[i] > arg[idx])
			idx = i;
	return idx;		
}

template <class T> 
double HMM<T>::observationProbabilitiy(int *obs, int t) {
	cout << endl << "=================================Forward:====================================" << endl;
	double f = Forward(obs,t);
	cout << endl << endl << "================================Backward:====================================" << endl;
	double b = Backward(obs,t);
	cout << endl << "==============================================================================" << endl;
	return b;
}

template <class T> 
void HMM<T>::printOptimalStateSequence(int *obs, int t) {
	int * seq = Viterbi(obs, t);
	double fp = Forward(obs,t);
	cout << endl << "Voor de observatie: " << endl;
	for (int i=0; i<t; i++)
		printObservations(obs[i]);
	cout << endl << endl << "is de optimale rij van states (p=" << fp << "): " << endl;
	for (int i=0; i<t; i++) 
		printStates(seq[i]);
	delete[] seq;
}

template <class T>
void HMM<T>::print() {
	cout << "+==============================+" << endl;
	printStates(S[2]);
	
	printf("Statess:\n");
	for (int i=0; i<Ns; i++)
		printStates(S[i]);
	printf("\nObservations:\n");
	for (int i=0; i<No; i++)
		printObservations(O[i]);
	printf("\n\nA:\n");
	A.printMatrix();
	printf("\n\nB:\n");
	B.printMatrix();
	printf("\n\nM:\n");
	M.printMatrix();
	printf("\n\nInit:\n");
	for (int i=0; i<Ns; i++)
		printf("%f ", p[i]);
	cout << endl << "+==============================+" << endl;
}

//main to test the model
int main() {
'
.,	const int NPAGES = 3;

	//description of states
	const int NSTATES = 3;
	int *States = new int[NSTATES];
	States[0] = 0; States[1] = 1; States[2] = 2;

	//description of observations
	const int NOBS = 3;
	int *Obs = new int[NOBS];
	Obs[0]=ONE;	Obs[1]=TWO;	Obs[2]=THREE;

	//init Transition & Emmision Matrix
	Matrix<double> Transition(3,3);			//2 cols, 2 rows
	Matrix<double> Emmision(3,3);			//3 cols, 2 rows

	//init Transition, Emmision & Init (A,B,p) data arrays
	double *init 		= new double[NSTATES];
	double *transition 	= new double[NSTATES*NSTATES];	//fill data for Transition
	double *emmision 	= new double[NSTATES*NOBS];		//fill data for Emmision	

	//set init values 
	init[0] = INITW; init[1] = INITR; init[2] = INITC;
	
	//set transition data probabilities
	transition[0] = WtoW; transition[3] = WtoR; transition[6] = WtoC;
	transition[1] = RtoW; transition[4] = RtoR; transition[7] = RtoC;
	transition[2] = CtoW; transition[5] = CtoR; transition[8] = CtoC;

	//set emmision data probabilities
	emmision[0] = WONE; emmision[3] = WTWO; emmision[6] = WTHREE;
	emmision[1] = RONE; emmision[4] = RTWO; emmision[7] = RTHREE;
	emmision[2] = CONE; emmision[5] = CTWO; emmision[8] = CTHREE;

	//fill Transition & Emmision matrices with appropriate probabilities.
	Transition.fillMatrix(transition, NSTATES*NSTATES);
	Emmision.fillMatrix(emmision, NSTATES*NOBS);
	
	//define HMM Model.
	HMM<double> hmm(States,NSTATES,Obs,NOBS,Transition,Emmision,init);	
	hmm.print();
	
	int *observations = new int[3];
	observations[0] = (int)ONE; observations[1] = (int)TWO; observations[2] = (int)THREE; //observations[3] = (int)THREE;
	hmm.printOptimalStateSequence(observations, 3);
	
	return 0;
}

/*

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

JacobEisner.h
#define WARM 	0
#define COLD 	1
//observations
#define ONE 	0
#define TWO 	1
#define THREE 	2
//state transitions
#define WtoW	.70
#define WtoC	.30
#define CtoC	.60
#define CtoW	.40
//emmisions
#define WONE	.2
#define WTWO	.4
#define WTHREE	.4
#define CONE	.5
#define CTWO	.4
#define CTHREE	.1

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

main.cpp	

int main() {
	//number of pages in Jacob Eisners Dairy.
	const int NPAGES = 3;

	//description of states
	const int NSTATES = 2;
	int *States = new int[NSTATES];
	States[0] = WARM; States[1] = COLD;

	//description of observations
	const int NOBS = 3;
	int *Obs = new int[NOBS];
	Obs[0]=ONE;	Obs[1]=TWO;	Obs[2]=THREE;

	//init Transition & Emmision Matrix
	Matrix<double> Transition(2,2);			//2 cols, 2 rows
	Matrix<double> Emmision(3,2);			//3 cols, 2 rows

	//init Transition, Emmision & Init (A,B,p) data arrays
	double *transition 	= new double[NSTATES*NSTATES];	//fill data for Transition
	double *emmision 	= new double[NSTATES*NOBS];		//fill data for Emmision
	double *init 		= new double[NSTATES];

	//set transition data probabilities
	transition[0] = WtoW; transition[2] = WtoC;
	transition[1] = CtoW; transition[3] = CtoC;

	//set emmision data probabilities
	emmision[0] = WONE; emmision[2] = WTWO; emmision[4] = WTHREE;
	emmision[1] = CONE; emmision[3] = CTWO; emmision[5] = CTHREE;

	//set init values 
	init[0] = INITW; init[1] = INITC;

	//fill Transition & Emmision matrices with appropriate probabilities.
	Transition.fillMatrix(transition, 4);
	Emmision.fillMatrix(emmision, 6);
	
	//define HMM Model.
	HMM<double> hmm(States,NSTATES,Obs,NOBS,Transition,Emmision,init);	
	hmm.print();
	
	int *page = new int[NPAGES];
	page[0]=THREE; page[1]=ONE; page[2]=THREE;
	
	double obsprob = hmm.observationProbabilitiy(page,3);
	
	cout << endl << "Probability for Observation Sequence: 3,1,3: " << obsprob << endl;
	
	return 0;
}

/*

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Output from main.cpp.


+==============================+
States:
 Warm  Cold
Observations:
 1 Ice  2 Ice  3 Ice

A:
| 0.7   0.3 |
| 0.4   0.6 |


B:
| 0.2   0.4     0.4 |
| 0.5   0.4     0.1 |

Init:
0.800000 0.200000
+==============================+

=================================Forward:====================================

fp:0.0448=fp-1: 0.32      *      A( Warm , Warm ): 0.7    *     B( Warm , 1 Ice ): 0.2
fp:0.0016=fp-1: 0.02      *      A( Cold , Warm ): 0.4    *     B( Warm , 1 Ice ): 0.2
F[s]: 0.0464

fp:0.048=fp-1: 0.32       *      A( Warm , Cold ): 0.3    *     B( Cold , 1 Ice ): 0.5
fp:0.006=fp-1: 0.02       *      A( Cold , Cold ): 0.6    *     B( Cold , 1 Ice ): 0.5
F[s]: 0.054

fp:0.012992=fp-1: 0.0464  *      A( Warm , Warm ): 0.7    *     B( Warm , 3 Ice ): 0.4
fp:0.00864=fp-1: 0.054    *      A( Cold , Warm ): 0.4    *     B( Warm , 3 Ice ): 0.4
F[s]: 0.021632

fp:0.001392=fp-1: 0.0464  *      A( Warm , Cold ): 0.3    *     B( Cold , 3 Ice ): 0.1
fp:0.00324=fp-1: 0.054    *      A( Cold , Cold ): 0.6    *     B( Cold , 3 Ice ): 0.1
F[s]: 0.004632

================================Backward:====================================

fp:0.28 = fp-1: 1         *      A( Warm , Warm ): 0.7    *     B( Warm , 3 Ice ): 0.4
fp:0.03 = fp-1: 1         *      A( Warm , Cold ): 0.3    *     B( Cold , 3 Ice ): 0.1
fp:0.16 = fp-1: 1         *      A( Cold , Warm ): 0.4    *     B( Warm , 3 Ice ): 0.4
fp:0.06 = fp-1: 1         *      A( Cold , Cold ): 0.6    *     B( Cold , 3 Ice ): 0.1
fp:0.0434 = fp-1: 0.31    *      A( Warm , Warm ): 0.7    *     B( Warm , 1 Ice ): 0.2
fp:0.033 = fp-1: 0.22     *      A( Warm , Cold ): 0.3    *     B( Cold , 1 Ice ): 0.5
fp:0.0248 = fp-1: 0.31    *      A( Cold , Warm ): 0.4    *     B( Warm , 1 Ice ): 0.2
fp:0.066 = fp-1: 0.22     *      A( Cold , Cold ): 0.6    *     B( Cold , 1 Ice ): 0.5
==============================================================================

Probability for Observation Sequence: 3,1,3: 0.026264

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*/
