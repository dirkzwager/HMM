//init
#define INITW 	.30
#define INITR	.40
#define INITC 	.30
//states
#define WARM 	0
#define RAIN	1
#define COLD 	2
//observations
#define ONE 	0
#define TWO 	1
#define THREE 	2
//state transitions
#define WtoW	.60
#define WtoR	.40
#define WtoC	.10
#define RtoW	.05
#define RtoR	.40	
#define RtoC	.45
#define CtoW	.30
#define CtoR	.20
#define CtoC	.40
//emmisions
#define WONE	.05
#define WTWO	.15
#define WTHREE	.80
#define RONE	.40
#define RTWO	.50
#define RTHREE	.10
#define CONE	.60
#define CTWO	.30
#define CTHREE	.10
	
//enum States { Warm, Rain,  Cold};
//enum Observations { One, Two, Three };

void printStates(int state) {
	switch (state) {
		case 0: printf(" Warm ");
	 	break;
		case 1: printf(" Rain ");
	 	break;
		case 2: printf(" Cold ");
		break;
	} 
}

void printObservations(int state) {
	switch (state) {
		case 0: printf(" 1 Ice ");
		break;
		case 1: printf(" 2 Ice ");
		break;
		case 2: printf(" 3 Ice ");
		break;
	} 
}

double randfrac() {
	return (double)(rand() % 100 + 1)/100.0;
}

//array containing probabilities for the different states, array length (n different states) and sampeling period.
int *hiddenStateGenerator(int npages) {	
	srand(time(NULL));
	int *states = new int[npages];
	double rand = randfrac();
	int i = 0;
	
	if (rand < INITW)
		states[i] = WARM;
	else
		states[i] = COLD;
	
	for (i=1; i<npages; i++) {
		rand = randfrac();
		switch (states[i-1]) {
			case WARM:
				if (rand < WtoW)
					states[i]  = WARM;
				else states[i] = COLD;
			break;
			case COLD:
				if (rand <= CtoC)
					states[i]  = COLD;
				else states[i] = WARM;
			break;		
		}
	}	
	return states;		
}

int *jacobsDiary(int *weather, int npages) {
	int *icecreams = new int[npages];
	double rand;
	
	for (int i=0; i<npages; i++) {
		rand = randfrac();
		switch(weather[i]) {
			case WARM:
				if (rand < WONE)
					icecreams[i] = ONE;
				else if (rand < WONE+WTWO && rand >= WONE)
					icecreams[i] = TWO;
				else 
					icecreams[i] = THREE;
			break;
			case COLD:
				if (rand < CONE)
					icecreams[i] = ONE;
				else if (rand < CONE+CTWO && rand >= CONE)
					icecreams[i] = TWO;
				else 
					icecreams[i] = THREE;
			break;
		} 
	}
	return icecreams;
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
