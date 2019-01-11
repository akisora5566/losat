////////////////////////////////////////////////////////////////////////////////
// 3-value Logic Simulator, written by Michael Hsiao
//   Began in 1992, revisions and additions of functions till 2013
/////////////////////////////////////////////////////////////////////////////////
//  Functions added to implement Genetic Algorithm and Controllability-based Partitioning
/////////////////////////////////////////////////////////////////////////////////
#include <cstdlib>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <ctime>
#include <vector>
#include <map>
#include <math.h>
using namespace std;

#define HZ 100
#define RETURN '\n'
#define EOS '\0'
#define COMMA ','
#define SPACE ' '
#define TAB '\t'
#define COLON ':'
#define SEMICOLON ';'

#define SPLITMERGE 'M'

#define T_INPUT 1
#define T_OUTPUT 2
#define T_SIGNAL 3
#define T_MODULE 4
#define T_COMPONENT 5
#define T_EXIST 9
#define T_COMMENT 10
#define T_END 11

#define TABLESIZE 5000
#define MAXIO 5000
#define MAXMODULES 5000
#define MAXDFF 10560

#define GOOD 1
#define FAULTY 2
#define DONTCARE -1
#define ALLONES 0xffffffff

#define MAXlevels 10000
#define MAXIOS 5120
#define MAXFanout 10192
#define MAXFFS 40048
#define MAXGATES 100000
#define MAXevents 100000

#define TRUE 1
#define FALSE 0

#define EXCITED_1_LEVEL 1
#define POTENTIAL 2
#define LOW_DETECT 3
#define HIGH_DETECT 4
#define REDUNDANT 5

enum
{
   JUNK,           /* 0 */
   T_input,        /* 1 */
   T_output,       /* 2 */
   T_xor,          /* 3 */
   T_xnor,         /* 4 */
   T_dff,          /* 5 */
   T_and,          /* 6 */
   T_nand,         /* 7 */
   T_or,           /* 8 */
   T_nor,          /* 9 */
   T_not,          /* 10 */
   T_buf,          /* 11 */
   T_tie1,         /* 12 */
   T_tie0,         /* 13 */
   T_tieX,         /* 14 */
   T_tieZ,         /* 15 */
   T_mux_2,        /* 16 */
   T_bus,          /* 17 */
   T_bus_gohigh,   /* 18 */
   T_bus_golow,    /* 19 */
   T_tristate,     /* 20 */
   T_tristateinv,  /* 21 */
   T_tristate1     /* 22 */
};


////////////////////////////////////////////////////////////////////////
// gateLevelCkt class
////////////////////////////////////////////////////////////////////////

class gateLevelCkt
{
    // circuit information
    int numgates;	// total number of gates (faulty included)
    int numFaultFreeGates;	// number of fault free gates
    int numout;		// number of POs
    int maxlevels;	// number of levels in gate level ckt
    int maxLevelSize;	// maximum number of gates in one given level
    int levelSize[MAXlevels];	// levelSize for each level
    int inputs[MAXIOS];
    int outputs[MAXIOS];
    int ff_list[MAXFFS];
    int *ffMap;
    unsigned char *gtype;	// gate type
    short *fanin;		// number of fanin, fanouts
    short *fanout;
    int *levelNum;		// level number of gate
    unsigned *po;
    int **inlist;		// fanin list
    int **fnlist;		// fanout list
    char *sched;		// scheduled on the wheel yet?
    unsigned int *value1;
    unsigned int *value2;	// value of gate
    int **predOfSuccInput;      // predecessor of successor input-pin list
    int **succOfPredOutput;     // successor of predecessor output-pin list

    // for simulator
    int **levelEvents;	// event list for each level in the circuit
    int *levelLen;	// evenlist length
    int numlevels;	// total number of levels in wheel
    int currLevel;	// current level
    int *activation;	// activation list for the current level in circuit
    int actLen;		// length of the activation list
    int *actFFList;	// activation list for the FF's
    int actFFLen;	// length of the actFFList

	// reset signal masking
	int* reset_power0;
	int* reset_power1;
	double* prob_maskto0;
	double* prob_maskto1;
	void FindResetInputs(); // calculate reset signal masking related values
	// called by constructor only

public:
    int numff;		// number of FF's
    int numpri;
    unsigned int *RESET_FF1;	// value of reset ffs read from *.initState
    unsigned int *RESET_FF2;	// value of reset ffs read from *.initState

    gateLevelCkt(string);	// constructor
    void setFaninoutMatrix();	// builds the fanin-out map matrix

    void applyVector(char *);	// apply input vector
    unsigned int* getValue1FFs();
	unsigned int* getValue2FFs();
    // simulator information
    void setupWheel(int, int);
    void insertEvent(int, int);
    int retrieveEvent();
    void goodsim();		// logic sim (no faults inserted)
	void reset();
    void setTieEvents();	// inject events from tied nodes

    void observeOutputs();	// print the fault-free outputs
    void observeFFs();
    void printGoodSig(ofstream, int);	// print the fault-free outputs to *.sig
    double* getProbMaskto0();
    double* getProbMaskto1();


    char *goodState;		// good state (without scan)
};

//////////////////////////////////////////////////////////////////////////
// Global Variables

char vectr[5120];
gateLevelCkt *circuit;
int OBSERVE, INIT0;
int vecNum=0;
int numTieNodes;
int TIES[512];
std::vector<float> count_zero_ffs, count_one_ffs, ff_bias;
std::vector<int> state_partition;
std::map< std::vector<char>, int > state_table0, state_table1, state_table2, state_table3, state_table4;

//////////////////////////////////////////////////////////////////////////
// Functions start here

std::vector<std::vector<int> > generate_random_test_sequence(int num_inputs, int num_vec){
    std::vector<std::vector<int> > test_sequence;
    std::vector<int> temp;
    int a;
    double probTo0, probTo1, val;
    double* probMaskTo0 = circuit->getProbMaskto0();
    double* probMaskTo1 = circuit->getProbMaskto1();
	for (int i =0; i < num_vec; i++){
        for (int j=0; j < num_inputs; j++){
        	//  UNCOMMENT THIS CODE FOR RESET SIGNAL MASKING////////////////////////////////
        	//  COMMENT THIS CODE TO REMOVE RESET SIGNAL MASKING //////////////////////////////
        	/*
			probTo0 = *(probMaskTo0 + j);
            probTo1 = *(probMaskTo1 + j);
            if(probTo0 > 0.5 || probTo1 > 0.5){
            	if(probTo0 > probTo1){
            		val = (double)rand()/ RAND_MAX;//probTo0;
            	if(val < probTo0){
            		a = 0;		
				}else{
					a = 1;
				}
			}else{
				val = rand()/ RAND_MAX;//probTo0;
            	if(val < probTo1){
            		a = 1;		
				}else{
					a = 0;
				}
			}	
			}else{
				a = rand()%2;
			}
			*/
			/////////////////////////////////////
			 a = rand()%2;     ////COMMENT THIS IF YOU ARE IMPLEMENRTING RESET SIGNAL MASKING
            						//uncomment this if you are running without reset signal masking
			 //////////////////////////////////////////////////////////
            temp.push_back(a);
        }
        test_sequence.push_back(temp);
        temp.clear();
    }
    return test_sequence;
}


//////////////////////////////////////////////////////////////////////////
// returns 1 if a regular vector, returns 2 if a scan
//////////////////////////////////////////////////////////////////////////
int getVector(ifstream &inFile, int vecSize)
{
    int i;
    char thisChar;

    inFile >> thisChar;
    while ((thisChar == SPACE) || (thisChar == RETURN))
	inFile >> thisChar;

    vectr[0] = thisChar;
    if (vectr[0] == 'E')
	return (0);

    for (i=1; i<vecSize; i++)
	inFile >> vectr[i];
    vectr[i] = EOS;

    return(1);
}

// logicSimFromFile() - logic simulates all vectors in vecFile on the circuit
int logicSimFromFile(ifstream &vecFile, int vecWidth)
{
    int moreVec;

    moreVec = 1;
    while (moreVec)
    {
        moreVec = getVector(vecFile, vecWidth);
        if (moreVec == 1)
        {
	        
        	circuit->applyVector(vectr);
            circuit->goodsim();      // simulate the vector
	        
            vecNum++;
        }  // if (moreVec == 1)
    }   // while (getVector...)

    return (vecNum);
}

// main()
void simulateRandomLogic(string vecName, int num_vectors){
	ifstream vecFile;
	ofstream vecFile_w;
	int totalNumVec;
	int vecWidth;
	int a;
	double probTo0, probTo1, val;
    double* probMaskTo0 = circuit->getProbMaskto0();
    double* probMaskTo1 = circuit->getProbMaskto1();
	vecFile_w.open(vecName.c_str(), ios::out);
    if (!vecFile_w)
    {
	cerr << "Can't open " << vecName << "\n";
	exit(-1);
    }
    vecFile_w << circuit->numpri;
    vecFile_w << "\n";
    for(int i =0; i < num_vectors; i++){
    	for(int j = 0; j < circuit->numpri; j++){
    		/*
    		probTo0 = *(probMaskTo0 + j);
            probTo1 = *(probMaskTo1 + j);
            if(probTo0 > 0.5 || probTo1 > 0.5){
            	if(probTo0 > probTo1){
            		val = (double)rand()/ RAND_MAX;//probTo0;
            	if(val < probTo0){
            		a = 0;		
				}else{
					a = 1;
				}
			}else{
				val = rand()/ RAND_MAX;//probTo0;
            	if(val < probTo1){
            		a = 1;		
				}else{
					a = 0;
				}
			}	
			}else{
				a = rand()%2;
			}
			vecFile_w << a;*/
							//COMMENT TO REMOVE RSM
    		vecFile_w << rand()%2;	//UNCOMMENT TO REMOVE RSM
		}
		vecFile_w << "\n";
	}
	vecFile_w << "END\n";
    vecFile_w.close();
    
	//CHANGE HERE-----------------------------------
    vecFile.open(vecName.c_str(), ios::in);
    if (!vecFile)
    {
	cerr << "Can't open " << vecName << "\n";
	exit(-1);
    }

    vecFile >> vecWidth;
	circuit->reset();
    circuit->setTieEvents();
    totalNumVec = logicSimFromFile(vecFile, vecWidth);
    vecFile.close();
	
}


void simulateLogic(string vecName, int num_vectors, std::vector< std::vector< int> > test_vec){
	ifstream vecFile;
	ofstream vecFile_w;
	int totalNumVec;
	int vecWidth;
	vecFile_w.open(vecName.c_str(), ios::out);
    if (!vecFile_w)
    {
	cerr << "Can't open " << vecName << "\n";
	exit(-1);
    }
    vecFile_w << circuit->numpri;
    vecFile_w << "\n";
    for(int i =0; i < num_vectors; i++){
    	for(int j = 0; j < circuit->numpri; j++){
    		vecFile_w << test_vec[i][j];
		}
		vecFile_w << "\n";
	}
	vecFile_w << "END\n";
    vecFile_w.close();
    
	//CHANGE HERE-----------------------------------
    vecFile.open(vecName.c_str(), ios::in);
    if (!vecFile)
    {
	cerr << "Can't open " << vecName << "\n";
	exit(-1);
    }

    vecFile >> vecWidth;
	circuit->reset();
    circuit->setTieEvents();
    totalNumVec = logicSimFromFile(vecFile, vecWidth);
    vecFile.close();
	
}


void partition(int num_vectors, string vecName){
	unsigned int* value1ffs;
    unsigned int* value2ffs;
    
    
	for(int count =0; count < 2000; count++){
    	
		simulateRandomLogic(vecName, num_vectors);
        value1ffs = circuit->getValue1FFs();
        value2ffs = circuit->getValue2FFs();
				 
        for(int i=0; i < circuit->numff; i++){
        	if(*value1ffs == *value2ffs){
        		if(*value1ffs == 0){
        			count_zero_ffs[i]++;
				}else{
        			count_one_ffs[i]++;
				} 
			}
			value1ffs++;	value2ffs++;
        }
    }	
	
	for(int i = 0; i < circuit->numff; i++){
		if(count_one_ffs[i] > count_zero_ffs[i]){
			ff_bias.push_back((count_one_ffs[i] - count_zero_ffs[i])/2000);
		}else{
			ff_bias.push_back((count_zero_ffs[i] - count_one_ffs[i])/2000);
		}
        //cout << ff_bias[i] <<"\n";
        if(ff_bias[i] <= 0.2){
            state_partition.push_back(0);
        }else if(ff_bias[i] <= 0.4){
            state_partition.push_back(1);
        }else if(ff_bias[i] <=0.6  ){
            state_partition.push_back(2);
        }else if(ff_bias[i] <= 0.8){
            state_partition.push_back(3);
        }else if(ff_bias[i] <= 1){
            state_partition.push_back(4);
        }
    }
	
	
	cout << "\n\n";
	for(int i =0; i < state_partition.size(); i++){
		cout << state_partition[i];
	}	
}

void increment_state_table(string vecName, std::vector < std::vector <int> > test_vec,int num_vectors){
	simulateLogic(vecName, num_vectors, test_vec);
	int num_inputs = circuit -> numpri;
	int num_ff = circuit -> numff;
	std::vector< char> curr_state_0, curr_state_1, curr_state_2, curr_state_3, curr_state_4;
	unsigned int* value1ffs = circuit->getValue1FFs();
	unsigned int* value2ffs = circuit->getValue2FFs();
	std::vector< char> states;
	
	for(int ff=0; ff < num_ff; ff++){
			if(*value1ffs == *value2ffs){
        		if(*value1ffs == 0){
        			states.push_back('0');cout <<"0";
				}else{
        		states.push_back('1');cout <<"1";
				} 
			}else{
				states.push_back('X');cout << "X";
				
			}
			value1ffs++; value2ffs++;
        }
    
	float f = 0;  
    for(int j =0; j < states.size(); j++){
        		//cout << j <<"\n";
            	if(state_partition[j] == 0){
                	curr_state_0.push_back(states[j]);
            	}else if(state_partition[j] == 1){
                	curr_state_1.push_back(states[j]);
            	}else if(state_partition[j] == 2){
                	curr_state_2.push_back(states[j]);
            	}else if(state_partition[j] == 3){
                	curr_state_3.push_back(states[j]);
            	}else if(state_partition[j] == 4){
                	curr_state_4.push_back(states[j]);
            	}
        	}
      
	  std::map< std::vector<char>, int >::iterator it;

        	if(state_table0.count(curr_state_0)>0){
            	it = state_table0.find(curr_state_0);
            	f += 2*(1/(float)it->second);
            	it->second++;
        	}else{
        		state_table0[curr_state_0] = 1;
        		f += 20;		//MAYBE THIS HELPS
        	}

        	if(state_table1.count(curr_state_1)>0){
        		it = state_table1.find(curr_state_1);
            	f += 4*(1/(float)it->second);
            	it->second++;
        	}else{
            	state_table1[curr_state_1] = 1;
            	f += 4*10;
        	}

        	if(state_table2.count(curr_state_2)>0){
        		it = state_table2.find(curr_state_2);
            	f += 6*(1/(float)it->second);
            	it->second++;
        	}else{
            	state_table2[curr_state_2] = 1;
            	f += 6*10;
        	}

        	if(state_table3.count(curr_state_3)>0){
        		it = state_table3.find(curr_state_3);
            	f += 8*(1/(float)it->second);
            	it->second++;
        	}else{
            	state_table3[curr_state_3] = 1;
            	f += 80;
        	}

        	if(state_table4.count(curr_state_4)>0){
        		it = state_table4.find(curr_state_4);
            	f += 15*(1/(float)it->second);
            	it->second++;
        	}else{
            	state_table4[curr_state_4] = 1;
            	f += 150;
        	}  	
        	cout <<"\n Increment done for this state\n";
        	cout << "Fitness of this vector is :" << f;
        //	int cont; cin>>cont;
    		
}

std::vector< float> compute_fitness(std::vector< std::vector< std::vector <int> > > curr_population, int num_vectors, string vecName){
	std::vector< float> fitness;
	int num_inputs = circuit -> numpri;
	int num_ff = circuit -> numff;
	int skip_flag;
	float f =0;
    std::vector< char> states;
    unsigned int* value1ffs; unsigned int* value2ffs;
	std::vector< char> curr_state_0, curr_state_1, curr_state_2, curr_state_3, curr_state_4;
    char** test_vec1;
    test_vec1 = new char* [num_vectors];
    for( int i =0; i < curr_population.size() ; i++){
    	skip_flag =0;
        
        simulateLogic(vecName, num_vectors, curr_population[i]);
        
		value1ffs = circuit->getValue1FFs();
        value2ffs = circuit->getValue2FFs();
        
		for(int ff=0; ff < num_ff; ff++){
			if(*value1ffs == *value2ffs){
        		if(*value1ffs == 0){
        			states.push_back('0');//cout <<"0";
				}else{
        		states.push_back('1');//cout <<"1";
				} 
			}else{
				states.push_back('X');
				skip_flag = 1;
			}
			value1ffs++; value2ffs++;
        }
        //cout << "\n";
        f = 0; 
        if(!skip_flag){
        	for(int j =0; j < states.size(); j++){
        		//cout << j <<"\n";
            	if(state_partition[j] == 0){
                	curr_state_0.push_back(states[j]);
            	}else if(state_partition[j] == 1){
                	curr_state_1.push_back(states[j]);
            	}else if(state_partition[j] == 2){
                	curr_state_2.push_back(states[j]);
            	}else if(state_partition[j] == 3){
                	curr_state_3.push_back(states[j]);
            	}else if(state_partition[j] == 4){
                	curr_state_4.push_back(states[j]);
            	}
        	}
        
        	//cout << i << "after divide states\t";

        	std::map< std::vector<char>, int >::iterator it;

        	if(state_table0.count(curr_state_0)>0){
            	it = state_table0.find(curr_state_0);
            	f += 2*(1/(float)it->second);
            	//it->second++;
        	}else{
        		//state_table0[curr_state_0] = 1;
        		f += 2*10;		
        	}

        	if(state_table1.count(curr_state_1)>0){
        		it = state_table1.find(curr_state_1);
            	f += 4*(1/(float)it->second);
            	//it->second++;
        	}else{
            	//state_table1[curr_state_1] = 1;
            	f += 4*10;
        	}

        	if(state_table2.count(curr_state_2)>0){
        		it = state_table2.find(curr_state_2);
            	f += 8*(1/(float)it->second);
            	//it->second++;
        	}else{
            	//state_table2[curr_state_2] = 1;
            	f += 8*10;
        	}

        	if(state_table3.count(curr_state_3)>0){
        		it = state_table3.find(curr_state_3);
            	f += 16*(1/(float)it->second);
            	//it->second++;
        	}else{
            	//state_table3[curr_state_3] = 1;
            	f += 16*10;
        	}

        	if(state_table4.count(curr_state_4)>0){
        		it = state_table4.find(curr_state_4);
            	f += 40*(1/(float)it->second);
            	//it->second++;
        	}else{
            	//state_table4[curr_state_4] = 1;
            	f += 40*10;
        	}
        
        	//cout << i << "after compute fitness\t";
        	//cout << "Fitness of input vector " << i << " is " << f <<"\n";
    	}
    	cout << "Fitness of input vector " << i << " is " << f <<"\n";
		curr_state_0.clear(); curr_state_1.clear(); curr_state_2.clear();
        curr_state_3.clear(); curr_state_4.clear(); states.clear();

        fitness.push_back(f);
    }
  
    return fitness;
}


std::vector< std::vector< std::vector< int > > > compute_using_ga(std::vector< std::vector< std::vector< int > > > curr_population, int num_vec){
	std::vector< std::vector< std::vector< int > > > next_population;
	int num_inputs = circuit -> numpri;
    for(int j = 0; j < curr_population.size()/2; j++){
//////////////////////////////////////////////////////////////////////////////
//              RANDOMLY SELECT TWO PARENTS
//////////////////////////////////////////////////////////////////////////////
    	int random_par1 = (rand())%curr_population.size();
        int random_par2 = (rand())%curr_population.size();

        std::vector< std::vector< int > > parent1 = curr_population[random_par1];
        std::vector< std::vector< int > > parent2 = curr_population[random_par2];

        std::vector< std::vector< int > > child1, child2;
        std::vector<int> temp1, temp2;
///////////////////////////////////////////////////////////////////////////////////////////
//                     SELECTED TWO PARENTS
///////////////////////////////////////////////////////////////////////////////////////////
//                      UNIFORM CROSSOVER
///////////////////////////////////////////////////////////////////////////////////////////
// Uniform crossover( random_par1, random_par2, child1, child2)
        std::vector< std::vector< int > > random_mask = generate_random_test_sequence(num_inputs, num_vec);
        for(int vec =0; vec < num_vec ; vec++){
            for(int in =0; in < num_inputs; in++){
                if(random_mask[vec][in] == 0){
                    temp1.push_back(parent1[vec][in]);
                    temp2.push_back(parent2[vec][in]);
                }else if(random_mask[vec][in] == 1){
                    temp1.push_back(parent2[vec][in]);
                    temp2.push_back(parent1[vec][in]);
                }
            }
            child1.push_back(temp1); child2.push_back(temp2);
            temp1.clear(); temp2.clear();
        }
//////////////////////////////////////////////////////////////////////////////////////////
//              CROSSOVER DONE
///////////////////////////////////////////////////////////////////////////////////////////
//                  MUTATION
///////////////////////////////////////////////////////////////////////////////////////////
//1% mutation of child1 and child2
// To determine 1% mutation we will first calculate how many bits should we flip
// Then we will randomly select said number of bits
    int num_bits_to_flip = 0.01*num_vec*num_inputs;
    
	if(num_bits_to_flip > 0){
    	int range = num_vec*num_inputs;
    	for(int flips =0; flips < num_bits_to_flip; flips++){
        	int bit = rand()%range;
        	child1[bit/num_inputs][bit%num_inputs] = 1 - child1[bit/num_inputs][bit%num_inputs];
        	child2[bit/num_inputs][bit%num_inputs] = 1 - child2[bit/num_inputs][bit%num_inputs];
    	}
    }
/////////////////////////////////////////////////////////////////////////////////////////
//                  MUTATION DONE
/////////////////////////////////////////////////////////////////////////////////////////

	next_population.push_back(child1); next_population.push_back(child2);
	}
	
	return next_population;
}

//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
    ifstream vecFile;
    ofstream vecFile_w;
    string cktName, vecName;
    int totalNumVec, vecWidth, i;
    int nameIndex;
    double ut;
    clock_t start, end;

    start = clock();

    if ((argc != 2) && (argc != 3))
    {
	cerr << "Usage: " << argv[0] << "[-io] <ckt>\n";
	cerr << "The -i option is to begin the circuit in a state in *.initState.\n";
	cerr << "and -o option is to OBSERVE fault-free outputs and FF's.\n";
	cerr << " Example: " << argv[0] << " s27\n";
	cerr << " Example2: " << argv[0] << " -o s27\n";
	cerr << " Example3: " << argv[0] << " -io s27\n";
	exit(-1);
    }

    if (argc == 3)
    {
	i = 1;
	nameIndex = 2;
	while (argv[1][i] != EOS)
	{
	    switch (argv[1][i])
	    {
		case 'o':
		    OBSERVE = 1;
		    break;
		case 'i':
		    INIT0 = 1;
		    break;
		default:
	    	    cerr << "Invalid option: " << argv[1] << "\n";
	    	    cerr << "Usage: " << argv[0] << " [-io] <ckt> <type>\n";
	    	    cerr << "The -i option is to begin the circuit in a state in *.initState.\n";
	    	    cerr << "and o option is to OBSERVE fault-free outputs.\n";
		    exit(-1);
		    break;
	    } 	// switch
	    i++;
	}	// while
    }
    else	// no option
    {
	nameIndex = 1;
	OBSERVE = 0;
	INIT0 = 0;
    }

    cktName = argv[nameIndex];
    vecName = argv[nameIndex];
    vecName += ".vec";
	
	int num_vectors = 10;
    circuit = new gateLevelCkt(cktName);
    
    int num_inputs = circuit -> numpri;
    int num_ff = circuit -> numff;
    
	for( int count = 0; count < num_ff; count++){
        count_zero_ffs.push_back(0);
        count_one_ffs.push_back(0);
    }
	cout<< "\nCONTROLLABILITY BASED PARTITIONING...\n";
	partition(num_vectors, vecName);
	cout << "CONTROLLABILITY BASED PARTITIONING COMPLETE...\n\n";
	
	std::vector< std::vector < std::vector < int> > > test_seq_population, curr_population, next_population;
    for( int times =0; times < 50; times++){
    	curr_population.clear();
    	cout << "\n ****************************************  TIMES  =" << times+1 << "*******************************************\n";
        //Create a current population of 1000 vectors
        for( int a = 0; a<100; a++) curr_population.push_back(generate_random_test_sequence(num_inputs, num_vectors));
        int num_test_seqs = curr_population.size();
        int index= 0; 
		float temp = 0;
        std::vector< float> fitness;
        std::vector< char> states;
        cout<< "Computing Fitness...\n\n";
        fitness = compute_fitness(curr_population, num_vectors, vecName);
        cout<< "Fitness Computed...\n\n";
		for( int gen_num = 0; gen_num < 500 ; gen_num++){
        	cout << "\nGENERATION NUMBER "<< gen_num +1;
            
            cout << "\n\tChoose maximum fitness element...\n";
            index =0; temp =0;
            
            for(int f_i =0; f_i < fitness.size(); f_i++){
            	if(fitness[f_i] > temp){
            		temp = fitness[f_i];
            		index = f_i;
				}	
			}
            test_seq_population.push_back(curr_population[index]);
            
            increment_state_table(vecName, curr_population[index], num_vectors);
            next_population.clear();
            cout << "\n\tCreate next population...\n";
			next_population = compute_using_ga(curr_population, num_vectors);
            fitness.clear();
            curr_population.clear();
            curr_population = next_population;
            cout<< "Compute fitness...\n\n";
            fitness = compute_fitness(curr_population, num_vectors, vecName);
        }
    }
    
    //ALL_VECTORS ARE STORED IN OUTPUT.VEC FILE
	ofstream output ("output.vec");
	for(int i =0; i < test_seq_population.size(); i++){
		output << circuit->numpri;
		output << "\n";
		for(int j=0; j< test_seq_population[i].size(); j++){
			for(int k =0; k<test_seq_population[i][j].size(); k++){
				output << test_seq_population[i][j][k];
			}
			output << "\n";
		}
		output << "END";
		output << "\n"; 
	}
	output.close();
    
	// output results
    end = clock();
    ut = (double) (end-start);
    cout << "Number of vectors: " << totalNumVec << "\n";
    cout << "Number of clock cycles elapsed: " << (double) ut << "\n";
}

////////////////////////////////////////////////////////////////////////
inline void gateLevelCkt::insertEvent(int levelN, int gateN)
{
    levelEvents[levelN][levelLen[levelN]] = gateN;
    levelLen[levelN]++;
}

////////////////////////////////////////////////////////////////////////
// gateLevelCkt class
////////////////////////////////////////////////////////////////////////
// constructor: reads in the *.lev file for the gate-level ckt
gateLevelCkt::gateLevelCkt(string cktName)
{
    ifstream yyin;
    string fName;
    int i, j, count;
    char c;
    int netnum, junk;
    int f1, f2, f3;
    int levelSize[MAXlevels];

    fName = cktName + ".lev";
    yyin.open(fName.c_str(), ios::in);
    if (!yyin)
    {
	cerr << "Can't open .lev file\n";
	exit(-1);
    }

    numpri = numgates = numout = maxlevels = numff = 0;
    maxLevelSize = 32;
    for (i=0; i<MAXlevels; i++)
	levelSize[i] = 0;

    yyin >>  count;	// number of gates
    yyin >> junk;

    // allocate space for gates
    gtype = new unsigned char[count+64];
    fanin = new short[count+64];
    fanout = new short[count+64];
    levelNum = new int[count+64];
    po = new unsigned[count+64];
    inlist = new int * [count+64];
    fnlist = new int * [count+64];
    sched = new char[count+64];
    value1 = new unsigned int[count+64];
    value2 = new unsigned int[count+64];

    // now read in the circuit
    numTieNodes = 0;
    for (i=1; i<count; i++)
    {
	yyin >> netnum;
	yyin >> f1;
	yyin >> f2;
	yyin >> f3;

	numgates++;
	gtype[netnum] = (unsigned char) f1;
	f2 = (int) f2;
	levelNum[netnum] = f2;
	levelSize[f2]++;

	if (f2 >= (maxlevels))
	    maxlevels = f2 + 5;
	if (maxlevels > MAXlevels)
	{
	    cerr << "MAXIMUM level (" << maxlevels << ") exceeded.\n";
	    exit(-1);
	}

	fanin[netnum] = (int) f3;
	if (f3 > MAXFanout)
	    cerr << "Fanin count (" << f3 << " exceeded\n";

	if (gtype[netnum] == T_input)
	{
	    inputs[numpri] = netnum;
	    numpri++;
	}
	if (gtype[netnum] == T_dff)
	{
	    if (numff >= (MAXFFS-1))
	    {
		cerr << "The circuit has more than " << MAXFFS -1 << " FFs\n";
		exit(-1);
	    }
	    ff_list[numff] = netnum;
	    numff++;
	}

	sched[netnum] = 0;

	// now read in the faninlist
	inlist[netnum] = new int[fanin[netnum]];
	for (j=0; j<fanin[netnum]; j++)
	{
	    yyin >> f1;
	    inlist[netnum][j] = (int) f1;
	}

	for (j=0; j<fanin[netnum]; j++)	  // followed by close to samethings
	    yyin >> junk;

	// read in the fanout list
	yyin >> f1;
	fanout[netnum] = (int) f1;

	if (gtype[netnum] == T_output)
	{
	    po[netnum] = TRUE;
	    outputs[numout] = netnum;
	    numout++;
	}
	else
	    po[netnum] = 0;

	if (fanout[netnum] > MAXFanout)
	    cerr << "Fanout count (" << fanout[netnum] << ") exceeded\n";

	fnlist[netnum] = new int[fanout[netnum]];
	for (j=0; j<fanout[netnum]; j++)
	{
	    yyin >> f1;
	    fnlist[netnum][j] = (int) f1;
	}

	if (gtype[netnum] == T_tie1)
        {
	    TIES[numTieNodes] = netnum;
	    numTieNodes++;
	    if (numTieNodes > 511)
	    {
		cerr, "Can't handle more than 512 tied nodes\n";
		exit(-1);
	    }
            value1[netnum] = ALLONES;
            value2[netnum] = ALLONES;
        }
        else if (gtype[netnum] == T_tie0)
        {
	    TIES[numTieNodes] = netnum;
	    numTieNodes++;
	    if (numTieNodes > 511)
	    {
		cerr << "Can't handle more than 512 tied nodes\n";
		exit(-1);
	    }
            value1[netnum] = 0;
            value2[netnum] = 0;
        }
        else
        {
	    // assign all values to unknown
	    value1[netnum] = 0;
	    value2[netnum] = ALLONES;
	}

	// read in and discard the observability values
        yyin >> junk;
        yyin >> c;    // some character here
        yyin >> junk;
        yyin >> junk;

    }	// for (i...)
    yyin.close();
    numgates++;
    numFaultFreeGates = numgates;

    // initialize reset signal masking related variables
    reset_power0 = new int[numpri];
    reset_power1 = new int[numpri];
    prob_maskto0 = new double[numpri];
    prob_maskto1 = new double[numpri];

    // now compute the maximum width of the level
    for (i=0; i<maxlevels; i++)
    {
	if (levelSize[i] > maxLevelSize)
	    maxLevelSize = levelSize[i] + 1;
    }

    // allocate space for the faulty gates
    for (i = numgates; i < numgates+64; i+=2)
    {
        inlist[i] = new int[2];
        fnlist[i] = new int[MAXFanout];
        po[i] = 0;
        fanin[i] = 2;
        inlist[i][0] = i+1;
	sched[i] = 0;
    }

    cout << "Successfully read in circuit:\n";
    cout << "\t" << numpri << " PIs.\n";
    cout << "\t" << numout << " POs.\n";
    cout << "\t" << numff << " Dffs.\n";
    cout << "\t" << numFaultFreeGates << " total number of gates\n";
    cout << "\t" << maxlevels / 5 << " levels in the circuit.\n";
    goodState = new char[numff];
    ffMap = new int[numgates];
    // get the ffMap
    for (i=0; i<numff; i++)
    {
	ffMap[ff_list[i]] = i;
	goodState[i] = 'X';
    }

    setupWheel(maxlevels, maxLevelSize);
    setFaninoutMatrix();

    if (INIT0)	// if start from a initial state
    {
	RESET_FF1 = new unsigned int [numff+2];
	RESET_FF2 = new unsigned int [numff+2];

	fName = cktName + ".initState";
	yyin.open(fName.c_str(), ios::in);
	if (!yyin)
	{	cerr << "Can't open " << fName << "\n";
		exit(-1);}

	for (i=0; i<numff; i++)
	{
	    yyin >>  c;

	    if (c == '0')
	    {
		RESET_FF1[i] = 0;
		RESET_FF2[i] = 0;
	    }
	    else if (c == '1')
	    {
		RESET_FF1[i] = ALLONES;
		RESET_FF2[i] = ALLONES;
	    }
	    else 
	    {
		RESET_FF1[i] = 0;
		RESET_FF2[i] = ALLONES;
	    }
	}

	yyin.close();
    }

    // calculate reset signal masking related values
    FindResetInputs();
}


////////////////////////////////////////////////////////////////////////
////********************** PRIVATE FUNCTIONS ***********************////
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// FindResetInputs()
// Calculate reset signal masking related values
////////////////////////////////////////////////////////////////////////
void gateLevelCkt::FindResetInputs()
{
    char* input_vector;
    input_vector = new char[numpri];
    for(int i = 0; i < numpri; i++)
    {
        input_vector[i] = 'X';
    }
    for(int i = 0; i < numpri; i++)
    {
        // calculate reset_power0
        input_vector[i] = '0';
        applyVector(input_vector);
        goodsim();
        applyVector(input_vector);
        goodsim();
        int reset_ff_count = 0;
        observeFFs();
        for (int j = 0; j < numff; j++)
        {
            if ((value1[ff_list[j]] && value2[ff_list[j]])       // normal 1
                || (!value1[ff_list[j]] && !value2[ff_list[j]])) // normal 0
            {
                reset_ff_count++;
            }
        }
        reset_power0[i] = reset_ff_count;
        reset_ff_count = 0;
        reset(); // reset the circuit to initial state


        //calculate reset_power1
        input_vector[i] = '1';
        applyVector(input_vector);
        goodsim();
        applyVector(input_vector);
        goodsim();
        observeFFs();

        for (int j = 0; j < numff; j++)
        {
            if ((value1[ff_list[j]] && value2[ff_list[j]])       // normal 1
                || (!value1[ff_list[j]] && !value2[ff_list[j]])) // normal 0
            {
                reset_ff_count++;
            }
        }
        reset_power1[i] = reset_ff_count;

        input_vector[i] = 'X';
        reset();

        prob_maskto0[i] = reset_power1[i] / (double)numff;
        prob_maskto1[i] = reset_power0[i] / (double)numff;
        prob_maskto0[i] = sqrt(prob_maskto0[i]);
        prob_maskto1[i] = sqrt(prob_maskto1[i]);

        //debug cout

        cout << "reset_power0[" << i << "] = " << reset_power0[i] << endl;
        cout << "reset_power1[" << i << "] = " << reset_power1[i] << endl;
        cout << "prob_maskto0[" << i << "] = " << prob_maskto0[i] << endl;
        cout << "prob_maskto1[" << i << "] = " << prob_maskto1[i] << endl;

    }
}



////////////////////////////////////////////////////////////////////////
// setFaninoutMatrix()
//	This function builds the matrix of succOfPredOutput and 
// predOfSuccInput.
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::setFaninoutMatrix()
{
    int i, j, k;
    int predecessor, successor;
    int checked[MAXFanout];
    int checkID;	// needed for gates with fanouts to SAME gate
    int prevSucc, found;

    predOfSuccInput = new int *[numgates+64];
    succOfPredOutput = new int *[numgates+64];
    for (i=0; i<MAXFanout; i++)
	checked[i] = 0;
    checkID = 1;

    prevSucc = -1;
    for (i=1; i<numgates; i++)
    {
	predOfSuccInput[i] = new int [fanout[i]];
	succOfPredOutput[i] = new int [fanin[i]];

	for (j=0; j<fanout[i]; j++)
	{
	    if (prevSucc != fnlist[i][j])
		checkID++;
	    prevSucc = fnlist[i][j];

	    successor = fnlist[i][j];
	    k=found=0;
	    while ((k<fanin[successor]) && (!found))
	    {
		if ((inlist[successor][k] == i) && (checked[k] != checkID))
		{
		    predOfSuccInput[i][j] = k;
		    checked[k] = checkID;
		    found = 1;
		}
		k++;
	    }
	}

	for (j=0; j<fanin[i]; j++)
	{
	    if (prevSucc != inlist[i][j])
		checkID++;
	    prevSucc = inlist[i][j];

	    predecessor = inlist[i][j];
	    k=found=0;
	    while ((k<fanout[predecessor]) && (!found))
	    {
		if ((fnlist[predecessor][k] == i) && (checked[k] != checkID))
		{
		    succOfPredOutput[i][j] = k;
		    checked[k] = checkID;
		    found=1;
		}
		k++;
	    }
	}
    }

    for (i=numgates; i<numgates+64; i+=2)
    {
	predOfSuccInput[i] = new int[MAXFanout];
	succOfPredOutput[i] = new int[MAXFanout];
    }
}

////////////////////////////////////////////////////////////////////////
// setTieEvents()
//	This function set up the events for tied nodes
////////////////////////////////////////////////////////////////////////
void gateLevelCkt::setTieEvents()
{
    int predecessor, successor;
    int i, j;

    for (i = 0; i < numTieNodes; i++)
    {
	  // different from previous time frame, place in wheel
	  for (j=0; j<fanout[TIES[i]]; j++)
	  {
	    successor = fnlist[TIES[i]][j];
	    if (sched[successor] == 0)
	    {
	    	insertEvent(levelNum[successor], successor);
		sched[successor] = 1;
	    }
	  }
    }	// for (i...)

    // initialize state if necessary
    if (INIT0 == 1)
    {
cout << "Initialize circuit to values in *.initState!\n";
	for (i=0; i<numff; i++)
	{
	    value1[ff_list[i]] = value2[ff_list[i]] = RESET_FF1[i];

	  for (j=0; j<fanout[ff_list[i]]; j++)
	  {
	    successor = fnlist[ff_list[i]][j];
	    if (sched[successor] == 0)
	    {
	    	insertEvent(levelNum[successor], successor);
		sched[successor] = 1;
	    }
	  }	// for j

	    predecessor = inlist[ff_list[i]][0];
	    value1[predecessor] = value2[predecessor] = RESET_FF1[i];

	  for (j=0; j<fanout[predecessor]; j++)
	  {
	    successor = fnlist[predecessor][j];
	    if (sched[successor] == 0)
	    {
	    	insertEvent(levelNum[successor], successor);
		sched[successor] = 1;
	    }
	  }	// for j

	}	// for i
    }	// if (INIT0)
}

////////////////////////////////////////////////////////////////////////
// applyVector()
//	This function applies the vector to the inputs of the ckt.
////////////////////////////////////////////////////////////////////////
void gateLevelCkt::applyVector(char *vec)
{
	
    unsigned int origVal1, origVal2;
    char origBit;
    int successor;
    int i, j;

    for (i = 0; i < numpri; i++)
    {
	origVal1 = value1[inputs[i]] & 1;
	origVal2 = value2[inputs[i]] & 1;
	if ((origVal1 == 1) && (origVal2 == 1))
	    origBit = '1';
	else if ((origVal1 == 0) && (origVal2 == 0))
	    origBit = '0';
	else
	    origBit = 'x';

	if ((origBit != vec[i]) && ((origBit != 'x') || (vec[i] != 'X')))
	{
	  switch (vec[i])
	  {
	    case '0':
		value1[inputs[i]] = 0;
		value2[inputs[i]] = 0;
		break;
	    case '1':
		value1[inputs[i]] = ALLONES;
		value2[inputs[i]] = ALLONES;
		break;
	    case 'x':
	    case 'X':
		value1[inputs[i]] = 0;
		value2[inputs[i]] = ALLONES;
		break;
	    default:
		cerr << vec[i] << ": error in the input vector.\n";
		exit(-1);
	  }	// switch

	  // different from previous time frame, place in wheel
	  for (j=0; j<fanout[inputs[i]]; j++)
	  {
	    successor = fnlist[inputs[i]][j];
	    if (sched[successor] == 0)
	    {
	    	insertEvent(levelNum[successor], successor);
		sched[successor] = 1;
	    }
	  }
	}	// if ((origBit...)
    }	// for (i...)
}

////////////////////////////////////////////////////////////////////////
// lowWheel class
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::setupWheel(int numLevels, int levelSize)
{
    int i;

    numlevels = numLevels;
    levelLen = new int[numLevels];
    levelEvents = new int * [numLevels];
    for (i=0; i < numLevels; i++)
    {
	levelEvents[i] = new int[levelSize];
	levelLen[i] = 0;
    }
    activation = new int[levelSize];
    
    actFFList = new int[numff + 1];
}

////////////////////////////////////////////////////////////////////////
int gateLevelCkt::retrieveEvent()
{
    while ((levelLen[currLevel] == 0) && (currLevel < maxlevels))
	currLevel++;

    if (currLevel < maxlevels)
    {
    	levelLen[currLevel]--;
        return(levelEvents[currLevel][levelLen[currLevel]]);
    }
    else
	return(-1);
}

////////////////////////////////////////////////////////////////////////
// goodsim() -
//	Logic simulate. (no faults inserted)
////////////////////////////////////////////////////////////////////////
void gateLevelCkt::goodsim()
{
    int sucLevel;
    int gateN, predecessor, successor;
    int *predList;
    int i;
    unsigned int val1, val2, tmpVal;

    currLevel = 0;
    actLen = actFFLen = 0;
    while (currLevel < maxlevels)
    {
    	gateN = retrieveEvent();
	if (gateN != -1)	// if a valid event
	{
	    sched[gateN]= 0;
    	    switch (gtype[gateN])
    	    {
	      case T_and:
    	    	val1 = val2 = ALLONES;
		predList = inlist[gateN];
    	    	for (i=0; i<fanin[gateN]; i++)
    	    	{
		    predecessor = predList[i];
		    val1 &= value1[predecessor];
		    val2 &= value2[predecessor];
    	    	}
	    	break;
	      case T_nand:
    	        val1 = val2 = ALLONES;
		predList = inlist[gateN];
    	    	for (i=0; i<fanin[gateN]; i++)
    	    	{
		    predecessor = predList[i];
		    val1 &= value1[predecessor];
		    val2 &= value2[predecessor];
    	    	}
	    	tmpVal = val1;
	    	val1 = ALLONES ^ val2;
	    	val2 = ALLONES ^ tmpVal;
	    	break;
	      case T_or:
    	        val1 = val2 = 0;
		predList = inlist[gateN];
    	        for (i=0; i<fanin[gateN]; i++)
    	        {
		    predecessor = predList[i];
		    val1 |= value1[predecessor];
		    val2 |= value2[predecessor];
    	    	}
	    	break;
	      case T_nor:
    	    	val1 = val2 = 0;
		predList = inlist[gateN];
		for (i=0; i<fanin[gateN]; i++)
    	    	{
		    predecessor = predList[i];
		    val1 |= value1[predecessor];
		    val2 |= value2[predecessor];
    	    	}
	    	tmpVal = val1;
	    	val1 = ALLONES ^ val2;
	    	val2 = ALLONES ^ tmpVal;
	    	break;
	      case T_not:
	    	predecessor = inlist[gateN][0];
	    	val1 = ALLONES ^ value2[predecessor];
	    	val2 = ALLONES ^ value1[predecessor];
	    	break;
	      case T_buf:
	    	predecessor = inlist[gateN][0];
	    	val1 = value1[predecessor];
	    	val2 = value2[predecessor];
		break;
	      case T_dff:
	    	predecessor = inlist[gateN][0];
	    	val1 = value1[predecessor];
	    	val2 = value2[predecessor];
		actFFList[actFFLen] = gateN;
		actFFLen++;
	    	break;
	      case T_xor:
		predList = inlist[gateN];
	    	val1 = value1[predList[0]];
	    	val2 = value2[predList[0]];

            	for(i=1;i<fanin[gateN];i++)
            	{
	    	    predecessor = predList[i];
                    tmpVal = ALLONES^(((ALLONES^value1[predecessor]) &
              		(ALLONES^val1)) | (value2[predecessor]&val2));
                    val2 = ((ALLONES^value1[predecessor]) & val2) |
                  	(value2[predecessor] & (ALLONES^val1));
                    val1 = tmpVal;
            	}
	    	break;
	      case T_xnor:
		predList = inlist[gateN];
		val1 = value1[predList[0]];
	    	val2 = value2[predList[0]];

            	for(i=1;i<fanin[gateN];i++)
            	{
	    	    predecessor = predList[i];
                    tmpVal = ALLONES^(((ALLONES^value1[predecessor]) &
                   	(ALLONES^val1)) | (value2[predecessor]&val2));
                    val2 = ((ALLONES^value1[predecessor]) & val2) |
                  	(value2[predecessor]& (ALLONES^val1));
                    val1 = tmpVal;
            	}
	    	tmpVal = val1;
		val1 = ALLONES ^ val2;
	    	val2 = ALLONES ^ tmpVal;
	    	break;
	      case T_output:
		predecessor = inlist[gateN][0];
	    	val1 = value1[predecessor];
	    	val2 = value2[predecessor];
	        break;
	      case T_input:
	      case T_tie0:
	      case T_tie1:
	      case T_tieX:
	      case T_tieZ:
	    	val1 = value1[gateN];
	    	val2 = value2[gateN];
	    	break;
	      default:
		cerr << "illegal gate type1 " << gateN << " " << gtype[gateN] << "\n";
		exit(-1);
    	    }	// switch

	    // if gate value changed
    	    if ((val1 != value1[gateN]) || (val2 != value2[gateN]))
	    {
		value1[gateN] = val1;
		value2[gateN] = val2;

		for (i=0; i<fanout[gateN]; i++)
		{
		    successor = fnlist[gateN][i];
		    sucLevel = levelNum[successor];
		    if (sched[successor] == 0)
		    {
		      if (sucLevel != 0)
			insertEvent(sucLevel, successor);
		      else	// same level, wrap around for next time
		      {
			activation[actLen] = successor;
			actLen++;
		      }
		      sched[successor] = 1;
		    }
		}	// for (i...)
	    }	// if (val1..)
	}	// if (gateN...)
    }	// while (currLevel...)
    // now re-insert the activation list for the FF's
    for (i=0; i < actLen; i++)
    {
	insertEvent(0, activation[i]);
	sched[activation[i]] = 0;

        predecessor = inlist[activation[i]][0];
        gateN = ffMap[activation[i]];
        if (value1[predecessor])
            goodState[gateN] = '1';
        else if (value2[predecessor] == 0)
            goodState[gateN] = '0';
        else
            goodState[gateN] = 'X';
    }
}

////////////////////////////////////////////////////////////////////////
// observeOutputs()
//	This function prints the outputs of the fault-free circuit.
////////////////////////////////////////////////////////////////////////
unsigned int* gateLevelCkt::getValue1FFs()
{
    unsigned int* value1_ffs;
    value1_ffs = new unsigned int[numff];
    for(int i =0; i < numff; i++){
    	*(value1_ffs + i) = value1[ff_list[i]];
	}
	return value1_ffs;
}

void gateLevelCkt::reset()
{
    // reset gate values
    //id_unknown_max = 0;
    for (int i = 0; i < numgates; i++)
    {
        value1[i] = 0;
        value2[i] = ALLONES;
        //id_unknown[i] = 0;
    }

    // reset stored ff values
    /*
    for (int i = 0; i < numff; i++)
    {
        value1[ff_list[i]] = 0;
        value2[ff_list[i]] = ALLONES;
        //id_unknown_ffs[i] = 0;
    }
     */
}

unsigned int* gateLevelCkt::getValue2FFs()
{
    unsigned int* value2_ffs;
    value2_ffs = new unsigned int[numff];
    for(int i =0; i < numff; i++){
    	*(value2_ffs + i) = value2[ff_list[i]];
	}
	return value2_ffs;
}

////////////////////////////////////////////////////////////////////////
// observeOutputs()
//	This function prints the outputs of the fault-free circuit.
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::observeOutputs()
{
    int i;

    cout << "Outs:\t";
    for (i=0; i<numout; i++)
    {
        if (value1[outputs[i]] && value2[outputs[i]])
            cout << "1";
        else if ((value1[outputs[i]] == 0) && (value2[outputs[i]] == 0))
            cout << "0";
        else
            cout << "X";
    }
    cout << endl;
}


////////////////////////////////////////////////////////////////////////
// observeFFs()
// prints the stored value of FFs
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::observeFFs()
{
    int i;
    cout << "FFs:\t";
    for (i=0; i<numff; i++)
    {
        if (value1[ff_list[i]] && value2[ff_list[i]])
            cout << "1";
        else if (!value1[ff_list[i]] && !value2[ff_list[i]])
            cout << "0";
        else
            cout << "X";
    }
    cout << endl;

}

////////////////////////////////////////////////////////////////////////
// getPorbMaskto0, getProbMaskto1
// return the reset signal masking related values
////////////////////////////////////////////////////////////////////////

double* gateLevelCkt::getProbMaskto0()
{
    return prob_maskto0;
}

double* gateLevelCkt::getProbMaskto1()
{
    return prob_maskto1;
}


