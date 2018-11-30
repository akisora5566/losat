////////////////////////////////////////////////////////////////////////
// File: main.cpp
// Author: Anikate Kaw, Victor Gan
// Acknowledgement: Dr. Michiael Shiao @ECE,Virginia Tech
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*****************************INSTRUCTIONS****************************/
/*
I. OVERVIEW
This program is based on the lsim_3val program written by Dr. Shiao.
It can evaluate 3-value sequential circuit with distinguished Xs.

II. USER INSTRUCTION
This program is designed to compile and execute on linux.
We recommend to use cmake to build the executable file
$ cd path/contain/CMakeList.txt
$ cmake --build desired/path/to/target --target losat -- -j 2
$ ./path/to/target/losat -o <circuit name>

for example:
$ cmake --build cmake-build-debug --target losat -- -j 2
$ ./cmake-build-debug/losat -o s382

<<<<<<<<<< TODO: Description of the main function

III. EXPLONATION OF FUNCTIONS & ALGORITHM
The following are the functions built by Dr. Shiao. They consist of almost
60% of the functions.
    a) data structure
    b) circuit read in
    c) event wheel

<<<<<<<<<< TODO: Elaberate the functions we added to this program


IV.CODE BY STUDENTS
<<<<<<<<<< TODO: Code by students
1) private variables
    int *id_unknown;
    int id_unknown_max;
2) public functions
    void resetIDunknown()
3) major modified functions
    constructor (initializing)
    applyVector
    goodsim (critical function changes in gates)

*/
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <sys/time.h>
#include <string>
#include <ctime>
#include <map>
#include "gateLevelCkt.h"
using namespace std;


//////////////////////////////////////////////////////////////////////////
// Global Variables
char vectr[5120];
gateLevelCkt *circuit;
int OBSERVE, INIT0;
int vecNum=0;
int numTieNodes;
int TIES[512];
char* input_vectors[10000];
std::vector<int> count_zero_ffs, count_one_ffs, ff_bias, state_partition;
std::map< std::vector<char>, int > state_table0, state_table1, state_table2, state_table3, state_table4;

//////////////////////////////////////////////////////////////////////////
// Functions start here
//////////////////////////////////////////////////////////////////////////

std::vector<std::vector<int> > generate_random_test_sequence(int num_inputs, int num_vec){
    std::vector<std::vector<int> > test_sequence;
    std::vector<int> temp;
    int a;
    double probTo0, probTo1, val;
    double* probMaskTo0 = circuit->getProbMaskto0();
    double* probMaskTo1 = circuit->getProbMaskto1();
	for (int i =0; i < num_vec; i++){
        for (int j=0; j < num_inputs; j++){
        	//    RESET SIGNAL MASKING////////////////////////////////
        	/*probTo0 = *(probMaskTo0 + j);
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
			}*/
			/////////////////////////////////////
			 a = rand()%2;     //-------------------WITH NO RESET SIGNAL MASKING------------
			 //////////////////////////////////////////////////////////
            temp.push_back(a);
        }
        test_sequence.push_back(temp);
        temp.clear();
    }
    return test_sequence;
}

char** generate_test_vector(int num_inputs, int num_vectors){
    char** test_vec;
    int a;
    double probTo0, probTo1, val;
    double* probMaskTo0 = circuit->getProbMaskto0();
    double* probMaskTo1 = circuit->getProbMaskto1();
    test_vec = new char* [num_vectors];
    for (int i = 0; i < num_vectors; i++)
    {
        test_vec[i] = new char[num_inputs];
        for (int j = 0; j < num_inputs; j++)
        {
        	//-----------------------RESET SIGNAL MASKING-----------------------
        	/*
        	probTo0 = *(probMaskTo0 + j);
            probTo1 = *(probMaskTo1 + j);
            if(probTo0 > 0.4 || probTo1 > 0.4){
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
			/////////////////////////////////////////////////////////////////////
            a = rand()%2;
            /////////////////////////////////////////////////////////////////////
            if(a == 0) test_vec[i][j] = '0'; else if( a ==1) test_vec[i][j] = '1';
        }
    }
    return test_vec;
}

std::vector< float> compute_fitness(std::vector< std::vector< std::vector <int> > > curr_population, int num_vectors){
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
        for(int elem =0; elem < num_vectors; elem++){
        	test_vec1[elem] = new char[num_inputs];
            for(int elem_i = 0; elem_i < num_inputs; elem_i++){
                if(curr_population[i][elem][elem_i] == 0)  test_vec1[elem][elem_i] = '0'; else if( curr_population[i][elem][elem_i] ==1)test_vec1[elem][elem_i] = '1';
            }
        }
        circuit->reset();
		circuit->setTieEvents();
        circuit->LogicSim(test_vec1, num_vectors);
		value1ffs = circuit->getValue1FFs();
        value2ffs = circuit->getValue2FFs();
		//circuit->reset();
		//cout<<"\n***************** ELEMENT " << i<< "**************************";
		//cout << "Final state of flip flop is ";
		for(int ff=0; ff < num_ff; ff++){
			if(*value1ffs == *value2ffs){
        		if(*value1ffs == 0){
        			states.push_back('0');//cout <<"0";
				}else{
        		states.push_back('1');//cout <<"1";
				} 
			}else{
				states.push_back('X');//cout << "X";
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
            	f += 0.2*(1/(float)it->second);
            	it->second++;
        	}else{
        		state_table0[curr_state_0] = 1;
        	}

        	if(state_table1.count(curr_state_1)>0){
        		it = state_table1.find(curr_state_1);
            	f += 0.4*(1/(float)it->second);
            	it->second++;
        	}else{
            	state_table1[curr_state_1] = 1;
        	}

        	if(state_table2.count(curr_state_2)>0){
        		it = state_table2.find(curr_state_2);
            	f += 0.6*(1/(float)it->second);
            	it->second++;
        	}else{
            	state_table2[curr_state_2] = 1;
        	}

        	if(state_table3.count(curr_state_3)>0){
        		it = state_table3.find(curr_state_3);
            	f += 0.8*(1/(float)it->second);
            	it->second++;
        	}else{
            	state_table3[curr_state_3] = 1;
        	}

        	if(state_table4.count(curr_state_4)>0){
        		it = state_table4.find(curr_state_4);
            	f += (1/(float)it->second);
            	it->second++;
        	}else{
            	state_table4[curr_state_4] = 1;
        	}
        
        	//cout << i << "after compute fitness\t";
        	cout << "Fitness of input vector " << i << " is " << f <<"\n";
    	}
    	//cout << "Fitness of input vector " << i << " is " << f <<"\n";
		curr_state_0.clear(); curr_state_1.clear(); curr_state_2.clear();
        curr_state_3.clear(); curr_state_4.clear(); states.clear();

        fitness.push_back(f);
    }
  
    return fitness;
}


std::vector< std::vector< std::vector< int > > > compute_using_ga(std::vector< std::vector< std::vector< int > > > curr_population, int num_vec){
	std::vector< std::vector< std::vector< int > > > next_population;
	int num_inputs = circuit -> numpri;
    for(int j = 1; j < curr_population.size()/2; j++){
//////////////////////////////////////////////////////////////////////////////
//              RANDOMLY SELECT TWO PARENTS
//////////////////////////////////////////////////////////////////////////////
    	int random_par1 = (curr_population.size() + rand())%curr_population.size();
        int random_par2 = (curr_population.size() + rand())%curr_population.size();

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


void partition(int num_vectors){
	unsigned int* value1ffs;
    unsigned int* value2ffs;
    char** test_vec;
	for(int count =0; count < 500; count++){
    	test_vec = generate_test_vector(circuit->numpri, num_vectors);
        //start = clock();
        circuit->setTieEvents();
		//cout << "Logic Sim "<< count <<"\n";
		circuit->reset();
        circuit->LogicSim(test_vec, num_vectors);
        value1ffs = circuit->getValue1FFs();
        value2ffs = circuit->getValue2FFs();
		//circuit->reset();

        for(int i=0; i < circuit->numff; i++){
            if(value1ffs[i] == 0 && value2ffs[i] == 0){
                count_zero_ffs[i]++;
            }else if(value1ffs[i] == 1 && value2ffs[i] == 1){
                count_one_ffs[i]++;
            }
        }
    }	

	for(int i = 0; i < circuit->numff; i++){
        ff_bias.push_back(abs(count_one_ffs[i] - count_zero_ffs[i])/500);

        if(ff_bias[i] <= 0.2){
            state_partition.push_back(0);
        }else if(ff_bias[i] <= 0.4){
            state_partition.push_back(1);
        }else if(ff_bias[i] <= 0.6){
            state_partition.push_back(2);
        }else if(ff_bias[i] <= 0.8){
            state_partition.push_back(3);
        }else if(ff_bias[i] <= 1){
            state_partition.push_back(4);
        }
    }	
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

    vectr[0] = toupper(thisChar); //in case input vec would contain 'x'
    if (vectr[0] == 'E')
        return (0);

    for (i=1; i<vecSize; i++)
    {
        inFile >> thisChar;
        vectr[i] = toupper(thisChar); //in case input vec would contain 'x'
    }
    vectr[i] = EOS;

    return(1);
}

// logicSimFromFile() - logic simulates all vectors in vecFile on the circuit
int logicSimFromFile(ifstream &vecFile, int vecWidth)
{
    int moreVec;
    moreVec = 1;
    vecNum = 0;
    while (moreVec)
    {
        moreVec = getVector(vecFile, vecWidth);
        if (moreVec == 1)
        {
            input_vectors[vecNum] = new char[vecWidth];
            for (int j = 0; j < vecWidth; j++)
            {
                input_vectors[vecNum][j] = vectr[j];
            }
            cout << "vector #" << vecNum << ": " << input_vectors[vecNum] << "\n";
            vecNum++;
        }
    }

    char** input_vectors_truesize;
    input_vectors_truesize = new char* [vecNum];
    cout << "vecNum = " << vecNum << endl;
    for (int i = 0; i < vecNum; i++)
    {
        input_vectors_truesize[i] = new char[vecWidth];
        for (int j = 0; j < vecWidth; j++)
        {
            input_vectors_truesize[i][j] = input_vectors[i][j];
        }
    }
    circuit->LogicSim(input_vectors_truesize, vecNum);

    return (vecNum);
}


// main()
int main(int argc, char *argv[])
{

    ifstream vecFile;
    string cktName, vecName;
    int totalNumVec, vecWidth, i, num_inputs, num_ff;
    unsigned int* value1ffs;
    unsigned int* value2ffs;
    int nameIndex;
    double ut;
    char** test_vec;
    clock_t start, end;

    int num_vectors = 10;

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

    circuit = new gateLevelCkt(cktName, INIT0);
    num_inputs = circuit -> numpri;
    num_ff = circuit -> numff;
	double* prob0; double* prob1;
	prob0 = circuit->getProbMaskto0();
    prob1 = circuit->getProbMaskto1();
    
    cout << "\nProbability Mask to 0:\t";
    for(int a =0; a < num_ff; a++){
    	cout << *prob0 << "\t";
    	prob0++;
	}
	
	
    cout << "\nProbability Mask to 1:\t";
    for(int a =0; a < num_ff; a++){
    	cout << *prob1 << "\t";
    	prob1++;
	}
	
	cout << "Press to continue";
	int continue_num;
	cin >> continue_num;
	
	for( int count = 0; count < num_ff; count++){
        count_zero_ffs.push_back(0);
        count_one_ffs.push_back(0);
    }
	cout<< "\nCONTROLLABILITY BASED PARTITIONING...\n";
	partition(num_vectors);
	cout << "CONTROLLABILITY BASED PARTITIONING COMPLETE...\n\n";
    std::vector< std::vector < std::vector < int> > > test_seq_population, curr_population, next_population;
    for( int times =0; times < 10; times++){
    	cout << "\n ****************************************  TIMES  =" << times+1 << "*******************************************\n";
        //Create a current population of 1000 vectors
        for( int a = 0; a<500; a++) curr_population.push_back(generate_random_test_sequence(num_inputs, num_vectors));
        int num_test_seqs = curr_population.size();
        std::vector< float> fitness;
        std::vector< char> states;
        cout<< "Computing Fitness...\n\n";
        fitness = compute_fitness(curr_population, num_vectors);
        cout<< "Fitness Computed...\n\n";
		for( int gen_num = 0; gen_num < 50 ; gen_num++){
        	cout << "\nGENERATION NUMBER "<< gen_num +1;
            //Getting test sequence with best fitness
            cout << "\n\tChoose maximum fitness element...\n";
            int index = std::max_element(fitness.begin(),fitness.end()) - fitness.begin();
            test_seq_population.push_back(curr_population[index]);
            next_population.clear();
            cout << "\n\tCreate next population...\n";
			next_population = compute_using_ga(curr_population, num_vectors);
            fitness.clear();
            curr_population = next_population;
            cout<< "Compute fitness...\n\n";
            fitness = compute_fitness(curr_population, num_vectors);
        }
    }
	
	
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
    /**********************OBSERVE OUTPUTS**********************************
    if (OBSERVE == 1)
    {
        circuit->observeOutputs();
    }
  //  vecFile.close();
	************************************************************************/
    // output results
    //end = clock();
    ut = (double) (end-start);
    cout << "Number of vectors: " << totalNumVec << "\n";
    cout << "Number of clock cycles elapsed: " << ut << "\n";
}
