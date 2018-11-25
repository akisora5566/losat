////////////////////////////////////////////////////////////////////////
// File: project1_906207701.cpp
// student: Victor Gan
// stu ID: 9062-07701
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*****************************INSTRUCTIONS****************************/
/*
I. OVERVIEW
This program is based on lsim_3val written by Dr. Shiao.
It can evaluate 3-value combinational circuit with distinguished Xs.

II. USER INSTRUCTION
This program is designed to compile and execute on linux.
$ g++ -o project1_906207701 project1_906207701.cpp
$ ./project1_906207701 -o <circuit name>

for example:
$ ls
c499a.lec   c499a.vec   project1_906207701  project1_906207701.cpp
$ ./project1_906207701 -o c499a

The program will read in the circuit described in XX.vec and the input
vectors in XX.lec. Then print the output regarding to each input vectors.
Also, when any unknown from the same source squashed at specific gate,
it will print "Xs from PI squashed at gate <gateID>!" on the screen.


III. EXPLONATION OF FUNCTIONS & ALGORITHM
The following are the functions built by Dr. Shiao. They consist of almost
80% of the functions.
    a) data structure
    b) circuit read in
    c) input vector read in
    d) event wheel

This main difference between this and the previous version is that now we
support distinguished unknowns. For example, in previous version, if two
unknowns meet at a 2-input AND gate, they would generate X. In this version,
if the two input are from the same source and complement each other, it
generate 0 instead of X.

To meet the requirement, I added a new variable: id_unknown into the circuit.
Now we have value1, value2 and id_unknown to represent an output value of a gate.
{value1, value2, id_unknown}
output 0: {0, 0, 0}
output 1: {1, 1, 0}
output X: {0, 1, id}, {1, 0, id}
As you can see, I used 0,1 and 1,0 to represent two different unknowns with the same id.


IV.CODE BY STUDENT
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



std::vector<std::vector<int> > generate_random_test_sequence(int num_inputs, int num_vec){
    std::vector<std::vector<int> > test_sequence;
    std::vector<int> temp;
    for (int i =0; i < num_vec; i++){
        for (int j=0; j < num_inputs; j++){
            int a = rand()%2;
            temp.push_back(a);
        }
        test_sequence.push_back(temp);
        temp.clear();
    }
    return test_sequence;
}


////////////////////////////////////////////////////////////////////////
// gateLevelCkt class
////////////////////////////////////////////////////////////////////////

class gateLevelCkt
{
private:
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

    // value of gates
    unsigned int *value1;   // each gate has two output values
    unsigned int *value2;	// 00==0, 01==unknown1, 10==unknown2, 11==1
    int *id_unknown;        // if id_unknown==0 => (output == 0 or 1)
                            // id_unknown && (value1,value2)==(10, 01)
                            // to determine different unknown sources
                            // and their complements
    unsigned int *value1_ffs; // additional value array for FFs
    unsigned int *value2_ffs; // same rule as value1 and value2
    int *id_unknown_ffs; // additional id_unknown array for FFs


    int id_unknown_max;     // record the max id_unknown throughout the circuit

    int **predOfSuccInput;  // predecessor of successor input-pin list
    int **succOfPredOutput; // successor of predecessor output-pin list

    // for simulator
    int **levelEvents;	// event list for each level in the circuit
    int *levelLen;	// evenlist length
    int numlevels;	// total number of levels in wheel
    int currLevel;	// current level
    int *activation;	// activation list for the current level in circuit
    int actLen;		// length of the activation list
    int *actFFList;	// activation list for the FF's
    int actFFLen;	// length of the actFFList
    int smallest_level; // the smallest level which has unevaluated events

public:
    int numff;		// number of FF's
    int numpri;		// number of PIs
    unsigned int *RESET_FF1;	// value of reset ffs read from *.initState
    unsigned int *RESET_FF2;	// value of reset ffs read from *.initState

    gateLevelCkt(string);	// constructor
    void setFaninoutMatrix();	// builds the fanin-out map matrix

    void applyVector(char *);	// apply input vector
    void applyFF(); // apply FFs from the last time frame
    void resetIDunknown(); // reset the value of gates to initial state
    void reset(); // reset the circuit to initial value

    // simulator information
    void setupWheel(int, int);
    void insertEvent(int, int);
    int retrieveEvent();
    void goodsim();		// logic sim (no faults inserted)
    void LogicSim(char** input_vectors, int timeframes); // goodsim the circuit with input vectors

    void setTieEvents();	// inject events from tied nodes

    void observeOutputs();	// print the fault-free outputs
    void observeFFs(); // print the newly reached state
    void printGoodSig(ofstream, int);	// print the fault-free outputs to *.sig
    void printGateValue(int gateN); // print the value of gateN

    char *goodState;		// good state (without scan)

    // get value
    unsigned int* getValue1FFs(); // return the array of value1_ffs
    unsigned int* getValue2FFs(); // return the array of value2_ffs
    int* getIDUnknownFFs(); // reutrn the array of id_unknown_ffs
};

//////////////////////////////////////////////////////////////////////////
// Global Variables

char** generate_test_vector(int num_inputs, int num_vectors){
    char** test_vec;
    test_vec = new char* [num_vectors];
    for (int i = 0; i < num_vectors; i++)
    {
        test_vec[i] = new char[num_inputs];
        for (int j = 0; j < num_inputs; j++)
        {
            int a = rand()%2;
            if(a == 0)  test_vec[i][j] = '0'; else if( a ==1)test_vec[i][j] = '1';
        }
    }
    return test_vec;
}


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

std::vector< int> compute_fitness(std::vector< std::vector< std::vector <int> > > curr_population, int num_vectors){
	std::vector< int> fitness;
	int num_inputs = circuit -> numpri;
	int num_ff = circuit -> numff;
    std::vector< char> states;
    unsigned int* value1ffs; unsigned int* value2ffs;
	std::vector< char> curr_state_0, curr_state_1, curr_state_2, curr_state_3, curr_state_4;
        char** test_vec1;
        test_vec1 = new char* [num_vectors];
        for( int i =0; i < curr_population.size() ; i++){
        	for(int elem =0; elem < num_vectors/*curr_population[i].size()*/; elem++){
        		test_vec1[elem] = new char[num_inputs];
                for(int elem_i = 0; elem_i < num_inputs/*curr_population[i][elem].size()*/; elem_i++){
                    if(curr_population[i][elem][elem_i] == 0)  test_vec1[elem][elem_i] = '0'; else if( curr_population[i][elem][elem_i] ==1)test_vec1[elem][elem_i] = '1';
                }
            }
        	circuit->setTieEvents();
        	circuit->LogicSim(test_vec1, num_vectors);
            value1ffs = circuit->getValue1FFs();
            value2ffs = circuit->getValue2FFs();
			circuit->reset();
			
            for(int i=0; i < num_ff; i++){
                if(value1ffs[i] == 0 && value2ffs[i] == 0){
                        states.push_back(0);
                    }else if(value1ffs[i] == 1 && value2ffs[i] == 1){
                        states.push_back('1');
                    }else if(value1ffs[i] == 0 && value2ffs[i] == 0){
                        states.push_back('0');
                    }else{
                        states.push_back('X');
                    }
                }
                //std::vector<int> states = circuit->logicsim(curr_population[i]);      //Flip flop states
                int f = 0; //compute_fit( state_partition, states);          //compute fitness
                for(int j =0; j < states.size(); j++){
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
                    f += 0.2*(1/it->second);
                    it->second++;
                }else{
                    state_table0[curr_state_0] = 1;
                }

                if(state_table1.count(curr_state_1)>0){
                    it = state_table1.find(curr_state_1);
                    f += 0.4*(1/it->second);
                    it->second++;
                }else{
                    state_table1[curr_state_1] = 1;
                }

                if(state_table2.count(curr_state_2)>0){
                    it = state_table2.find(curr_state_2);
                    f += 0.6*(1/it->second);
                    it->second++;
                }else{
                    state_table2[curr_state_2] = 1;
                }

                if(state_table3.count(curr_state_3)>0){
                    it = state_table3.find(curr_state_3);
                    f += 0.8*(1/it->second);
                    it->second++;
                }else{
                    state_table3[curr_state_3] = 1;
                }

                if(state_table4.count(curr_state_4)>0){
                    it = state_table4.find(curr_state_4);
                    f += (1/it->second);
                    it->second++;
                }else{
                    state_table4[curr_state_4] = 1;
                }
                curr_state_0.clear(); curr_state_1.clear(); curr_state_2.clear();
                curr_state_3.clear(); curr_state_4.clear();

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
    if(num_bits_to_flip == 0) num_bits_to_flip++;

    int range = num_vec*num_inputs;
    for(int flips =0; flips < num_bits_to_flip; flips++){
        int bit = rand()%range;
        child1[bit/num_inputs][bit%num_inputs] = 1 - child1[bit/num_inputs][bit%num_inputs];
        child2[bit/num_inputs][bit%num_inputs] = 1 - child2[bit/num_inputs][bit%num_inputs];
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
        //circuit->setTieEvents();
		cout << "Logic Sim "<< count <<"/n";
        circuit->LogicSim(test_vec, num_vectors);
        value1ffs = circuit->getValue1FFs();
        value2ffs = circuit->getValue2FFs();
		circuit->reset();

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

    circuit = new gateLevelCkt(cktName);
    num_inputs = circuit -> numpri;
    num_ff = circuit -> numff;

    for( int count = 0; count < num_ff; count++){
        count_zero_ffs.push_back(0);
        count_one_ffs.push_back(0);
    }
	cout<< "\nCONTROLLABILITY BASED PARTITIONING...\n\n";
	partition(num_vectors);
	cout << "\nCONTROLLABILITY BASED PARTITIONING COMPLETE...\n\n";
    std::vector< std::vector < std::vector < int> > > test_seq_population, curr_population, next_population;
    for( int times =0; times < 10; times++){
    	cout << "\n ****************************************  TIMES  =" << times << "*******************************************\n";
        //Create a current population of 1000 vectors
        for( int a = 0; a<10; a++) curr_population.push_back(generate_random_test_sequence(num_inputs, num_vectors));
        int num_test_seqs = curr_population.size();
        std::vector< int> fitness;
        std::vector< char> states;
        cout<< "Compute fitness...\n\n";
        fitness = compute_fitness(curr_population, num_vectors);
        for( int gen_num = 0; gen_num < 10 ; gen_num++){
        	cout << "\nGENERATION NUMBER "<< gen_num;
            //Getting test sequence with best fitness
            cout<< "\nChoose maximum fitness element...\n";
            int index = std::max_element(fitness.begin(),fitness.end()) - fitness.begin();
            test_seq_population.push_back(curr_population[index]);
            next_population.clear();
            cout << "\nCreate next population...\n";
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
		for(int j=0; j< test_seq_population[i].size(); j++){
			for(int k =0; k<test_seq_population[i][j].size(); k++){
				output << test_seq_population[i][j][k];
			}
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


////////////////////////////////////////////////////////////////////////
inline void gateLevelCkt::insertEvent(int levelN, int gateN)
{
    levelEvents[levelN][levelLen[levelN]] = gateN;
    levelLen[levelN]++;
    if (levelN < smallest_level)
    {
        smallest_level = levelN;
    }
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
    for (i=0; i<MAXlevels; i++) {
        levelSize[i] = 0;
    }
    yyin >> count;	// number of gates
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
    id_unknown = new int[count+64];
    value1_ffs = new unsigned int[count+64];
    value2_ffs = new unsigned int[count+64];
    id_unknown_ffs = new int[count+64];
    id_unknown_max = 0;
    smallest_level = 0;

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

        // initializing gate values
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
            id_unknown[netnum] = 0;
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
            id_unknown[netnum] = 0;
        }
        else
        {
            // assign all values to unknown
            value1[netnum] = 0;
            value2[netnum] = ALLONES;
            id_unknown[netnum] = 0;
            id_unknown_ffs[netnum] = 0;
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

    // initializing stored ff values
    for (i=0; i<numff; i++)
    {
        value1_ffs[i] = 0;
        value2_ffs[i] = ALLONES;
        id_unknown_ffs[i] = 0;
    }

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
}

////////////////////////////////////////////////////////////////////////
// setFaninoutMatrix()
// This function builds the matrix of succOfPredOutput and
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
            k = found = 0;
            while ((k < fanin[successor]) && (!found))
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
        if (id_unknown[inputs[i]] == 0) //is 1/0 or the initial unknown
        {
            if ((origVal1 == 1) && (origVal2 == 1))
                origBit = '1';
            else if ((origVal1 == 0) && (origVal2 == 0))
                origBit = '0';
            else
                origBit = 'X'; // the initial unknown
        }
        else
        {
            origBit = 'X';
        }

        // if new input is X, we must reapply
        // or if new input is 0 or 1, we reapply only if it changed
        if ((vec[i] == 'X') || (origBit != vec[i]))
        {
            switch (vec[i])
            {
            case '0':
                    value1[inputs[i]] = 0;
                    value2[inputs[i]] = 0;
                    id_unknown[inputs[i]] = 0;
                    //debuging cout
                    //cout << "vec[" << i << "] = " << vec[i] << ", id = " << id_unknown[inputs[i]] << endl;
                    break;
            case '1':
                    value1[inputs[i]] = ALLONES;
                    value2[inputs[i]] = ALLONES;
                    id_unknown[inputs[i]] = 0;
                    //debuging cout
                    //cout << "vec[" << i << "] = " << vec[i] << ", id = " << id_unknown[inputs[i]] << endl;
                    break;
            case 'X':
                    value1[inputs[i]] = 0;
                    value2[inputs[i]] = ALLONES;
                    id_unknown[inputs[i]] = id_unknown_max;
                    id_unknown_max++;
                    //debuging cout
                    //cout << "vec[" << i << "] = " << vec[i] << ", id = " << id_unknown[inputs[i]] << endl;
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
// applyFF()
// Apply FFs from the last time frame
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::applyFF()
{
    unsigned int original_value1, original_value2;
    int original_id_unknown;
    int current_ff = 0;
    int successor = 0;
    /****you already have the ff_list and Map****/
    for (int i = 0; i < numff; i++)
    {
        original_value1 = value1_ffs[i];
        original_value2 = value2_ffs[i];
        original_id_unknown = id_unknown_ffs[i];

        // if not equal to reset value
        if (original_value1 || !original_value2 || original_id_unknown)
        {
            // set value of ffs
            current_ff = ff_list[i];
            value1[current_ff] = value1_ffs[i];
            value2[current_ff] = value2_ffs[i];
            id_unknown[current_ff] = id_unknown_ffs[i];
            //debuging cout
            /*
            char char_val1 = value1[current_ff] ? '1':'0';
            char char_val2 = value2[current_ff] ? '1':'0';
            cout << "ff[" << i << "] = " << char_val1 << " " << char_val2 << ", id = " << id_unknown[current_ff] << endl;
            */

            //put sucessors on the event wheel
            for (int j = 0; j < fanout[current_ff]; j++)
            {
                successor = fnlist[current_ff][j];
                if (sched[successor] == 0)
                {
                    insertEvent(levelNum[successor], successor);
                    sched[successor] = 1;
                }
            }
        }
        //dubuging cout
        /*
        else
        {
            cout << "ff[" << i << "] has no new value." << endl;
        }
        */
    }
}


////////////////////////////////////////////////////////////////////////
// resetIDunknown()
//This function reset the id_unknowns to initial state
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::resetIDunknown()
{
    for (int i = 0; i < numgates; i++)
    {
        id_unknown[i] = 0;
    }
    id_unknown_max = 0;
}


////////////////////////////////////////////////////////////////////////
// reset()
//This function reset the value of all gates to initial state
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::reset()
{
    // reset gate values
    id_unknown_max = 0;
    for (int i = 0; i < numgates; i++)
    {
        value1[i] = 0;
        value2[i] = ALLONES;
        id_unknown[i] = 0;
    }

    // reset stored ff values
    for (int i = 0; i < numff; i++)
    {
        value1_ffs[i] = 0;
        value2_ffs[i] = ALLONES;
        id_unknown_ffs = 0;
    }
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

    while ((levelLen[smallest_level] == 0) && (smallest_level < maxlevels))
        smallest_level++;

    currLevel = smallest_level;
    if (currLevel < maxlevels)
    {
    	levelLen[currLevel]--;
        while ((levelLen[smallest_level] == 0) && (smallest_level < maxlevels))
            smallest_level++;
        //cout << "now smallest_level: " << smallest_level << endl;
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
    int gateN;
    int predecessor, successor; // the gateN of pred and succ
    int predecessor2; // to reference another predecessor if required
    int ff_number; // to distinguish ff1, ff2, ff3 ...
    bool evaluated = false; // a flag to record if output of this gate is evaluated
    int id_temp = 0; // temp id_unknown
    int *predList;
    int i, j; // for loops
    unsigned int val1, val2, tmpVal;
    int fanin_count; // the fanin count of current gate
    int activated_ff_count = 0;

    currLevel = 0;
    actLen = actFFLen = 0;
    while (currLevel < maxlevels || activated_ff_count)
    {
        gateN = retrieveEvent();
        fanin_count = fanin[gateN];
        if (gateN != -1)	// if a valid event
        {
            sched[gateN]= 0;
            evaluated = false;
            id_temp = 0;
            // variables for T_xor and T_xnor
            int xi_count = 0;              //num of xs in predecessors' output
            int xi_bar_count = 0;          //num of xs_bar in predecessors' output
            int one_of_the_id_unknown = 0; //of predecessors
            bool has_unknown = 0;
            bool has_real_unknown = 0;
            bool has_controlling = 0;
            switch (gtype[gateN])
            {
            case T_and:
                /************* X distinguished**************/
                // If two predecessors are unknowns
                // 1) have the same id and complement each other => they squash!
                // 2) have different id => generate a new unknown
                for (i = 0; i < fanin_count; i++)
                {
                    for (j = i+1; j < fanin_count; j++)
                    {
                        predecessor = inlist[gateN][i];
                        predecessor2 = inlist[gateN][j];
                        if ((id_unknown[predecessor] != 0) && (id_unknown[predecessor2] !=0))
                        {
                            // 1) have the same id and complement each other
                            if (id_unknown[predecessor] == id_unknown[predecessor2])
                            {
                                if (value1[predecessor] != value1[predecessor2])
                                {
                                    val1 = 0;
                                    val2 = 0;
                                    id_temp = 0;
                                    //cout << "Xs from PI squashed on gate " << gateN << "!\n";
                                    evaluated = true;
                                }
                            }
                            // 2) have different id
                            else
                            {
                                val1 = 0;
                                val2 = ALLONES;
                                id_unknown_max++;
                                id_temp = id_unknown_max;
                                evaluated = true;
                            }
                        }
                    }
                    if (value1[predecessor] != value2[predecessor])
                    {
                        has_unknown = 1;
                    }
                    if (id_unknown[predecessor] != 0)
                    {
                        has_real_unknown = 1;
                    }
                    if (value1[predecessor] && value1[predecessor])
                    {
                        has_controlling = 1;
                    }
                }
                if (evaluated)
                    break;
                else if (has_unknown && !has_real_unknown && !has_controlling)
                {
                    // new X with max id
                    val1 = 0;
                    val2 = ALLONES;
                    id_unknown_max++;
                    id_temp = id_unknown_max;
                }
                else
                // could be:
                // a) no unknown => output = 0 or 1
                // b) all unknowns have the same id, no complements => output = 0 or u
                {
                    val1 = val2 = ALLONES;
                    for (i = 0; i < fanin[gateN]; i++)
                    {
                        predecessor = inlist[gateN][i];
                        val1 &= value1[predecessor];
                        val2 &= value2[predecessor];
                        if (id_unknown[predecessor] != 0)
                        {
                            id_temp = id_unknown[predecessor];
                        }
                    }

                    //if ((val1 == 0) && (val2 == 0))
                    //    id_temp = 0;
                    if ((val1 == 0) && (val2 == 0)){
                        //debug
                        //if (id_temp)
                        //    cout << "X, id = " << id_temp << " is blocked at gate " << gateN << endl;
                        id_temp = 0;
                    }
                }
                break;
            case T_nand:
                /************* X distinguished**************/
                // If two predecessors are unknowns
                // 1) have the same id and complement each other => they squash!
                // 2) have different id => generate a new unknown
                for (i = 0; i < fanin_count; i++)
                {
                    for (j = i+1; j < fanin_count; j++)
                    {
                        predecessor = inlist[gateN][i];
                        predecessor2 = inlist[gateN][j];
                        if ((id_unknown[predecessor] != 0) && (id_unknown[predecessor2] !=0))
                        {
                            // 1) have the same id and complement each other
                            if (id_unknown[predecessor] == id_unknown[predecessor2])
                            {
                                //debug
                                //cout << "same id_unknown: " << id_unknown[predecessor] << " ";
                                //cout << id_unknown[predecessor2] << endl;
                                if (value1[predecessor] != value1[predecessor2])
                                {
                                    val1 = ALLONES;
                                    val2 = ALLONES;
                                    id_temp = 0;
                                    //cout << "Xs from PI squashed on gate " << gateN << "!\n";
                                    evaluated = true;
                                }
                            }
                            // 2) have different id
                            else
                            {
                                val1 = 0;
                                val2 = ALLONES;
                                id_unknown_max++;
                                id_temp = id_unknown_max;
                                //debug
                                //cout << "new unknown at gate " << gateN << ", id = " << id_temp << endl;
                                evaluated = true;
                            }
                        }
                    }
                    if (value1[predecessor] != value2[predecessor])
                    {
                        has_unknown = 1;
                    }
                    if (id_unknown[predecessor] != 0)
                    {
                        has_real_unknown = 1;
                    }
                    if (value1[predecessor] && value1[predecessor])
                    {
                        has_controlling = 1;
                    }
                }
                if (evaluated)
                    break;
                else if (has_unknown && !has_real_unknown && !has_controlling)
                {
                    // new X with max id
                    val1 = 0;
                    val2 = ALLONES;
                    id_unknown_max++;
                    id_temp = id_unknown_max;
                }
                else
                // a) no unknown => output = 0 or 1
                // b) all unknowns have the same id, no complements => output = 1 or u
                {
                    val1 = val2 = ALLONES;
                    for (i = 0; i < fanin[gateN]; i++)
                    {
                        predecessor = inlist[gateN][i];
                        val1 &= value1[predecessor];
                        val2 &= value2[predecessor];
                        if (id_unknown[predecessor] != 0)
                            id_temp = id_unknown[predecessor];
                    }
                    val1 = ALLONES ^ val1;
                    val2 = ALLONES ^ val2;
                    if ((val1 == ALLONES) && (val2 == ALLONES)){
                        //debug
                        /*
                        if (id_temp)
                            cout << "X, id = " << id_temp << " is blocked at gate " << gateN << endl;
                        */
                        id_temp = 0;
                    }
                }
                break;
            case T_or:
                /************* X distinguished**************/
                // If two predecessors are unknowns
                // 1) have the same id and complement each other => they squash!
                // 2) have different id => generate a new unknown
                for (i = 0; i < fanin_count; i++)
                {
                    for (j = i+1; j < fanin_count; j++)
                    {
                        predecessor = inlist[gateN][i];
                        predecessor2 = inlist[gateN][j];
                        if ((id_unknown[predecessor] != 0) && (id_unknown[predecessor2] !=0))
                        {
                            // 1) have the same id and complement each other
                            if (id_unknown[predecessor] == id_unknown[predecessor2])
                            {
                                if (value1[predecessor] != value1[predecessor2])
                                {
                                    val1 = 1;
                                    val2 = 1;
                                    id_temp = 0;
                                    //cout << "Xs from PI squashed on gate " << gateN << "!\n";
                                    evaluated = true;
                                }
                            }
                            // 2) have different id
                            else
                            {
                                val1 = 0;
                                val2 = ALLONES;
                                id_unknown_max++;
                                id_temp = id_unknown_max;
                                evaluated = true;
                            }
                        }
                    }
                    if (value1[predecessor] != value2[predecessor])
                    {
                        has_unknown = 1;
                    }
                    if (id_unknown[predecessor] != 0)
                    {
                        has_real_unknown = 1;
                    }
                    if (value1[predecessor] && value1[predecessor])
                    {
                        has_controlling = 1;
                    }
                }
                if (evaluated)
                    break;
                else if (has_unknown && !has_real_unknown && !has_controlling)
                {
                    // new X with max id
                    val1 = 0;
                    val2 = ALLONES;
                    id_unknown_max++;
                    id_temp = id_unknown_max;
                }
                else
                // a) no unknown => output = 0 or 1
                // b) all unknowns have the same id, no complements => output = 1 or u
                {
                    val1 = val2 = 0;
                    for (i = 0; i < fanin[gateN]; i++)
                    {
                        predecessor = inlist[gateN][i];
                        val1 |= value1[predecessor];
                        val2 |= value2[predecessor];
                        if (id_unknown[predecessor] != 0)
                            id_temp = id_unknown[predecessor];
                    }
                    if ((val1 == ALLONES) && (val2 == ALLONES))
                        id_temp = 0;
                }
                break;
            case T_nor:
                /************* X distinguished**************/
                // If two predecessors are unknowns
                // 1) have the same id and complement each other => they squash!
                // 2) have different id => generate a new unknown
                for (i = 0; i < fanin_count; i++)
                {
                    predecessor = inlist[gateN][i];
                    printGateValue(predecessor);
                    for (j = i+1; j < fanin_count; j++)
                    {
                        predecessor2 = inlist[gateN][j];
                        if ((id_unknown[predecessor] != 0) && (id_unknown[predecessor2] !=0))
                        {
                            // 1) have the same id and complement each other
                            if (id_unknown[predecessor] == id_unknown[predecessor2])
                            {
                                if (value1[predecessor] != value1[predecessor2])
                                {
                                    val1 = 0;
                                    val2 = 0;
                                    id_temp = 0;
                                    //cout << "Xs from PI squashed on gate " << gateN << "!\n";
                                    evaluated = true;
                                }
                            }
                            // 2) have different id
                            else
                            {
                                val1 = 0;
                                val2 = ALLONES;
                                id_unknown_max++;
                                id_temp = id_unknown_max;
                                evaluated = true;
                            }
                        }
                    }
                    if (value1[predecessor] != value2[predecessor])
                    {
                        has_unknown = 1;
                    }
                    if (id_unknown[predecessor] != 0)
                    {
                        has_real_unknown = 1;
                    }
                    if (value1[predecessor] && value1[predecessor])
                    {
                        has_controlling = 1;
                    }
                }
                if (evaluated)
                    break;
                else if (has_unknown && !has_real_unknown && !has_controlling)
                {
                    // new X with max id
                    val1 = 0;
                    val2 = ALLONES;
                    id_unknown_max++;
                    id_temp = id_unknown_max;
                }
                else
                // a) no unknown => output = 0 or 1
                // b) all unknowns have the same id, no complements => output = 0 or u
                {
                    val1 = val2 = 0;
                    for (i = 0; i < fanin[gateN]; i++)
                    {
                        predecessor = inlist[gateN][i];
                        val1 |= value1[predecessor];
                        val2 |= value2[predecessor];
                        if (id_unknown[predecessor] != 0)
                            id_temp = id_unknown[predecessor];
                    }
                }
                val1 = ALLONES ^ val1;
                val2 = ALLONES ^ val2;
                if ((val1 == 0) && (val2 == 0))
                        id_temp = 0;
                break;
            case T_not:
                predecessor = inlist[gateN][0];
                val1 = ALLONES ^ value1[predecessor];
                val2 = ALLONES ^ value2[predecessor];
                id_temp = id_unknown[predecessor];
                // debug
                /*
                cout << "not " << gateN << " has input gate " << predecessor << " ";
                cout << value1[predecessor] << " & " << value2[predecessor];
                cout << " id = " << id_unknown[predecessor] << endl;
                */
                break;
            case T_buf:
                predecessor = inlist[gateN][0];
                val1 = value1[predecessor];
                val2 = value2[predecessor];
                id_temp = id_unknown[predecessor];
                break;
            case T_dff:
                predecessor = inlist[gateN][0];
                ff_number = ffMap[gateN]; //ff1, ff2, or ff3? ...
                val1 = value1[gateN];
                val2 = value2[gateN];
                id_temp = id_unknown[gateN];
                value1_ffs[ff_number] = value1[predecessor];
                value2_ffs[ff_number] = value2[predecessor];
                id_unknown_ffs[ff_number] = id_unknown[predecessor];
                activated_ff_count--;

                //debug
                // output values
                /*
                if (id_temp == 0)
                {
                    if (val1 && val2)
                    {
                        cout << "ff[" << ff_number << "] output 1" << endl;
                    }
                    else if (!val1 && !val2)
                    {
                        cout << "ff[" << ff_number << "] output 0" << endl;
                    }
                    else
                    {
                        cout << "ff[" << ff_number << "] output initial unknown" << endl;
                    }
                }
                else
                {
                    cout << "ff[" << ff_number << "] output unknown, id " << id_temp << endl;
                }
                // stored values
                if (id_unknown_ffs[ff_number] == 0)
                {
                    if (value1_ffs[ff_number] && value2_ffs[ff_number])
                    {
                        cout << "ff[" << ff_number << "] stored 1" << endl;
                    }
                    else if (!value1_ffs[ff_number] && !value2_ffs[ff_number])
                    {
                        cout << "ff[" << ff_number << "] stored 0" << endl;
                    }
                    else
                    {
                        cout << "ff[" << ff_number << "] stored initial unknown" << endl;
                    }
                }
                else
                {
                    cout << "ff[" << ff_number << "] stored unknown, id " << id_unknown_ffs[ff_number] << endl;
                }
                */

                // This part was done by Dr. Hsiao
                // Functionally unknown
                actFFList[actFFLen] = gateN;
                actFFLen++;
                break;
            case T_xor:
                predList = inlist[gateN];
                // if two predecessors has different id => generate a new unknown
                for (i = 0; i < fanin_count; i++)
                {
                    for (j = i+1; j < fanin_count; j++)
                    {
                        predecessor = inlist[gateN][i];
                        predecessor2 = inlist[gateN][j];
                        if ((id_unknown[predecessor] != 0) && (id_unknown[predecessor2] !=0))
                        {
                            if (id_unknown[predecessor] != id_unknown[predecessor2]) {
                                // generate a new unknown
                                val1 = 0;
                                val2 = ALLONES;
                                id_unknown_max++;
                                id_temp = id_unknown_max;
                                break;
                            }
                        }
                    }
                }

                // count Xi and Xi_bar
                for (i = 0; i < fanin_count; i++)
                {
                    predecessor = inlist[gateN][i];
                    if (id_unknown[predecessor] != 0)
                    {
                        if (value1[predecessor] == 0)
                        {
                            xi_count++;
                        }
                        else
                        {
                            xi_bar_count++;
                        }
                        one_of_the_id_unknown = id_unknown[predecessor];
                    }
                }

                if ((xi_count + xi_bar_count) == 0)
                {
                    // no x => evaluate normally
                    for(i = 1; i < fanin[gateN]; i++)
                    {
                        predecessor = predList[i];
                        val1 = val1 ^ value1[predecessor];
                        val2 = val2 ^ value2[predecessor];
                        id_temp = 0;

                        // Dr. Shiao's evaluation equation
                        /*
                        tmpVal = ALLONES^(((ALLONES^value1[predecessor]) &
                                           (ALLONES^val1)) | (value2[predecessor]&val2));
                        val2 = ((ALLONES^value1[predecessor]) & val2) |
                               (value2[predecessor] & (ALLONES^val1));
                        val1 = tmpVal;
                        */
                    }
                    break;
                }
                else if (!(xi_count%2) && !(xi_bar_count%2))
                {
                    // xi and xi_bar are both even => squashed, 0
                    //cout << "Xs from PI squashed on gate " << gateN << "!\n";
                    val1 = 0;
                    val2 = 0;
                    id_temp = 0;
                    break;
                }
                else if ((xi_count%2) && (xi_bar_count%2))
                {
                    // xi and xi_bar are both odd => squashed, 1
                    //cout << "Xs from PI squashed on gate " << gateN << "!\n";
                    val1 = ALLONES;
                    val2 = ALLONES;
                    id_temp = 0;
                    break;
                }
                else if (!(xi_count%2) && (xi_bar_count%2))
                {
                    // xi even, xi_bar odd => x_bar, same id
                    val1 = 0;
                    val2 = ALLONES;
                    id_temp = one_of_the_id_unknown;
                    break;
                }
                else if ((xi_count%2) && !(xi_bar_count%2))
                {
                    // xi odd, xi_bar even => x, same id
                    val1 = 0;
                    val2 = ALLONES;
                    id_temp = one_of_the_id_unknown;
                    break;
                }
                break;

            case T_xnor:
                predList = inlist[gateN];
                // if two predecessors has different id => generate a new unknown
                for (i = 0; i < fanin_count; i++)
                {
                    for (j = i+1; j < fanin_count; j++)
                    {
                        predecessor = inlist[gateN][i];
                        predecessor2 = inlist[gateN][j];
                        if ((id_unknown[predecessor] != 0) && (id_unknown[predecessor2] !=0))
                        {
                            if (id_unknown[predecessor] != id_unknown[predecessor2]) {
                                // generate a new unknown
                                val1 = 0;
                                val2 = ALLONES;
                                id_unknown_max++;
                                id_temp = id_unknown_max;
                                break;
                            }
                        }
                    }
                }

                // count Xi and Xi_bar
                for (i = 0; i < fanin_count; i++)
                {
                    predecessor = inlist[gateN][i];
                    if (id_unknown[predecessor] != 0)
                    {
                        if (value1[predecessor] == 0)
                        {
                            xi_count++;
                        }
                        else
                        {
                            xi_bar_count++;
                        }
                        one_of_the_id_unknown = id_unknown[predecessor];
                    }
                }

                if ((xi_count + xi_bar_count) == 0)
                {
                    // no x => evaluate normally
                    for(i = 1; i < fanin[gateN]; i++)
                    {
                        predecessor = predList[i];
                        val1 = val1 ^ value1[predecessor];
                        val2 = val2 ^ value2[predecessor];
                        id_temp = 0;

                    }
                    // not at the end
                    val1 ^= ALLONES;
                    val2 ^= ALLONES;
                    break;
                }
                else if (!(xi_count%2) && !(xi_bar_count%2))
                {
                    // xi and xi_bar are both even => squashed, 1
                    cout << "Xs from PI squashed on gate " << gateN << "!\n";
                    val1 = ALLONES;
                    val2 = ALLONES;
                    id_temp = 0;
                    break;
                }
                else if ((xi_count%2) && (xi_bar_count%2))
                {
                    // xi and xi_bar are both odd => squashed, 0
                    cout << "Xs from PI squashed on gate " << gateN << "!\n";
                    val1 = 0;
                    val2 = 0;
                    id_temp = 0;
                    break;
                }
                else if (!(xi_count%2) && (xi_bar_count%2))
                {
                    // xi even, xi_bar odd => x, same id
                    val1 = ALLONES;
                    val2 = 0;
                    id_temp = one_of_the_id_unknown;
                    break;
                }
                else if ((xi_count%2) && !(xi_bar_count%2))
                {
                    // xi odd, xi_bar even => x_bar, same id
                    val1 = ALLONES;
                    val2 = 0;
                    id_temp = one_of_the_id_unknown;
                    break;
                }
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
            } // switch

            //debug
            /*
            char char_val1 = val1 ? '1':'0';
            char char_val2 = val2 ? '1':'0';
            cout << "gate" << gateN << " has output: " << char_val1 << " & " << char_val2;
            cout << " id = " << id_temp << endl;
            */

            // if gate value changed
            if ((val1 != value1[gateN]) || (val2 != value2[gateN]) || (id_temp != id_unknown[gateN]))
            {
                value1[gateN] = val1;
                value2[gateN] = val2;
                id_unknown[gateN] = id_temp;

                for (i=0; i<fanout[gateN]; i++)
                {
                    successor = fnlist[gateN][i];
                    sucLevel = levelNum[successor];
                    if (sched[successor] == 0)
                    {
                        if (gtype[successor] == T_dff)
                        {
                            activated_ff_count++;
                            insertEvent(sucLevel, successor);
                        }
                        else if (sucLevel != 0)
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


    /******** old FF handle function by Dr. Hsiao********/
    /************** temporarily blocked******************/
    // now re-insert the activation list for the FF's
    /*
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
    */
}

////////////////////////////////////////////////////////////////////////
// LogicSim(char** input_vectors)
// Do the sequential logic simulation with input vectors
////////////////////////////////////////////////////////////////////////
void gateLevelCkt::LogicSim(char** input_vectors, int timeframes)
{
    //cout << "sizeof(input_vectors) = " << sizeof(input_vectors) << endl;
    //cout << "sizeof(input_vectors[0]) = " << sizeof(input_vectors[0]) << endl;
    //int timeframes = sizeof(input_vectors);
    //cout << "timeframes =" << timeframes << endl;

    for (int i = 0; i < timeframes; i++)
    {
        //cout << "vector for timeframe #" << i << ": " << input_vectors[i] << endl;
        applyFF();
        //cout << "applyFF succeeded" << endl;
        applyVector(input_vectors[i]);
        //cout << "applyVector succeeded" << endl;
        goodsim();
        //cout << "goodsim succeeded" << endl;
        observeFFs();
        observeOutputs();


    }
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

void gateLevelCkt::observeFFs()
{
    int i;
    cout << "FFs:\t";
    for (i=0; i<numff; i++)
    {
        if (value1_ffs[i] && value2_ffs[i])
            cout << "1";
        else if (!value1_ffs[i] && !value2_ffs[i])
            cout << "0";
        else
            cout << "X";
    }
    cout << endl;

}

void gateLevelCkt::printGateValue(int gateN)
{
    unsigned int val1 = value1[gateN];
    unsigned int val2 = value2[gateN];
    int id = id_unknown[gateN];
/*
    if (val1 && val2)
    {
        cout << "gate[" << gateN << "] output 1, id = " << id << endl;
    }
    else if (!val1 && !val2)
    {
        cout << "gate[" << gateN << "] output 0, id = " << id << endl;
    }
    else
    {
        cout << "gate[" << gateN << "] output unknown, id = " << id << endl;
    }
*/
}

unsigned int* gateLevelCkt::getValue1FFs()
{
    return value1_ffs;
}

unsigned int* gateLevelCkt::getValue2FFs()
{
    return value2_ffs;
}

int* gateLevelCkt::getIDUnknownFFs(){
    return id_unknown_ffs;
}
