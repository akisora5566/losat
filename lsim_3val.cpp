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
#include <sys/times.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
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
    int numpri;		// number of PIs
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
public:
    int numff;		// number of FF's
    unsigned int *RESET_FF1;	// value of reset ffs read from *.initState
    unsigned int *RESET_FF2;	// value of reset ffs read from *.initState

    gateLevelCkt(string);	// constructor
    void setFaninoutMatrix();	// builds the fanin-out map matrix

    void applyVector(char *);	// apply input vector
    void resetIDunknown(); // reset the value of gates to initial state

    // simulator information
    void setupWheel(int, int);
    void insertEvent(int, int);
    int retrieveEvent();
    void goodsim();		// logic sim (no faults inserted)

    void setTieEvents();	// inject events from tied nodes

    void observeOutputs();	// print the fault-free outputs
    void printGoodSig(ofstream, int);	// print the fault-free outputs to *.sig

    char *goodState;		// good state (without scan)
};

//////////////////////////////////////////////////////////////////////////
// Global Variables

char vector[5120];
gateLevelCkt *circuit;
int OBSERVE, INIT0;
int vecNum=0;
int numTieNodes;
int TIES[512];


//////////////////////////////////////////////////////////////////////////
// Functions start here


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

    vector[0] = toupper(thisChar); //in case input vec would contain 'x'
    if (vector[0] == 'E')
	return (0);

    for (i=1; i<vecSize; i++)
    {
        inFile >> thisChar;
        vector[i] = toupper(thisChar); //in case input vec would contain 'x'
    }
    vector[i] = EOS;

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
            //buggy statement here!!
            circuit->resetIDunknown();  // reset value before simulation
            cout << "vector #" << vecNum << ": " << vector << "\n";
            circuit->applyVector(vector);
            circuit->goodsim();      // simulate the vector
            if (OBSERVE) circuit->observeOutputs();
                vecNum++;
        } // if (moreVec == 1)
    } // while (getVector...)

    return (vecNum);
}


// main()
int main(int argc, char *argv[])
{
    ifstream vecFile;
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

    circuit = new gateLevelCkt(cktName);

    vecFile.open(vecName.c_str(), ios::in);
    if (!vecFile)
    {
        cerr << "Can't open " << vecName << "\n";
        exit(-1);
    }

    vecFile >> vecWidth;

    circuit->setTieEvents();
    totalNumVec = logicSimFromFile(vecFile, vecWidth);
    vecFile.close();

    // output results
    end = clock();
    ut = (double) (end-start);
    cout << "Number of vectors: " << totalNumVec << "\n";
    cout << "Number of clock cycles elapsed: " << ut << "\n";
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
    id_unknown_max = 0;

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
    int max_pi_id_unknown = 1;
    
    for (i = 0; i < numpri; i++)
    {
        origVal1 = value1[inputs[i]] & 1;
        origVal2 = value2[inputs[i]] & 1;
        if (id_unknown[inputs[i]] == 0)
        {
            if ((origVal1 == 1) && (origVal2 == 1))
                origBit = '1';
            else if ((origVal1 == 0) && (origVal2 == 0))
                origBit = '0';
        }
        else
            origBit = 'X';
        
        // if new input is X, we must reapply
        // else if new input is 0 or 1, we reapply only if it changed   
        if ((vec[i] != 'X') || (origBit != vec[i]))
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
                    id_unknown[inputs[i]] = max_pi_id_unknown;
                    max_pi_id_unknown++;
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
// resetIDunknown()
//	This function reset the value of gates to initial state
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::resetIDunknown()
{
    for (int i = 0; i < numgates; i++)
    {
        id_unknown[i] = 0;
        id_unknown_max = 0;
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
    int gateN; 
    int predecessor, successor; // the gateN of pred and succ
    int predecessor2; // to reference another predecessor if required
    bool evaluated = false; // a flag to record if output of this gate is evaluated
    int id_temp = 0; // temp id_unknown
    int *predList;
    int i, j; // for loops
    unsigned int val1, val2, tmpVal;
    int fanin_count; // the fanin count of current gate

    currLevel = 0;
    actLen = actFFLen = 0;
    while (currLevel < maxlevels)
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
                                    cout << "Xs from PI squashed on gate " << gateN << "!\n";
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
                }
                if (evaluated)
                    break;
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
                                    cout << "Xs from PI squashed on gate " << gateN << "!\n";
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
                }
                if (evaluated)
                    break;
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
                                    cout << "Xs from PI squashed on gate " << gateN << "!\n";
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
                }
                if (evaluated)
                    break;
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
                                    cout << "Xs from PI squashed on gate " << gateN << "!\n";
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
                }
                if (evaluated)
                    break;
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
                val1 = value1[predecessor];
                val2 = value2[predecessor];
                id_temp = id_unknown[predecessor];
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
                    cout << "Xs from PI squashed on gate " << gateN << "!\n";
                    val1 = 0;
                    val2 = 0;
                    id_temp = 0;
                    break;
                }
                else if ((xi_count%2) && (xi_bar_count%2))
                {
                    // xi and xi_bar are both odd => squashed, 1
                    cout << "Xs from PI squashed on gate " << gateN << "!\n";
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
            cout << "gate" << gateN << " has output: " << val1 << " & " << val2;
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

void gateLevelCkt::observeOutputs()
{
    int i;

    cout << "\t";
    for (i=0; i<numout; i++)
    {
	if (value1[outputs[i]] && value2[outputs[i]])
	    cout << "1";
	else if ((value1[outputs[i]] == 0) && (value2[outputs[i]] == 0))
	    cout << "0";
	else
	    cout << "X";
    }

    cout << "\n";
    for (i=0; i<numff; i++)
    {
	if (value1[ff_list[i]] && value2[ff_list[i]])
	    cout << "1";
	else if ((value1[ff_list[i]] == 0) && (value2[ff_list[i]] == 0))
	    cout << "0";
	else
	    cout << "X";
    }

    cout << "\n";
}

