///////////////////////////////////////////////////////////////////////
// Fault Simulator, written by Michael Hsiao
//   Began in 1992, revisions and additions of functions till 2013
///////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <strings.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <iostream>
//#include "header.h"

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
#define MAXFanout 101920
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
    int numFaultsInserted;	// number of faults inserted in circuit
    int fgates[32];	// list of faulty gate indices in gates[]
    int fIndices[32];	// list of the fault indices in circuit

    // fault information
    int *fGate;		// actual faults
    int *fIO;		// actual faults' input on gate
    int *fStuck;	// actual fault s-a-value
    int *fprelist;	// predecessor fault list
    int *fsuclist;	// successor fault list
    int *fstatus;	// status of faults
    int *numExcited;	// number of times a fault is excited
    int *no_faultdrop_fstatus;	// count of number of times fault detected
int *tmpDropTable;	// tmporary table for dropped faults (NO_FAULTDROP op)
int *tmpDropDet;	// tmporary dropped types (NO_FAULTDROP option)
    int *stfList;	// state faults list per fault
    int *stfCount;	// length of state faults list
    int **predOfSuccInput;	// predecessor of successor input-pin list
    int **succOfPredOutput;	// successor of predecessor output-pin list
    int faultId;	// fault id for this set of faults

    // fault state info
    int *fStateFF;
    char *fStateVal;
    int *fStateNext;
    int fStateSize;	// size of fState array
    int fStateFree;	// first free location in array

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
    char *inFFsched;		// faulty ff scheduled ?
    int *currId;		// fault id for each gate
    unsigned int *value1;
    unsigned int *value2;	// value of gate
    unsigned int *fvalue1;
    unsigned int *fvalue2;	// faulty value of gate
    int faultActive(int,int,int);	// this fault active?

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
    int totalNumFaults;	// number of faults in ckt
    int numFaults;	// total number of undetected faults
    int startfLoc; 	// start indix for undetected faults
    int numff;		// number of FF's
    unsigned int *RESET_FF1;	// value of reset ffs read from *.initState
    unsigned int *RESET_FF2;	// value of reset ffs read from *.initState

    gateLevelCkt(char *);	// constructor
    int insertFaults(int, int);	// insert faulty gate(s)
    void removeFaults();	// remove faulty gate(s)
    void setupFaults(char *);	// set up fault data
    void dropFault(int);	// drop a fault from fault list
    int tmpDropIndex;
    void tmpDrop(int);		// tmporary drop a fault (NO_FAULTDROP option)
    void restoreDrop();		// restore the tmp dropped faults (NO_FAULTDROP)
    void setFaninoutMatrix();	// builds the fanin-out map matrix
    void setupfState();		// set up fault state array
    int fStateNewLoc();		// get a new location from fState array
    void fStateFreeLoc(int);	// free a location in the fState array

    void applyVector(char *);	// apply input vector
    void clearFF(int);		// reset all the FF's fvalues to good values
    void saveStateFaults();	// saves all the faults propgated to FF's
    void setStateFaults();	// sets the faults at ff's for current frame

    // simulator information
    void setupWheel(int, int);
    void insertEvent(int, int);
    int retrieveEvent();
    void goodsim();		// logic sim (no faults inserted)
    void globalFsim();		// fault sim (for all faults)

    void setTieEvents();	// inject events from tied nodes

    void observeOutputs();	// print the fault-free outputs
    void printGoodSig(FILE *, int);	// print the fault-free outputs to *.sig
    void printSig(FILE *);	// print the faulty outputs to *.sig
    void writeOutput(FILE *, FILE*, FILE*);	// write the det and ufl files
    void writeOutputNDET(FILE*);	// write the ndet file

    char *goodState;		// good state (without scan)
};

//////////////////////////////////////////////////////////////////////////
// Global Variables

char vector[5120];
gateLevelCkt *circuit;
int OUTPUT_DET, NO_FAULTDROP, OBSERVE, INIT0;
int SCANID;
int vecNum=0;
int numTieNodes;
int TIES[512];


//////////////////////////////////////////////////////////////////////////
// Functions start here


//////////////////////////////////////////////////////////////////////////
// returns 1 if a regular vector, returns 2 if a scan
//////////////////////////////////////////////////////////////////////////
int getVector(FILE *inFile, int vecSize)
{
    int i;
    char thisChar;

    fscanf(inFile, "%c", &thisChar);
    while ((thisChar == SPACE) || (thisChar == RETURN))
	fscanf(inFile, "%c", &thisChar);

    vector[0] = thisChar;
    if (vector[0] == 'E')
	return (0);

    for (i=1; i<vecSize; i++)
	fscanf(inFile, "%c", &vector[i]);
    vector[i] = EOS;

    fscanf(inFile, "%c", &thisChar);

    // read till end of line
    while ((thisChar != RETURN) && (thisChar != EOF))
    	fscanf(inFile, "%c", &thisChar);

    return(1);
}

void faultSim(FILE *sigFile)
{
    int numfaults, fullScaleIndex;

      // applies vector to the PIs and sets up initial events to wheel
      circuit->applyVector(vector);
      // reset the fault index
      fullScaleIndex = circuit->startfLoc;

      // for all faults in the circuit
      circuit->goodsim();
if (NO_FAULTDROP) circuit->printGoodSig(sigFile, vecNum);
      circuit->clearFF(0);	// first reset the FF's to the good ckt values
      while (fullScaleIndex != -1)
      {
    	fullScaleIndex = circuit->insertFaults(31, fullScaleIndex);

	circuit->setStateFaults();	// set state faults for current faults
    	circuit->globalFsim();
	numfaults = circuit->numFaults;

if (NO_FAULTDROP) circuit->printSig(sigFile);

    	circuit->removeFaults();	// stuck at faults
	circuit->clearFF(1);	// reset the FF's to the good ckt values
      }	// while (fullScaleIndex...)

}

// faultSimFromFile() - fault simulates all vectors in vecFile on the circuit
int faultSimFromFile(FILE *vecFile, int vecWidth, FILE *sigFile)
{
    int moreVec;

    moreVec = 1;
    while (moreVec && (circuit->numFaults > 0))
    {
     moreVec = getVector(vecFile, vecWidth);
     if (moreVec == 1)
     {
	faultSim(sigFile);	// simulate the vector

if (NO_FAULTDROP)
{
printf("vector #%d detected %d faults.\n", vecNum, circuit->tmpDropIndex);
circuit->restoreDrop();
}
else
printf("vector #%d detected %d faults.\n", vecNum, circuit->totalNumFaults - circuit->numFaults);
if (OBSERVE) circuit->observeOutputs();

      vecNum++;
     }	// if (moreVec == 1)
    }	// while (getVector...)

    return (vecNum);
}

// this function writes the faults to output files
void gateLevelCkt::writeOutputNDET(FILE *detFile)
{
    int i;

    for (i=0; i<totalNumFaults; i++)
    {
	if (no_faultdrop_fstatus[i] > 0)	// detected
	    fprintf(detFile, "%d %d %d; detected %d times\n", fGate[i], fIO[i], fStuck[i], no_faultdrop_fstatus[i]);
	else
	{
	    if (fstatus[i] == POTENTIAL)	// potentially detected
	    {
	    	fprintf(detFile, "%d %d %d;\t*** potentially detected ***\n",
			fGate[i], fIO[i], fStuck[i]);
	    }
	    else	// undetected
	    	fprintf(detFile, "%d %d %d; not detected\n", fGate[i], fIO[i], fStuck[i]);
	}
    }
}

// this function writes the faults to output files
void gateLevelCkt::writeOutput(FILE *detFile, FILE *uflFile, FILE *eFile)
{
    int i;

    for (i=0; i<totalNumFaults; i++)
    {
	if (fstatus[i] == LOW_DETECT)	// detected
	    fprintf(detFile, "%d %d %d;\n", fGate[i], fIO[i], fStuck[i]);
	else
	{
	    if (fstatus[i] == POTENTIAL)	// potentially detected
	    {
	    	fprintf(uflFile, "%d %d %d;\t*** potentially detected ***\n",
			fGate[i], fIO[i], fStuck[i]);
	    }
	    else	// undetected
	    	fprintf(uflFile, "%d %d %d;\n", fGate[i], fIO[i], fStuck[i]);

	    if (fstatus[i] == EXCITED_1_LEVEL)	// excited but not detected
		fprintf(eFile, "%d %d %d; *** excited %d times\n", fGate[i], fIO[i], fStuck[i], numExcited[i]);
	}
    }
}

// main()
main(int argc, char *argv[])
{
    FILE *vecFile, *detFile, *ndetFile, *uflFile, *frsFile, *exFile, *sigFile;
    char cktName[256], vecName[256], nDetName[256];
    char uflName[256], detName[256], exName[256];
    char frsName[256], sigName[256];
    int totalNumVec, vecWidth, i;
    int nameIndex;
    float ut, st;
    struct tms t_buf;

    if ((argc != 2) && (argc != 3))
    {
	fprintf(stderr, "Usage: %s [-dino] <ckt>\n", argv[0]);
	fprintf(stderr, "The -d option is for dumping detected faults on the screen.\n");
	fprintf(stderr, "The -n option is for NO FAULT DROPPING. (.ndet & .sig file generated)\n");
	fprintf(stderr, "The -i option is to begin the circuit in a state in *.initState.\n");
	fprintf(stderr, "and -o option is to OBSERVE fault-free outputs and FF's.\n");
	fprintf(stderr, "Vectors (*.vec) may have SCAN FFs as well.\n");
	fprintf(stderr, "\tResults are written to ckt.det and ckt.ufl files\n");
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
		case 'd':
	    	    OUTPUT_DET = 1;
		    break;
		case 'n':
	    	    NO_FAULTDROP = 1;
		    break;
		case 'o':
		    OBSERVE = 1;
		    break;
		case 'i':
		    INIT0 = 1;
		    break;
		default:
	    	    fprintf(stderr, "Invalid option: %s\n", argv[1]);
	    	    fprintf(stderr, "Usage: %s [-dnio] <ckt> <type>\n", argv[0]);
		    fprintf(stderr, "The -d option is for dumping detected faults on the screen.\n");
	    	    fprintf(stderr, "The -n option is for NO FAULT DROPPING.\n");
	    	    fprintf(stderr, "The -i option is to begin the circuit in a state in *.initState.\n");
	    	    fprintf(stderr, "and o option is to OBSERVE fault-free outputs.\n");
		    exit(-1);
		    break;
	    } 	// switch
	    i++;
	}	// while
    }
    else	// no option
    {
	nameIndex = 1;
	OUTPUT_DET = 0;
	NO_FAULTDROP = 0;
	OBSERVE = 0;
	INIT0 = 0;
    }

    strcpy(cktName, argv[nameIndex]);
    strcpy(vecName, argv[nameIndex]);
    strcpy(uflName, argv[nameIndex]);
    strcpy(detName, argv[nameIndex]);
    strcpy(frsName, argv[nameIndex]);
    strcpy(nDetName, argv[nameIndex]);
    strcpy(exName, argv[nameIndex]);
    strcat(vecName, ".vec");
    strcat(uflName, ".ufl");
    strcat(detName, ".det");
    strcat(frsName, ".frs");
    strcat(nDetName, ".ndet");
    strcat(exName, ".excited_but_undetected");
	
	std::string line;
    int temp_i;
    std::ifstream myfile ("output.vec");
    //std::ofstream vec_file (vecName);
    std::vector < std::vector <std::string> > vectors;
    std::vector< std::string> temp; 
  //vec_file.open();
  
    //////////////////Create a array of all vectors//////////////////////////////////////
    if (myfile.is_open())
    {
    	while(getline (myfile,line) )
      	{
    		temp_i =0;
      		if(line == "END"){
      			temp_i =1;
	  		}
	  		
	  		//std::cout << line << "\n";
	  		
	  		if(!temp_i)temp.push_back(line);
	  		else{
	  			vectors.push_back(temp);
	  			/*
				  std::cout << "\n The temp is :";
	  			for(int re =0; re < temp.size(); re++){
	  				std::cout << temp[re];	
				}
				std::cout <<"\n";
	  			*/
				//std::cin  >> temp_i;
	  			temp.clear();
	  		}
    	}
    	myfile.close();
    	//vec_file.close();
  	}
	
	//-----------------------------------PRINT VECTORS-------------------------
	/*
	for(int o = 0; o < vectors.size(); o++){
		for(int u = 0; u < vectors[o].size(); u++){
			std::cout<< vectors[o][u];
		}
		std::cout << "\n";
	}
	*/
	//std::cout << "PRINTED ALL VECTORS...............................\n";
	//std::cin >> temp_i;
	/*
	std::cout << "All vectors are ready!!\n";
	std::cout << vectors.size() << " vectors are generated. \nPress a number to continue..............\n";
	//std::cin >> temp_i;
	std::cout << "\n";
	*/
	SCANID = 1;
	
	circuit = new gateLevelCkt(cktName);
	for(int i =0; i < vectors.size(); i++){			//fault simulate all vectors
		//circuit = new gateLevelCkt(cktName);
		
		std::ofstream vec_file (vecName);
		for(int j=0; j < vectors[i].size(); j++){
			vec_file << vectors[i][j];
			vec_file << "\n";
		}
		vec_file << "END\n";
		vec_file.close();

		vecFile = fopen(vecName, "r");
    	if (vecFile == NULL)
    	{
			fprintf(stderr, "Can't open %s\n", vecName);
			exit(-1);
    	}

    	fscanf(vecFile, "%d", &vecWidth);

    	if (NO_FAULTDROP)
    	{
      		strcpy(sigName, argv[nameIndex]);
      		strcat(sigName, ".sig");
			sigFile = fopen(sigName, "w");
			if (sigFile == NULL)
    		{
	  			fprintf(stderr, "Can't open %s\n", sigName);
	  			exit(-1);
    		}
    	}
    	circuit->setTieEvents();
    	totalNumVec = faultSimFromFile(vecFile, vecWidth, sigFile);
    	fclose(vecFile);

    	// output results
    	frsFile = fopen(frsName, "w");
    	if (frsFile == NULL)
    	{
			fprintf(stderr, "Can't open %s\n", frsName);
			exit(-1);
    	}
    	times(&t_buf);         // get the time
    	ut = (float) t_buf.tms_utime / HZ; // get the user time
    	st = (float)t_buf.tms_stime / HZ; // get the system time
    	fprintf(frsFile, "Number of vectors: %d\n", totalNumVec);
    	fprintf(frsFile, "Detected %d in total %d faults, FC = %f\n",
		circuit->totalNumFaults-circuit->numFaults, circuit->totalNumFaults,
		(float) (circuit->totalNumFaults-circuit->numFaults) /
		circuit->totalNumFaults); 
    	fprintf(frsFile, "Time: user %f, system %f\n", ut, st);

    	detFile = fopen(detName, "w");
    	if (detFile == NULL)
    	{
			fprintf(stderr, "Can't open %s\n", detName);
			exit(-1);
    	}
    	uflFile = fopen(uflName, "w");
    	if (uflFile == NULL)
    	{
			fprintf(stderr, "Can't open %s\n", uflName);
			exit(-1);
    	}
    	
		exFile = fopen(exName, "w");
    	if (exFile == NULL)
    	{
			fprintf(stderr, "Can't open %s\n", exName);
			exit(-1);
    	}
    	
		circuit->writeOutput(detFile, uflFile, exFile);

    	if (NO_FAULTDROP)
    	{
      		fclose(sigFile);
      		ndetFile = fopen(nDetName, "w");
      		if (ndetFile == NULL)
      		{
				fprintf(stderr, "Can't open %s\n", nDetName);
				exit(-1);
      		}
      		circuit->writeOutputNDET(ndetFile);
    	}

    	fclose(frsFile);
    	fclose(detFile);
    	fclose(uflFile);
    	fclose(exFile);
		//delete circuit;
		circuit->clearFF(1);
		//std::cout << "Press a number to continue.................\n";
		//std::cin >> temp_i;		
	}
}


////////////////////////////////////////////////////////////////////////
inline void gateLevelCkt::insertEvent(int levelN, int gateN)
{
    levelEvents[levelN][levelLen[levelN]] = gateN;
    levelLen[levelN]++;
}

////////////////////////////////////////////////////////////////////////
// fStateFreeLoc()
//	This function frees up a location in fState Array.
////////////////////////////////////////////////////////////////////////

inline void gateLevelCkt::fStateFreeLoc(int loc)
{
    int oldFree = fStateFree;

    fStateFree = loc;
    fStateNext[fStateFree] = oldFree;
}

////////////////////////////////////////////////////////////////////////
// setupfState()
//	This function sets up the fState array.
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::setupfState()
{
    int i, j;

    // fault state size = 12 * cube_root(totalNumFaults) * cube_root(numff)
    for (i=1; i*i*i < totalNumFaults; i *= 2)
	;
    for (j=1; j*j*j < numff; j *= 2)
	;
    fStateSize = 12 * i * j;

    fStateVal = new char[fStateSize];
    fStateFF = new int[fStateSize];
    fStateNext = new int[fStateSize];
    fStateFree = 0;
    for (i=0; i<fStateSize; i++)
    {
	fStateFF[i] = -1;
 	fStateVal[i] = 2;	// don't care
	fStateNext[i] = i+1;
    }
    fStateNext[fStateSize-1] = -1;
}

////////////////////////////////////////////////////////////////////////
// fStateNewLoc()
//	This function returns a new address in the fState Array.
////////////////////////////////////////////////////////////////////////

int gateLevelCkt::fStateNewLoc()
{
    int i, loc = fStateFree;
    int oldSize = fStateSize;

    fStateFree = fStateNext[loc];
    if (fStateFree != -1)
    {
	return(loc);
    }
    else	// allocate more space
    {
fprintf(stderr, "DOUBLING fState Array %d\n", fStateSize);
fflush(stderr);
	fStateSize *= 2;
	fStateFF = (int *) realloc(fStateFF, fStateSize*sizeof(int));
        fStateVal = (char *) realloc(fStateVal, fStateSize*sizeof(char));
        fStateNext = (int *) realloc(fStateNext, fStateSize*sizeof(int));

        fStateFree = oldSize;
        for (i=oldSize; i<fStateSize; i++)
        {
	    fStateFF[i] = -1;
 	    fStateVal[i] = 2;	// don't care
            fStateNext[i] = i+1;
        }

	fStateNext[fStateSize-1] = -1;

	return(loc);
    }	// else allocate more space
}

////////////////////////////////////////////////////////////////////////
// saveStateFaults()
//	This function saves all faults propagated to the ff's by the
// current fault group.
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::saveStateFaults()
{
    int i, j;
    int mask, fIndex, ffIndex;
    int predecessor;
    unsigned int goodVal1, goodVal2, gateVal1, gateVal2, faultVal1, faultVal2;

    for (i=0; i < actLen; i++)
    {
	sched[activation[i]] = 0;
	predecessor = inlist[activation[i]][0];	// take input of the FF
	if (currId[predecessor] == faultId)
	{
	  mask = 1;
	  gateVal1 = fvalue1[predecessor];
	  gateVal2 = fvalue2[predecessor];
	  goodVal1 = value1[predecessor] & mask;
	  goodVal2 = value2[predecessor] & mask;

	  for (j=0; j < numFaultsInserted; j++)
	  {
	    mask = mask << 1;
	    goodVal1 = goodVal1 << 1;
	    goodVal2 = goodVal2 << 1;
	    faultVal1 = gateVal1 & mask;
	    faultVal2 = gateVal2 & mask;
	    fIndex = fIndices[j];
	    if (fstatus[fIndex] < LOW_DETECT)	// if fault isn't detected
	    {
	      if ((goodVal1 != faultVal1) || (goodVal2 != faultVal2))
	      {	// different from good ckt values
		  ffIndex = fStateNewLoc();
		  fStateFF[ffIndex] = activation[i];
		  if ((!faultVal1) && (faultVal2))	// don't care
		    fStateVal[ffIndex] = 2;
		  else if (faultVal1)
		    fStateVal[ffIndex] = 1;
		  else
		    fStateVal[ffIndex] = 0;
		  fStateNext[ffIndex] = stfList[fIndex];
		  stfList[fIndex] = ffIndex;
		  stfCount[fIndex]++;
	      }
	    }
	  }	// for (j...)
	}	// if (currId...)
    }	// for (i...)
}

////////////////////////////////////////////////////////////////////////
// gateLevelCkt class
////////////////////////////////////////////////////////////////////////

// constructor: reads in the *.lev file for the gate-level ckt
gateLevelCkt::gateLevelCkt(char *cktName)
{
    FILE *yyin;
    char fName[256];
    int i, j, count;
    char c;
    int netnum, junk;
    int f1, f2, f3;
    int levelSize[MAXlevels];

    strcpy(fName, cktName);
    strcat(fName, ".lev");
    yyin = fopen(fName, "r");
    if (yyin == NULL)
    {
	fprintf(stderr, "Can't open .lev file\n");
	exit(-1);
    }

    numpri = numgates = numout = maxlevels = numff = 0;
    maxLevelSize = 32;
    for (i=0; i<MAXlevels; i++)
	levelSize[i] = 0;

    fscanf(yyin, "%d", &count);	// number of gates
    fscanf(yyin, "%d", &junk);

    // allocate space for gates
    gtype = new unsigned char[count+64];
    fanin = new short[count+64];
    fanout = new short[count+64];
    levelNum = new int[count+64];
    po = new unsigned[count+64];
    inlist = new int * [count+64];
    fnlist = new int * [count+64];
    sched = new char[count+64];
    inFFsched = new char[count];	// only for the ff's
    currId = new int[count+64];
    value1 = new unsigned int[count+64];
    value2 = new unsigned int[count+64];
    fvalue1 = new unsigned int[count+64];
    fvalue2 = new unsigned int[count+64];

    // now read in the circuit
    numTieNodes = 0;
    for (i=1; i<count; i++)
    {
	fscanf(yyin, "%d", &netnum);
	fscanf(yyin, "%d", &f1);
	fscanf(yyin, "%d", &f2);
	fscanf(yyin, "%d", &f3);

	numgates++;
	gtype[netnum] = (unsigned char) f1;
	f2 = (int) f2;
	levelNum[netnum] = f2;
	levelSize[f2]++;

	if (f2 >= (maxlevels))
	    maxlevels = f2 + 5;
	if (maxlevels > MAXlevels)
	{
	    fprintf(stderr, "MAXIMUM level (%d) exceeded.\n", maxlevels);
	    exit(-1);
	}

	fanin[netnum] = (int) f3;
	if (f3 > MAXFanout)
	    fprintf(stderr, "Fanin count (%d) exceeded\n", fanin[netnum]);

	if (gtype[netnum] == T_input)
	{
	    inputs[numpri] = netnum;
	    numpri++;
	}
	if (gtype[netnum] == T_dff)
	{
	    if (numff >= (MAXFFS-1))
	    {
		fprintf(stderr, "The circuit has more than %d FFs\n", MAXFFS-1);
		exit(-1);
	    }
	    ff_list[numff] = netnum;
	    numff++;
	}

	sched[netnum] = inFFsched[netnum] = 0;

	// now read in the faninlist
	inlist[netnum] = new int[fanin[netnum]];
	for (j=0; j<fanin[netnum]; j++)
	{
	    fscanf(yyin, "%d", &f1);
	    inlist[netnum][j] = (int) f1;
	}

	for (j=0; j<fanin[netnum]; j++)	  // followed by close to samethings
	    fscanf(yyin, "%d", &junk);

	// read in the fanout list
	fscanf(yyin, "%d", &f1);
	fanout[netnum] = (int) f1;
	currId[netnum] = 0;

	if (gtype[netnum] == T_output)
	{
	    po[netnum] = TRUE;
	    outputs[numout] = netnum;
	    numout++;
	}
	else
	    po[netnum] = 0;

	if (fanout[netnum] > MAXFanout)
	    fprintf(stderr, "Fanout count (%d) exceeded\n", fanout[netnum]);

	fnlist[netnum] = new int[fanout[netnum]];
	for (j=0; j<fanout[netnum]; j++)
	{
	    fscanf(yyin, "%d", &f1);
	    fnlist[netnum][j] = (int) f1;
	}

	if (gtype[netnum] == T_tie1)
        {
	    TIES[numTieNodes] = netnum;
	    numTieNodes++;
	    if (numTieNodes > 511)
	    {
		fprintf(stderr, "Can't handle more than 512 tied nodes\n");
		exit(-1);
	    }
            value1[netnum] = fvalue1[netnum] = ALLONES;
            value2[netnum] = fvalue2[netnum] = ALLONES;
        }
        else if (gtype[netnum] == T_tie0)
        {
	    TIES[numTieNodes] = netnum;
	    numTieNodes++;
	    if (numTieNodes > 511)
	    {
		fprintf(stderr, "Can't handle more than 512 tied nodes\n");
		exit(-1);
	    }
            value1[netnum] = fvalue1[netnum] = 0;
            value2[netnum] = fvalue2[netnum] = 0;
        }
        else
        {
	    // assign all values to unknown
	    value1[netnum] = fvalue1[netnum] = 0;
	    value2[netnum] = fvalue2[netnum] = ALLONES;
	}

	// read in and discard the observability values
/*
	fscanf(yyin, "%d", &junk);
	fscanf(yyin, "%c", &c);		// ';' or 'O'
	fscanf(yyin, "%d", &junk);	// 0 controllability
	fscanf(yyin, "%d", &junk);	// 1 controllability
*/

	// read till end of line
	while ((c = getc(yyin)) != '\n' && c != EOF)
	    ;
    }	// for (i...)
    fclose(yyin);
    numgates++;
    numFaultFreeGates = numgates;
    numFaultsInserted = 0;

    // now compute the maximum width of the level
    for (i=0; i<maxlevels; i++)
    {
	if (levelSize[i] > maxLevelSize)
	    maxLevelSize = levelSize[i] + 1;
    }

    // initialize all nodes to unknowns
    for (i=numgates; i<numgates+64; i++)
    {
	value1[i] = fvalue1[i] = 0;
	value2[i] = fvalue2[i] = ALLONES;
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

    printf("Successfully read in circuit:\n");
    printf("\t%d PIs.\n", numpri);
    printf("\t%d POs.\n", numout);
    printf("\t%d Dffs.\n", numff);
    printf("\t%d total number of gates\n", numFaultFreeGates);
    printf("\t%d levels in the circuit.\n", maxlevels / 5);
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
    setupFaults(cktName);

    totalNumFaults = numFaults;	// numFaults = number of faults left
    printf("\t%d total number of faults\n", totalNumFaults);

    setupfState();	// set up the fault list for faults with FF's

    if (INIT0)	// if start from a initial state
    {
	RESET_FF1 = new unsigned int [numff+2];
	RESET_FF2 = new unsigned int [numff+2];

	strcpy(fName, cktName);
	strcat(fName, ".initState");
	yyin = fopen(fName, "r");
	if (yyin == NULL)
	{	fprintf(stderr, "Can't open %s\n", fName); exit(-1);}

	for (i=0; i<numff; i++)
	{
	    fscanf(yyin, "%c", &c);

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

	fclose(yyin);
    }
}

////////////////////////////////////////////////////////////////////////
// clearFF()
//	This function resets all fvalues of FF's to the good ckt's values
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::clearFF(int faulty)
{
    int i, ffNum;

    if (faulty)		// faulty ckt
    {
      for (i=0; i<actFFLen; i++)
      {
	ffNum = actFFList[i];
	inFFsched[ffNum] = 0;
	sched[ffNum] = 0;		// reset the schedule bit
	fvalue1[ffNum] = value1[ffNum];	// reset the faulty
	fvalue2[ffNum] = value2[ffNum];	// value to be good
      }
    }
    else		// good ckt
    {
      for (i=0; i<actFFLen; i++)
      {
	ffNum = actFFList[i];
	fvalue1[ffNum] = value1[ffNum];	// reset the faulty
	fvalue2[ffNum] = value2[ffNum];	// value to be good
      }
    }
}

////////////////////////////////////////////////////////////////////////
// setStateFaults()
//	This function sets up the faults at the FF's for the current
// time frame.
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::setStateFaults()
{
    int n, i, j;
    int stfLoc, fStateNode, successor;
    int fscaleIndex;
    unsigned int faultyVal;

    actFFLen = 0;
    // now insert the affected ff's by this fault if any
    for (n = 0; n < numFaultsInserted; n++)
    {
      fscaleIndex = fIndices[n];
      for (i=0; i<stfCount[fscaleIndex]; i++)
      {
	faultyVal = (1 << (n+1));
	stfLoc = stfList[fscaleIndex];
	fStateNode = fStateFF[stfLoc];
	if (inFFsched[fStateNode] == 0)
	{
	    actFFList[actFFLen] = fStateNode;
	    actFFLen++;
	    inFFsched[fStateNode] = 1;
	}

	if (fStateVal[stfLoc] == 0)	// 0
	{
	    faultyVal = ALLONES ^ faultyVal;
	    fvalue1[fStateNode] &= faultyVal;
	    fvalue2[fStateNode] &= faultyVal;
	}
	else if (fStateVal[stfLoc] == 1)	// 1
	{
	    fvalue1[fStateNode] |= faultyVal;
	    fvalue2[fStateNode] |= faultyVal;
	}
	else				// X
	{
	    fvalue2[fStateNode] |= faultyVal;
	    faultyVal = ALLONES ^ faultyVal;
	    fvalue1[fStateNode] &= faultyVal;
 	}

	currId[fStateNode] = faultId;

	// schedule the successor to the event queue
	for (j=0; j<fanout[fStateFF[stfLoc]]; j++)
	{
	    successor = fnlist[fStateNode][j];
	    if (sched[successor] == 0)
	    {
		insertEvent(levelNum[successor], successor);
		sched[successor] = 1;
	    }
	}
	stfList[fscaleIndex] = fStateNext[stfLoc];
	fStateFreeLoc(stfLoc);
      }	// for (i...)
      stfCount[fscaleIndex] = 0;
    }	// for (n...)
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
// faultActive()
//	This function returns one if the fault will be active for the
// current vector (given already the good ckt values around the fault)
////////////////////////////////////////////////////////////////////////

int gateLevelCkt::faultActive(int gateN, int ioN, int stuckN)
{
    int successor, predecessor;
    int i, j, propagate;

    if (po[gateN])
    {
	if ((stuckN && !value2[gateN]) || (!stuckN && value1[gateN]))
	    return(1);
	else
	    return(0);
    }

    propagate = 1;
    if (ioN > 0)	// fault at input of gate
    {
	// no propagation if stuckN and goodckt value are the same
	if ((stuckN == 0) && (value2[inlist[gateN][ioN-1]] == 0) || 
	    (stuckN && value1[inlist[gateN][ioN-1]]))
	    return(0);

	i = 0;
	switch (gtype[gateN])
	{
	    case T_and:
	    case T_nand:
		while ((i < fanin[gateN]) && (propagate))
		{
		    if (ioN != (i+1))
		    {
		    	predecessor = inlist[gateN][i];
			// no propagate if one of the inputs is zero
			if ((value1[predecessor] == 0) &&
			    (value2[predecessor] == 0))
				propagate = 0;
		    }
		    i++;
		}	// while (i...)
		break;
	    case T_or:
	    case T_nor:
		while ((i < fanin[gateN]) && (propagate))
		{
		    if (ioN != (i+1))
		    {
		    	predecessor = inlist[gateN][i];
			// no propagate if one of the inputs is zero
			if ((value1[predecessor]) && (value2[predecessor]))
				propagate = 0;
		    }
		    i++;
		}	// while (i...)
		break;
	    default:
		break;
	}	// switch
    }	// if fault on input
    else	// fault stuck at output
    {
	// check if stuck value the same as good value
	if ((stuckN == 0) && (value2[gateN] == 0) || (stuckN && value1[gateN]))
	    return(0);
    }

    // second level and also for faults at output of gates
    if ((propagate) && (gtype[gateN] != T_dff))
    {
	i = 0;
	propagate = 0;
	while ((i < fanout[gateN]) && (!propagate))
	{
	    successor = fnlist[gateN][i];
	    switch (gtype[successor])
	    {
		case T_and:
		case T_nand:
		    // check if all the inputs besides the connection
		    j = 0;
		    propagate = 1;
		    while ((j < fanin[successor]) && (propagate))
		    {
			predecessor = inlist[successor][j];
			if (predecessor != gateN)
			{
			    // no propagate if one of the inputs is zero
			    if ((value1[predecessor] == 0) &&
				(value2[predecessor] == 0))
				propagate = 0;
			}
			j++;
		    }	// while (j..)
		    break;
		case T_or:
		case T_nor:
		    // check if all the inputs besides the connection
		    j = 0;
		    propagate = 1;
		    while ((j < fanin[successor]) && (propagate))
		    {
			predecessor = inlist[successor][j];
			if (predecessor != gateN)
			{
			    // no propagate if one of the inputs is one
			    if ((value1[predecessor]) && (value2[predecessor]))
				propagate = 0;
			}
			j++;
		    }	// while (j...)
		    break;
		default:
		    propagate = 1;
	    }
	    i++;
	}
    }	// if (propagate)

    return(propagate);
}

////////////////////////////////////////////////////////////////////////
// insertFaults()
// 	This function insert numF number of faults into the circuit.  It
// makes sure that all numF faults are excited by the current good input
// vector.
// Returns the index of the next fault.
////////////////////////////////////////////////////////////////////////
int gateLevelCkt::insertFaults(int numF, int fscaleIndex)
{
    int n, i;			// counters
    int successor, predecessor;		// predecessor & successor gate
    int gateN, ioN, stuckN;
    unsigned int faultyVal;
    int tmpLevel;

    n = 0;
    faultId++;		// increment the fault id
    if (numF > 31)
    {
	fprintf(stderr, "Can't insert more than 31 faults.\n");
	exit(-1);
    }

    while ((n < numF) && (fscaleIndex != -1))
    {
      gateN = fGate[fscaleIndex];
      ioN = fIO[fscaleIndex];
      stuckN = fStuck[fscaleIndex];
      if ((stfCount[fscaleIndex] > 0) || faultActive(gateN, ioN, stuckN))
      {
	if (fstatus[fscaleIndex] == 0)
	   fstatus[fscaleIndex] = EXCITED_1_LEVEL;
	numExcited[fscaleIndex]++;

      	fgates[numFaultsInserted] = numgates;
      	fIndices[numFaultsInserted] = fscaleIndex;

      	currId[numgates+1] = faultId;

      	if ((ioN == 0) && (po[gateN] == 0))	// fault on output of the gate
      	{					// but not primary output
	  successor = fnlist[gateN][0];
	  if (levelNum[successor] == levelNum[gateN] + 1)
	  {	// already a fault inserted at this net
	    levelNum[numgates] = levelNum[gateN] + 2;
	    gateN = successor;
	  }
	  else	// first fault at this net
	    levelNum[numgates] = levelNum[gateN] + 1;

    	  fanout[numgates] = fanout[gateN];
	  inlist[numgates][1] = gateN;
	  if (stuckN == 1)
	  {
	    gtype[numgates] = T_or;
	    faultyVal = (1 << (n+1));
	    fvalue1[numgates+1] = fvalue2[numgates+1] = faultyVal;
	  }
	  else
	  {
	    gtype[numgates] = T_and;
	    faultyVal = ALLONES ^ (1 << (n+1));
	    fvalue1[numgates+1] = fvalue2[numgates+1] = faultyVal;
	  }
	  value1[numgates] = value1[gateN];	// set the good value of gate
	  value2[numgates] = value2[gateN];

	  for (i=0; i<fanout[gateN]; i++)
	  {
	    successor = fnlist[numgates][i] = fnlist[gateN][i];
	    if (successor < numFaultFreeGates)	// successor not faulty
		predOfSuccInput[numgates][i] = predOfSuccInput[gateN][i];
	    else	// successor a faulty gate
		predOfSuccInput[numgates][i] = 1;

	    inlist[successor][predOfSuccInput[numgates][i]] = numgates;
	  }
          succOfPredOutput[numgates][1] = 0;
	  if (gateN >= numFaultFreeGates) // double faults, gateN is 1st fault
	    predOfSuccInput[gateN][0] = 1;

	  fanout[gateN] = 1;
	  fnlist[gateN][0] = numgates;
	}
      	else	// fault on input of gate
        {
	  if ((po[gateN] == 1) && (ioN == 0))	//primary output
	   ioN = 1;
	  predecessor = inlist[numgates][1] = inlist[gateN][ioN-1];

	  //check if fault on input of DFF
	  if (gtype[gateN] == T_dff)
	    tmpLevel = maxlevels;
	  else
	    tmpLevel = levelNum[gateN];
	  if (levelNum[predecessor] == tmpLevel - 1)
	  {	// already a fault inserted at this net
	    gateN = predecessor;
	    levelNum[numgates] = tmpLevel - 2;
	    predecessor = inlist[numgates][1] = inlist[gateN][1];
	    ioN = 2;
	  }
	  else
	  {
	    levelNum[numgates] = tmpLevel - 1;
	  }	// else

	  fanout[numgates] = 1;

          succOfPredOutput[numgates][1] = succOfPredOutput[gateN][ioN-1];
	  predOfSuccInput[numgates][0] = ioN-1;
	  if (gateN >= numFaultFreeGates)	// double faults at a site, 
            succOfPredOutput[gateN][1] = 0;
	  if (predecessor >= numFaultFreeGates)
	    predOfSuccInput[predecessor][succOfPredOutput[numgates][1]] = 1;
	  fnlist[predecessor][succOfPredOutput[numgates][1]] = numgates;

	  if (stuckN == 1)
	  {
	    gtype[numgates] = T_or;
	    faultyVal = (1 << (n+1));
	    fvalue1[numgates+1] = fvalue2[numgates+1] = faultyVal;
	  }
	  else	// stuck at 0
	  {
	    gtype[numgates] = T_and;
	    faultyVal = ALLONES ^ (1 << (n+1));
	    fvalue1[numgates+1]= fvalue2[numgates+1] = faultyVal;
	  }
	  value1[numgates] = value1[predecessor]; // set the good value of gate
	  value2[numgates] = value2[predecessor];

	  fnlist[numgates][0] = gateN;
	  inlist[gateN][ioN-1] = numgates;
        }	// else fault on input

        // insert the fault into the time wheel
	fvalue1[numgates] = 0;
	fvalue2[numgates] = ALLONES;
       	insertEvent(levelNum[numgates], numgates);
      	sched[numgates] = 1;

        n++;
        numgates += 2;
        numFaultsInserted++;
      }		// if (stfCount > 0 || faultActive)

      fscaleIndex = fsuclist[fscaleIndex];
    }	// outer while (n...)

    return(fscaleIndex);
}

////////////////////////////////////////////////////////////////////////
// removeFaults()
//	gateN = gate number for the inserted faulty gate
////////////////////////////////////////////////////////////////////////
void gateLevelCkt::removeFaults()
{
    int gateN;
    int predecessor;
    int successor;
    int fIndex, fio, fgate;
    int n, j;			// counters

    for (n = 0; n < numFaultsInserted; n++)
    {
      numgates -= 2;
      fIndex = fIndices[n];

      gateN = fgates[n];
      predecessor = inlist[gateN][1];	// on the inserted fault gate

      fio = fIO[fIndex];
      fgate = fGate[fIndex];

      // reconnections if fault at input of the gate
      if ((fio > 0) || (po[fgate] == 1))
      {
	successor = fnlist[gateN][0];
	fnlist[predecessor][succOfPredOutput[gateN][1]] = successor;
	inlist[successor][predOfSuccInput[gateN][0]] = predecessor;
	if (predecessor >= numFaultFreeGates)
	    predOfSuccInput[predecessor][succOfPredOutput[gateN][1]] = 
		predOfSuccInput[gateN][0];
      }
      else	// fault inserted at output
      {
        fanout[predecessor] += (fanout[gateN]-1);
	for (j=0; j<fanout[predecessor]; j++)
	{
	    successor = fnlist[gateN][j];
	    fnlist[predecessor][j] = successor;
	    inlist[successor][predOfSuccInput[gateN][j]] = predecessor;
	}	// for (j...)
      }	// else
    }	// for (n...)
    numFaultsInserted = 0;
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
printf("Initialize circuit to values in *.initState!\n");
	for (i=0; i<numff; i++)
	{
	    value1[ff_list[i]] = value2[ff_list[i]] = RESET_FF1[i];
	    fvalue1[ff_list[i]] = fvalue2[ff_list[i]] = RESET_FF2[i];

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
	    fvalue1[predecessor] = fvalue2[predecessor] = RESET_FF2[i];

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
		value1[inputs[i]] = fvalue1[inputs[i]] = 0;
		value2[inputs[i]] = fvalue2[inputs[i]] = 0;
		break;
	    case '1':
		value1[inputs[i]] = fvalue1[inputs[i]] = ALLONES;
		value2[inputs[i]] = fvalue2[inputs[i]] = ALLONES;
		break;
	    case 'x':
	    case 'X':
		value1[inputs[i]] = fvalue1[inputs[i]] = 0;
		value2[inputs[i]] = fvalue2[inputs[i]] = ALLONES;
		break;
	    default:
		fprintf(stderr, "%c: error in the input vector.\n", vec[i]);
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
		fprintf(stderr, "illegal gate type1 %d %d\n", gateN, gtype[gateN]);
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
// globalFsim() -
//	Global fault simulate (faults have been inserted.)
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::globalFsim()
{
    int sucLevel;
    int gateN, predecessor, successor;
    int *predList;
    int i, mask;
    unsigned int val1, val2, tmpVal;
    unsigned goodVal1, goodVal2, faultVal1, faultVal2;

    currLevel = 1;	// start from everyone node after PI's and dff's
    actLen = 0;
    while (currLevel < maxlevels)
    {
    	gateN = retrieveEvent();
	if (gateN != -1)	// if a valid event
	{
	    sched[gateN] = 0;
    	    switch (gtype[gateN])
    	    {
	      case T_and:
    	    	val1 = val2 = ALLONES;
		predList = inlist[gateN];
    	    	for (i=0; i<fanin[gateN]; i++)
    	    	{
		    predecessor = predList[i];
		    if (currId[predecessor] == faultId)
		    {
		    	val1 &= fvalue1[predecessor];
		    	val2 &= fvalue2[predecessor];
		    }
		    else
		    {
		    	val1 &= value1[predecessor];
		    	val2 &= value2[predecessor];
		    }
    	    	}
	    	break;
	      case T_nand:
    	        val1 = val2 = ALLONES;
		predList = inlist[gateN];
    	    	for (i=0; i<fanin[gateN]; i++)
    	    	{
		    predecessor = predList[i];
		    if (currId[predecessor] == faultId)
		    {
		    	val1 &= fvalue1[predecessor];
		    	val2 &= fvalue2[predecessor];
		    }
		    else
		    {
		    	val1 &= value1[predecessor];
		    	val2 &= value2[predecessor];
		    }
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
		    if (currId[predecessor] == faultId)
		    {
		    	val1 |= fvalue1[predecessor];
		    	val2 |= fvalue2[predecessor];
		    }
		    else
		    {
		    	val1 |= value1[predecessor];
		    	val2 |= value2[predecessor];
		    }
    	    	}
	    	break;
	      case T_nor:
    	    	val1 = val2 = 0;
		predList = inlist[gateN];
		for (i=0; i<fanin[gateN]; i++)
    	    	{
		    predecessor = predList[i];
		    if (currId[predecessor] == faultId)
		    {
		    	val1 |= fvalue1[predecessor];
		    	val2 |= fvalue2[predecessor];
		    }
		    else
		    {
		    	val1 |= value1[predecessor];
		    	val2 |= value2[predecessor];
		    }
    	    	}
	    	tmpVal = val1;
	    	val1 = ALLONES ^ val2;
	    	val2 = ALLONES ^ tmpVal;
	    	break;
	      case T_not:
	    	predecessor = inlist[gateN][0];
		if (currId[predecessor] == faultId)
		{
	    	    val1 = ALLONES ^ fvalue2[predecessor];
	    	    val2 = ALLONES ^ fvalue1[predecessor];
		}
		else
		{
	    	    val1 = ALLONES ^ value2[predecessor];
	    	    val2 = ALLONES ^ value1[predecessor];
		}
	    	break;
	      case T_buf:
	      case T_dff:
	    	predecessor = inlist[gateN][0];
		if (currId[predecessor] == faultId)
		{
	    	    val1 = fvalue1[predecessor];
	    	    val2 = fvalue2[predecessor];
		}
		else
		{
	    	    val1 = value1[predecessor];
	    	    val2 = value2[predecessor];
		}
	    	break;
	      case T_xor:
		predList = inlist[gateN];
		if (currId[predList[0]] == faultId)
		{
	    	    val1 = fvalue1[predList[0]];
	    	    val2 = fvalue2[predList[0]];
		}
		else
		{
	    	    val1 = value1[predList[0]];
	    	    val2 = value2[predList[0]];
		}

            	for(i=1;i<fanin[gateN];i++)
            	{
	    	    predecessor = predList[i];
		    if (currId[predecessor] == faultId)
		    {
                      tmpVal = ALLONES^(((ALLONES^fvalue1[predecessor]) &
              		(ALLONES^val1)) | (fvalue2[predecessor]&val2));
                      val2 = ((ALLONES^fvalue1[predecessor]) & val2) |
                  	(fvalue2[predecessor] & (ALLONES^val1));
		    }
		    else
		    {
                      tmpVal = ALLONES^(((ALLONES^value1[predecessor]) &
              		(ALLONES^val1)) | (value2[predecessor]&val2));
                      val2 = ((ALLONES^value1[predecessor]) & val2) |
                  	(value2[predecessor] & (ALLONES^val1));
		    }
                    val1 = tmpVal;
            	}
	    	break;
	      case T_xnor:
		predList = inlist[gateN];
		if (currId[predList[0]] == faultId)
		{
		    val1 = fvalue1[predList[0]];
	    	    val2 = fvalue2[predList[0]];
		}
		else
		{
		    val1 = value1[predList[0]];
	    	    val2 = value2[predList[0]];
		}

            	for(i=1;i<fanin[gateN];i++)
            	{
	    	    predecessor = predList[i];
		    if (currId[predecessor] == faultId)
		    {
                      tmpVal = ALLONES^(((ALLONES^fvalue1[predecessor]) &
                   	(ALLONES^val1)) | (fvalue2[predecessor]&val2));
                      val2 = ((ALLONES^fvalue1[predecessor]) & val2) |
                  	(fvalue2[predecessor]& (ALLONES^val1));
		    }
		    else
		    {
                      tmpVal = ALLONES^(((ALLONES^value1[predecessor]) &
                   	(ALLONES^val1)) | (value2[predecessor]&val2));
                      val2 = ((ALLONES^value1[predecessor]) & val2) |
                  	(value2[predecessor]& (ALLONES^val1));
		    }
                    val1 = tmpVal;
            	}
	    	tmpVal = val1;
		val1 = ALLONES ^ val2;
	    	val2 = ALLONES ^ tmpVal;
	    	break;
	      case T_output:
		predecessor = inlist[gateN][0];
		if (currId[predecessor] == faultId)
		{
	    	    val1 = fvalue1[predecessor];
	    	    val2 = fvalue2[predecessor];
		}
		else
		{
	    	    val1 = value1[predecessor];
	    	    val2 = value2[predecessor];
		}
	        break;
	      case T_input:
	      case T_tie0:
	      case T_tie1:
	      case T_tieX:
	      case T_tieZ:
		if (currId[gateN] == faultId)
		{
	    	    val1 = fvalue1[gateN];
	    	    val2 = fvalue2[gateN];
		}
		else
		{
	    	    val1 = value1[gateN];
	    	    val2 = value2[gateN];
		}
	    	break;
	      default:
		fprintf(stderr, "illegal gate type3 %d\n", gateN);
		exit(-1);
    	    }	// switch

	    // if gate value changed
    	    if ((val1 != value1[gateN]) || (val2 != value2[gateN]))
	    {
		currId[gateN] = faultId;
		fvalue1[gateN] = val1;
		fvalue2[gateN] = val2;

		if (po[gateN])		// if primary output node
		{
		    mask = 1;
		    //goodVal1 = val1 & mask;
		    //goodVal2 = val2 & mask;
		    goodVal1 = value1[gateN] & mask;
		    goodVal2 = value2[gateN] & mask;

		    if ((goodVal1) || (!goodVal2)) // good value not don't care
		    {
			for (i=0; i < numFaultsInserted; i++)
			{
	    		    mask = mask << 1;
	    		    faultVal1 = val1 & mask;
    	    		    faultVal2 = val2 & mask;

	    		    if ((!faultVal1) && (faultVal2)) // fault don't care
	    		    {				// potential detect
				if (fstatus[fIndices[i]] == 0)
				{
			    	    fstatus[fIndices[i]] = POTENTIAL;
if (OUTPUT_DET)
printf("Potential detect output #%d, fault #%d\n", gateN, i);
				}
	    		    }
	    		    else if (((goodVal1) && (!faultVal1)) || 
				((!goodVal1) && (faultVal1)))
	    		    {				// detect
			    	if (fstatus[fIndices[i]] < LOW_DETECT)
{
int fIndex = fIndices[i];
if (OUTPUT_DET)
{
	std::ofstream log ("log.txt", std::ios::app);
	log << fGate[fIndex] ;
	log << " ";
	log << fGate[fIndex] ;
	log << " ";
	log << fStuck[fIndex]; 
	log << "\r\n" ;
	log.close();
printf("%d %d %d;\n", fGate[fIndex], fIO[fIndex], fStuck[fIndex]);
}
if (!NO_FAULTDROP)
{
    dropFault(fIndex);
    fstatus[fIndex] = LOW_DETECT;
}
else
tmpDrop(fIndex);
}
			    }
		    	}	// for (i...)
    		    }	// if (!((good...)
		}	// if (po...)

		// schedule the successors into timewheel
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

    // now insert the faulty effects from the activation list for the FF's
    if (actLen > 0)
    {
	saveStateFaults();
    }
}

////////////////////////////////////////////////////////////////////////
// lowFault member functions
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::setupFaults(char *cktName)
{
    FILE *fFile;
    char fFileName[256];
    char **poFaults;
    int i, preIndex, sucIndex;
    char thisChar;

    strcpy(fFileName, cktName);
    strcat(fFileName, ".eqf");
    fFile = fopen(fFileName, "r");
    if (fFile == NULL)
    {
	fprintf(stderr, "Can't open fault file %s\n", fFileName);
	exit(-1);
    }

    // first count the number of faults
    numFaults = 0;
    while ((thisChar = getc(fFile)) != EOF)
    {
	if (thisChar == ';')
	    numFaults++;
    }
    fclose(fFile);

    // now allocate space and read in the faults
    fGate = new int[numFaults];
    fIO = new int[numFaults];
    fStuck = new int[numFaults];
    fprelist = new int[numFaults];
    fsuclist = new int[numFaults];
    fstatus = new int[numFaults];
    numExcited = new int[numFaults];
    no_faultdrop_fstatus = new int[numFaults];
// the next two lines are for storing tmp drop faults
tmpDropTable = new int[numFaults];
tmpDropDet = new int[numFaults];
tmpDropIndex = 0;
    stfList = new int [numFaults];
    stfCount = new int [numFaults];

    // to avoid quad faults at PO's
    poFaults = new char *[numgates];
    for (i=0; i<numgates; i++)
    {
	poFaults[i] = new char[2];
	poFaults[i][0] = poFaults[i][1] = 0;
    }

    startfLoc = 0;
    faultId = 0;

    fFile = fopen(fFileName, "r");

    for (i=0; i<numFaults; i++)
    {
    	fscanf(fFile, "%d", &fGate[i]);
    	fscanf(fFile, "%d", &fIO[i]);
    	fscanf(fFile, "%d", &fStuck[i]);
	fstatus[i] = no_faultdrop_fstatus[i] = 0;
	numExcited[i] = 0;

	fprelist[i] = i-1;	// initialize the links
	fsuclist[i] = i+1;

	stfCount[i] = 0;

	// now to avoid counting too many faults at PO's
	if (gtype[fGate[i]] == T_output)
	{
	    if (poFaults[fGate[i]][fStuck[i]] == 1)
	    {
		i--;
		numFaults--;
	    }
	    else
		poFaults[fGate[i]][fStuck[i]] = 1;
	}

	// pass the rest of line
	while ((thisChar = getc(fFile)) != ';')
	    ;
	while ((thisChar = getc(fFile)) != RETURN)
	{
	    if (thisChar == 'R')        // redundant fault
		fstatus[i] = REDUNDANT;
	}
    }
    // boundary conditions
    fprelist[0] = -1;
    fsuclist[numFaults-1] = -1;

    fclose(fFile);

    // now drop RED faults
    for (i=0; i<numFaults; i++)
    {
	if (fstatus[i] == REDUNDANT)
	{
	    preIndex = fprelist[i];
	    sucIndex = fsuclist[i];

	    if (sucIndex != -1)
	    	fprelist[sucIndex] = preIndex;

	    if (preIndex != -1)
	    	fsuclist[preIndex] = sucIndex;
	    else
		startfLoc = sucIndex;
	}
    }
}

////////////////////////////////////////////////////////////////////////
// dropFault()
//	This function removes a fault from the lowFault class and
// rebuilds the links of the fault list.
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::dropFault(int fIndex)
{
    int preIndex = fprelist[fIndex];
    int sucIndex = fsuclist[fIndex];

    numFaults--;
    if (sucIndex != -1)
    	fprelist[sucIndex] = preIndex;

    if (preIndex != -1)
    	fsuclist[preIndex] = sucIndex;
    else
	startfLoc = sucIndex;
}

////////////////////////////////////////////////////////////////////////
// tmpDrop()
//	This function drops the fault temporarily.
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::tmpDrop(int fIndex)
{
    tmpDropTable[tmpDropIndex] = fIndex;
    tmpDropDet[tmpDropIndex] = fstatus[fIndex];
    fstatus[fIndex] = LOW_DETECT;
    tmpDropIndex++;
    no_faultdrop_fstatus[fIndex]++;
}

////////////////////////////////////////////////////////////////////////
// restoreDrop()
//	This function restores the temporarily dropped faults
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::restoreDrop()
{
    int i;

    for (i=0; i<tmpDropIndex; i++)
	fstatus[tmpDropTable[i]] = tmpDropDet[i];

    tmpDropIndex = 0;
}

////////////////////////////////////////////////////////////////////////
// observeOutputs()
//	This function prints the outputs of the fault-free circuit.
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::observeOutputs()
{
    int i;

    printf("\t");
    for (i=0; i<numout; i++)
    {
	if (value1[outputs[i]] && value2[outputs[i]])
	    printf("1");
	else if ((value1[outputs[i]] == 0) && (value2[outputs[i]] == 0))
	    printf("0");
	else
	    printf("X");
    }

    printf("\n");
    for (i=0; i<numff; i++)
    {
	if (value1[ff_list[i]] && value2[ff_list[i]])
	    printf("1");
	else if ((value1[ff_list[i]] == 0) && (value2[ff_list[i]] == 0))
	    printf("0");
	else
	    printf("X");
    }

    printf("\n");
}

void gateLevelCkt::printGoodSig(FILE *sigFile, int v)
{
    int i;

    fprintf(sigFile, "vec %d: ", v);
    for (i=0; i<numout; i++)
    {
	if (value1[outputs[i]] && value2[outputs[i]])
	    fprintf(sigFile, "1");
	else if ((value1[outputs[i]] == 0) && (value2[outputs[i]] == 0))
	    fprintf(sigFile, "0");
	else
	    fprintf(sigFile, "X");
    }

    fprintf(sigFile, "\n");
}

void gateLevelCkt::printSig(FILE *sigFile)
{
    int i, j;
    unsigned int val1, val2, mask = 1;

  for (j=0; j<numFaultsInserted; j++)
  {
    mask = mask << 1;
    if (fstatus[fIndices[j]] >= LOW_DETECT)
    {
	fprintf(sigFile, "%d: ", fIndices[j]);

	// now print the faulty output
      for (i=0; i<numout; i++)
      {
	if (currId[outputs[i]] == faultId)
	{
	  val1 = fvalue1[outputs[i]] & mask;
	  val2 = fvalue2[outputs[i]] & mask;
	}
	else
	{
	  val1 = value1[outputs[i]] & mask;
	  val2 = value2[outputs[i]] & mask;
	}

	if (val1 && val2)
	    fprintf(sigFile, "1");
	else if ((val1 == 0) && (val2 == 0))
	    fprintf(sigFile, "0");
	else
	    fprintf(sigFile, "X");
      }	 // for (i...)
      fprintf(sigFile, "\n");
    }	// if (fstatus...)
  }	// for (j...)
}

