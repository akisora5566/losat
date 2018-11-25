//
// Created by akisora5566 on 11/24/18.
//
#ifndef LOSAT_GATELEVELCKT_H
#define LOSAT_GATELEVELCKT_H

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

using namespace std;

//////////////////////////////////////////////////////////////////////////
// Global Variables

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
private:
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
    int INIT0;
    int numTieNodes;
    int TIES[512];

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
    unsigned int *RESET_FF1;	// value of reset ffs read from *.initState
    unsigned int *RESET_FF2;	// value of reset ffs read from *.initState

    gateLevelCkt(string, int);	// constructor
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

#endif //LOSAT_GATELEVELCKT_H
