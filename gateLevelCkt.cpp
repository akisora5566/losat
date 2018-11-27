#include <cstdlib>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <cmath>
#include "gateLevelCkt.h"
using namespace std;


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
gateLevelCkt::gateLevelCkt(string cktName, int INIT0)
{
    ifstream yyin;
    string fName;
    int i, j, count;
    char c;
    int netnum, junk;
    int f1, f2, f3;
    int levelSize[MAXlevels];
    this->INIT0 = INIT0;

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

    // initializing stored ff values
    value1_ffs = new unsigned int[count+64];
    value2_ffs = new unsigned int[count+64];
    id_unknown_ffs = new int[count+64];
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
        int reset_ff_count = 0;
        for (int j = 0; j < numff; j++)
        {
            if ((value1_ffs[j] && value2_ffs[j])       // normal 1
                || (!value1_ffs[j] && !value2_ffs[j])) // normal 0
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
        for (int j = 0; j < numff; j++)
        {
            if ((value1_ffs[j] && value2_ffs[j])       // normal 1
                || (!value1_ffs[j] && !value2_ffs[j])) // normal 0
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
        /*
        cout << "reset_power0[" << i << "] = " << reset_power0[i] << endl;
        cout << "reset_power1[" << i << "] = " << reset_power1[i] << endl;
        cout << "prob_maskto0[" << i << "] = " << prob_maskto0[i] << endl;
        cout << "prob_maskto1[" << i << "] = " << prob_maskto1[i] << endl;
        */
    }
}



////////////////////////////////////////////////////////////////////////
////********************** PUBLIC FUNCTIONS ************************////
////////////////////////////////////////////////////////////////////////

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
    //cout << "Initialize circuit to values in *.initState!\n";
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
                    id_unknown_max++;
                    id_unknown[inputs[i]] = id_unknown_max;
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
        id_unknown_ffs[i] = 0;
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
                    //printGateValue(predecessor);
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
                    //cout << "Xs from PI squashed on gate " << gateN << "!\n";
                    val1 = ALLONES;
                    val2 = ALLONES;
                    id_temp = 0;
                    break;
                }
                else if ((xi_count%2) && (xi_bar_count%2))
                {
                    // xi and xi_bar are both odd => squashed, 0
                    //cout << "Xs from PI squashed on gate " << gateN << "!\n";
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
    /*int i;

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
    cout << endl;*/
}


////////////////////////////////////////////////////////////////////////
// observeFFs()
// prints the stored value of FFs
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::observeFFs()
{
    /*int i;
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
    cout << endl;*/

}


////////////////////////////////////////////////////////////////////////
// printGateValue
// prints the value of every gates with 0,1,X and id
////////////////////////////////////////////////////////////////////////

void gateLevelCkt::printGateValue(int gateN)
{
    /*
	unsigned int val1 = value1[gateN];
    unsigned int val2 = value2[gateN];
    int id = id_unknown[gateN];

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


////////////////////////////////////////////////////////////////////////
// getValue1FFs, getValue2FFs, getIDUunknownFFs
// return all the stored value of FFs
////////////////////////////////////////////////////////////////////////

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

