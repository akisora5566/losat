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
#include <sys/times.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include "gateLevelCkt.h"
using namespace std;


//////////////////////////////////////////////////////////////////////////
// Global Variables
char vector[5120];
int OBSERVE, INIT0;
int vecNum=0;
char* input_vectors[100000];
gateLevelCkt *circuit;


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
    vecNum = 0;
    while (moreVec)
    {
        moreVec = getVector(vecFile, vecWidth);
        if (moreVec == 1)
        {
            input_vectors[vecNum] = new char[vecWidth];
            for (int j = 0; j < vecWidth; j++)
            {
                input_vectors[vecNum][j] = vector[j];
            }
            //cout << "vector #" << vecNum << ": " << input_vectors[vecNum] << "\n";
            vecNum++;
        }
    }

    char** input_vectors_truesize;
    input_vectors_truesize = new char* [vecNum];
    //cout << "vecNum = " << vecNum << endl;
    for (int i = 0; i < vecNum; i++)
    {
        input_vectors_truesize[i] = new char[vecWidth];
        for (int j = 0; j < vecWidth; j++)
        {
            input_vectors_truesize[i][j] = input_vectors[i][j];
        }
    }
    circuit->reset();
    circuit->LogicSim(input_vectors_truesize, vecNum);

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

    circuit = new gateLevelCkt(cktName, INIT0);


    vecFile.open(vecName.c_str(), ios::in);
    if (!vecFile)
    {
        cerr << "Can't open " << vecName << "\n";
        exit(-1);
    }

    vecFile >> vecWidth;

    circuit->setTieEvents();
    totalNumVec = logicSimFromFile(vecFile, vecWidth);
    if (OBSERVE == 1)
    {
        circuit->observeOutputs();
    }
    vecFile.close();

    // output results
    end = clock();
    ut = (double) (end-start);
    cout << "Number of vectors: " << totalNumVec << "\n";
    cout << "Number of clock cycles elapsed: " << ut << "\n";
}
