# losat
Step 0: CHange the value of num_vectors in main.cpp according to the circuit you are testing to avoid X in final state.<br/>
Step 1: Run main.cpp on sequential circuit (main.exe -o iscas89/sXXX)<br/>
    - Comment and uncomment code to run simulation with or without reset signal masking<br/>
    -This will generate a file called output.vec<br/>
Step 2: Run my_fsim.cpp on sequential circuit(my_fsim.exe -d iscas89/sXXX)<br/>
    -When the code pauses for the first time note the number of total faults<br/>
    -Continue to get total number of faults.<br/>
