# losat
Step 1: Run main.cpp on sequential circuit (main.exe -o iscas89/sXXX)
    - Comment and uncomment code to run simulation with or without reset signal masking
    -This will generate a file called output.vec
Step 2: Run my_fsim.cpp on sequential circuit(my_fsim.exe -d iscas89/sXXX)
    -When the code pauses for the first time note the number of total faults
    -Continue to get total number of faults.
