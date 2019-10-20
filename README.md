Proving Expected Sensitivity of Probabilistic Programs with Randomized Variable-Dependent Termination Time
===========
Guide to the Artifact Accompanying the Paper
--------
# Overview
This artifact contains our implementation of the RSM synthesis algorithm as described in the paper. In this work, we consider probabilistic while loops where the number of iterations is not fixed, but randomized and depends on the initial input values. The aim of this artifact is to synthesize linear RSM-maps for such programs, together with other conditions such as lipschitz continuity, bounded update, RSM-continuity and difference bound (which helps prove expected sensitivity of such programs).

# Note
The artifact currently supports simple assignment statemnts of the form x:=x+r (where r is a sampling variable with a predefined probability distribution), so it may fail for complex custom inputs as well as these programs which have no linear RSM maps.

# File Structure, Prerequisites and Compilation
## File Structure
There are several folders in the artifact. Here is a description of each:
* Inputs: contains all the twenty experimental examples (see Section 9.2 in the paper) and the related txt files about their initial valuations
* Outputs: contains two folders named "Table2" and "Table3". When an experiment is run, the final output is automatically put in "Table2" or "Table3" (which corresponds to the results in Section 9.2)
* Expected Outputs: contains the final outputs of the experimental examples as obtained by us (see Table 2 and Table 3 in Section 9.2)
* CFG: contains a C++ tool for obtaining Control Flow Graphs of programs
* Synthesis: contains a Matlab implementation for obtaining the results in our paper
* Custom: contains two folders named "Custom_inputs" and "Custom_outputs", the users can put the custom inputs into "Custom_inputs", and the final outputs will be put into "Custom_outputs" automatically by our bash scripts.

## Prerequisites
Our codes are written in C++ and Matlab. C++ is used for creating the Control-Flow Graphs of input programs and Matlab is used for the main RSM synthesis algorithm (which involves symbolic computations and an optimization problem). Moreover, we use Bash scripts for connecting the different components. Hence, we suggest that you run the artifact in a Linux environment or the provided virtual machine.
#### Complete List of Prerequisites
* g++: for compiling C++ sources
*	An environment capable of running Bash
*	Matlab (make sure that the “matlab” command is enabled in the terminal. In Ubuntu this can be done by installing the “matlab-support” package)

## Compilation
To compile the codes, open a terminal in the “Tool” folder and run this command: <br>
./compile.sh <br>
You need to have a C++-11-compliant version of g++ installed.


# Correspondence between the Paper and the Artifact
The file “example name.pdf” contains a description of the correspondence between the example programs in the artifact and those in the paper (see Table 2 and 3 in Section 9.2). 

# Obtaining the Results Reported in the Paper
To obtain all the results reported in the paper, simply open a terminal in the “Tool” folder and run the following command:<br>
./main-experiment.sh <br>
When the execution ends, all obtained results will be available in the “Outputs” folder (in "Outputs/Table2" folder or "Outputs/Table3" folder). These results match Tables 2 and 3 in our paper. <br>
To obtain Table 2 or table 3 separately, you can use the following command:<br>
./main-experiment.sh table2 <br>
or <br>
./main-experiment.sh table3 <br>





