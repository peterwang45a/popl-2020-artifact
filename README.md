Proving Expected Sensitivity of Probabilistic Programs with Randomized Variable-Dependent Termination Time
===========
Guide to the Artifact Accompanying the Paper
--------
# Overview
This artifact contains our implementation of the RSM synthesis algorithm as described in the paper. In this work, we consider probabilistic while loops where the number of iterations is not fixed, but randomized and depends on the initial input values. The aim of this artifact is to synthesize linear RSM-maps for such programs, together with other conditions such as lipschitz continuity, bounded update, RSM-continuity and difference bound (which helps prove expected sensitivity of such programs).

# Assumptions
We assume that the input programs have simple assignment statemnts of the form x:=x+r (where r is a sampling variable with a predefined probability distribution), so it may fail for complex custom instances as well as these programs which have no linear RSM maps.

# File Structure, Prerequisites and Compilation
## File Structure
There are several folders in the artifact. Here is a description of each:
* Inputs: contains all the twenty experimental examples (see Section 9.2 in the paper) and the related txt files about their initial valuations
* Outputs: contains two folders named "Table2" and "Table3". When an experiment is run, the final output (which is a txt file named "*_log.txt") is automatically put in "Table2" or "Table3" (which corresponds to the results in Section 9.2)
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
The file “example name.pdf” contains a description of the correspondence between the example programs in the artifact and those in the paper (see Table 2 and 3 in Section 9.2). For example, in the paper there is an example named "mini-Roulette", we call it "mini" for brevity in the artifact.

# Obtaining the Results Reported in the Paper
To obtain all the results reported in the paper, simply open a terminal in the “Tool” folder and run the following command:<br>
./main-examples.sh <br>
When the execution ends, all obtained results will be available in the “Outputs” folder (actually in "Outputs/Table2" folder or "Outputs/Table3" folder). These results match Tables 2 and 3 in our paper. <br>
<br>
To obtain Table 2 or table 3 separately, you can use the following command:<br>
./main-examples.sh table2 <br>
or <br>
./main-examples.sh table3 <br>
The obtained results will also be put in "Outputs/Table2" folder or "Outputs/Table3" folder.<br>
## Warning
The runtime may be a little longer than those in our experimental results (we made a test in the virtual machine and it was within half a minute for each example). This is because we are including a trial version of Matlab in the virtual machine (because Matlab is proprietary software). This trial version does not use multithreading and is significantly slower than a full-fledged Matlab. <br>

**Note:** You can also run the experiment on a single example program or your own custom input. See below for more details.

# Matlab License
We have pre-installed a trial version of Matlab on the Virtual Machine. If you are running the experiments on your own machine, or if Matlab asks for a license, you can obtain a free 30-day trial license at https://www.mathworks.com/campaigns/products/trials.html.

# Running on a Single Example
We have listed all the twenty experimental examples. See “example name.pdf” for the correspondence between the example names in the artifact and those in the paper. To run the experiment on a single example program, e.g. program x, you can use the following command:<br>
./run-example.sh x <br>
For example, <br>
./run-example.sh mini <br>
runs the experiment on the "mini-Roulette" only. <br>
./run-example.sh vmini <br>
runs the experiment on the "mini-roulette variant" only. <br>
**Note:** Please use the `simplified example names` in the artifact to run the above commands.

# Running on a Custom Input
To run the algorithm on a custom program, please follow these steps: <br>
1.**Create an input file**  whose filename suffix is ".program". This file contains a program, its preconditions, loop gurads and post conditions, and the distributions used for its sampling variables. See the example programs whose filename suffix are ".program" in the “Inputs” folder. Basically, the file should begin with the distributions. This is followed by one # sign. Next, put the preconditions and loop guards in one square bracket. That is followed by one # sign. Then, put the preconditions and post conditions in the bracket followed by one # sign. The main program starts after the third #. We assume this input file is called “input.program” and put it in the "Custom/Custom_inputs" folder.
2.**Create a configuration file** named “inputconfig.txt” in the “Custom/Custom_inputs” folder. The configuration file contains the initial valuations for the custom program. See the example configuration files (whose suffix are "config.txt") in the “Inputs” folder.
3.Run the Script by opening a terminal in the “Tool” folder and running this command:<br>
./run-example.sh input <br>
4.The final output will be put in the "Custom/Custom_outputs" folder automatically.
**Note:** see "readme.txt" in the "Custom" folder, we also put one sample program and its output respectively in the "Custom/Custom_inputs" folder and "Custom/Custom_outputs" folder. (The jpg file describes this sample program.)




