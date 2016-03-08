% =====                                                              ==== 
%  This file is a typical user input file that is used to run the code.
% =====                                                              ==== 


clear all; clc; format short g; format compact; close all;

addpath ./FitnessFcts
addpath ./src
addpath ./src/StiffnessOpt
addpath ./GUI

% --- bottom [ 45   -45    90     0    45    90     0    45] Top   
E1   = 13.0e9;
E2   = 72.0e9;
G12  = 26.9e9;
v12  = 0.33;

tply = 0.000127;  % ply thickness
h    = 8*tply;

A2Match ={[
   1.0874e+11   5.8225e+10  -9.2917e+09
   5.8225e+10   1.0874e+11  -9.2917e+09
  -9.2917e+09  -9.2917e+09   2.5255e+10]};
B2Match ={[
  -9.7029e+09   4.1122e+08  -6.9687e+09
   4.1122e+08   8.8804e+09  -6.9687e+09
  -6.9687e+09  -6.9687e+09   4.1122e+08]};
D2Match ={[
   1.0602e+11   5.7454e+10   -1.626e+10
   5.7454e+10   1.1299e+11   -1.626e+10
   -1.626e+10   -1.626e+10   2.4484e+10]};

Objectives.mat = [E1 E2 G12 v12 h];
 
IndexAStiff = ones(3,3);
IndexBStiff = ones(3,3);
IndexDStiff = ones(3,3);


Objectives.Table   = [{'Laminate #'}     {'Nplies'}   {'A2Match'}  {'B2Match'} {'D2Match'}  {'A Scaling'} {'B Scaling'} {'D Scaling'} ;
                            {1}          {[8 8]}       A2Match       B2Match     D2Match   {IndexAStiff} {IndexBStiff} {IndexDStiff}];

                        
Objectives.Type       = 'ABD';
Objectives.UserFct    = true;
Objectives.FitnessFct = @(A,B,D) UserDefined_ABDFitness(A,B,D);


%% === Design Guidelines 
%                        [Symmetry,  Balanced,  Damtol,   Rule10percent,  Disorientation,  Contiguity,  InternalContinuity,  Covering];
Constraints.Vector     = [false   ,    false ,  false ,      false     ,      false     ,     false  ,      false         ,     false];
Constraints.DeltaAngle = 5;                       
Constraints.ply_t      = 0.000127;          % ply thickness    

Constraints.NContiguity   = 3;  % optional (only needed if Contiguity = true)
Constraints.NInternalCont = 3;  % optional (only needed if InternalContinuity = true)



%% === Options 
GAoptions.Npop    = 100;                    % Population size
GAoptions.Ngen    = 2000;                    % Number of generations
GAoptions.NgenMin = 1000;                    % Minimum number of generation calculated
GAoptions.Elitism = 0.01;                   % Percentage of elite passing to the next Gen. (from 0 to 1)
GAoptions.PC      = 0.75;                   % Percentage of crossover (from 0.1 to 1)
GAoptions.IniPopFEASIBLE = 1;               % Either (1 or 2), Ensure the initial population 1:respect all design guidelines, 2:and addition respect user function constraints
GAoptions.FitnessLimit = 1e-5;               

GAoptions.PlotInterval = [10];              % Refresh plot every X itterations         
GAoptions.SaveInterval = [2];               % Save Data every    X itterations (in Results.txt)   
GAoptions.PlotFct      = @gaplotbestf;      % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;   % Custom ouput function (can be changed to output any information regarding the evolution of the population)


%% === Run
[Output]  = OptiBLESS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)
