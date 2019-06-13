% =====                                                              ====== 
%
%       Simple Single Patch Example with Stiffness Matrix Matching
%       This is the same example as Learn2useExample1.m except that
%       sitffness matrices are used instead of lamination parameters
% 
% %  In this Example:
%  A single laminate with stacking sequence [45/-45/90/0/45/90/0/45] is
%  used. The stiffness matrices corresponding to the given stacking 
%  sequence are pre-calculated (i.e. A2Match). The optimiser tries to 
%  match the values of 'A2Match' ,'B2Match' and 'D2Match' by evaluating 
%  various stacking sequences and comparing their stifness values against
%  'A2Match' ,'B2Match' and 'D2Match'. A minimisation approach is used 
%  (i.e. lower fitness = better match). 
%
% =====                                                              ====== 

% ----------------------------------------------------------------------- %

clear all; clc; format short g; format compact; close all;

addpath ./FitnessFcts
addpath ./src
addpath ./src/StiffnessOpt      % ----- Supplementary path for stiffness optimisation ------
addpath ./GUI


%% Objective
% In comparison to the first example, material properties are required
% in order to calculate the laminate stiffness matrices. These values
% are stored in the objective structure:
% Objectives.mat = [E1 E2 G12 v12 h] 
% 
% The Scaling coefficient are also replaced by scaling matrices for each of
% the stiffness matrices (IndexAStiff, IndexBStiff and IndexDStiff)
% 
% --- bottom [ 45   -45    90     0    45    90     0    45] Top   
E1   = 13.0e9;  % principal axis young's modulus
E2   = 72.0e9;  % Second axis young's modulus
G12  = 26.9e9;  % Torsional rigidity
v12  = 0.33;    % Poisson's ratio

tply = 0.000127;  % a ply thickness
h    = 8*tply;    % Laminate thickness ( eight time because we assume to know that the laminate contains 8 plies)

% Stiffness matrices to retrieve using OptiBLESS
A2Match ={[
   1.0874e+11   5.8225e+10  -9.2917e+09
   5.8225e+10   1.0874e+11  -9.2917e+09
  -9.2917e+09  -9.2917e+09   2.5255e+10]}; % in-plane stiffness matrix

B2Match ={[
  -9.7029e+09   4.1122e+08  -6.9687e+09
   4.1122e+08   8.8804e+09  -6.9687e+09
  -6.9687e+09  -6.9687e+09   4.1122e+08]}; % Coupled stiffness matrix


D2Match ={[
   1.0602e+11   5.7454e+10   -1.626e+10
   5.7454e+10   1.1299e+11   -1.626e+10
   -1.626e+10   -1.626e+10   2.4484e+10]}; % Out-of-plane stiffness matrix


Objectives.mat = [E1 E2 G12 v12 h];
 

IndexAStiff = ones(3,3); % scaling coefficient for A (note that you may want to use only 6 ones to avoid repeating symmetrical value here (i.e. [1 1 1;0 1 1;0 0 1])
IndexBStiff = ones(3,3); % scaling coefficient for B
IndexDStiff = ones(3,3); % scaling coefficient for D


Objectives.Table   = [{'Laminate #'}     {'Nplies'}   {'A2Match'}  {'B2Match'} {'D2Match'}  {'A Scaling'} {'B Scaling'} {'D Scaling'} ;
                            {1}          {[8 8]}       A2Match       B2Match     D2Match   {IndexAStiff} {IndexBStiff} {IndexDStiff}];

                        
Objectives.Type       = 'ABD';
Objectives.UserFct    = false;
Objectives.FitnessFct = @(A,B,D) RMSE_ABD(A,B,D,Objectives);


%% === Design Guidelines 
%                        [Symmetry,  Balanced,  Damtol,   Rule10percent,  Disorientation,  Contiguity,  InternalContinuity,  Covering];
Constraints.Vector     = [false   ,    false ,  false ,      false     ,      false     ,     false  ,      false         ,     false];
Constraints.DeltaAngle = 45;                       
Constraints.ply_t      = 0.000127;          % ply thickness    

Constraints.NContiguity   = 3;  % optional (only needed if Contiguity = true)
Constraints.NInternalCont = 3;  % optional (only needed if InternalContinuity = true)
Constraints.Implicit10PercentRule = false; 


%% === Options 
GAoptions.Npop    = 100;                    % Population size
GAoptions.Ngen    = 200;                    % Number of generations
GAoptions.NgenMin = 200;                    % Minimum number of generation calculated
GAoptions.Elitism = 0.01;                   % Percentage of elite passing to the next Gen. (from 0 to 1)
GAoptions.PC      = 0.75;                   % Percentage of crossover (from 0.1 to 1)
GAoptions.IniPopFEASIBLE = 1;               % Either (1 or 2), Ensure the initial population 1:respect all design guidelines, 2:and addition respect user function constraints


GAoptions.FitnessLimit = 1e5;               % ----- Note the value here is much higher due the high value of stiffness matrices ------


GAoptions.PlotInterval = [10];              % Refresh plot every X itterations         
GAoptions.SaveInterval = [2];               % Save Data every    X itterations (in Results.txt)   
GAoptions.PlotFct      = @gaplotbestf;      % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;   % Custom ouput function (can be changed to output any information regarding the evolution of the population)


%% === Run
[Output]  = OptiBLESS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)

% Note that the fitness value is not zero because of truncation errors of
% the input stiffness matrices

%% --- Compare the retrieved stiffness matrices with the one given as input
[AOpt,BOpt,DOpt] = Convert_SS2ABD(E1,E2,v12,G12,tply,cell2mat(Output.SS_Patch),true);                  

fprintf('\n')
display('--- A Matrix ---')
display(A2Match{1})     % input matrix
display(AOpt)           % Retrieved Matrix
fprintf('\n')
display('--- B Matrix ---')
display(B2Match{1})
display(BOpt)
fprintf('\n')
display('--- D Matrix ---')
display(D2Match{1})
display(DOpt)
