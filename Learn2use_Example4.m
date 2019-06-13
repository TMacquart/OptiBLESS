% =====                                                              ====== 
%
%  Simple Three Symmetric Patch Example with Lamination Parameter Matching
%
% =====                                                              ====== 


clear all; clc; format short g; format compact; close all;

addpath ./src
addpath ./FitnessFcts


% GuideLamDv = [-45 0 45 90 0 -45  90 45 0 -45 0 45]s;    
% Lam1       = [-45   45 90 0 -45  90 45 0 -45   45]s;   
% Lam2       = [-45   45 90 0      90 45 0 -45     ]s;   

Lp2Match = [
    % Guide     % Lam 1     % Lam 2
        0.16667    0          0             % V1A
        0          0          0             % V2A
        0         -0.2        0             % V3A
        0          0          0             % V4A
        0          0          0             % V1B
        0          0          0             % V2B
        0          0          0             % V3B
        0          0          0             % V4B
        0.13657   -0.084     -0.11719       % V1D
       -0.12153   -0.114     -0.046875      % V2D
       -0.013889  -0.248     -0.23438       % V3D      
        0          0          0];           % V4D
    

Objectives.Type   = 'LP';
ScalingCoef       = [1 1 1 1, 1 1 1 1, 1 1 1 1]';
Objectives.Table  = [{'Laminate #'}  {'Nplies [LB UB]'}     {'LP2Match'}     {'Scaling Coefficient'} ;
                            {1}           {2*[12 12]}         Lp2Match(:,1)     {ScalingCoef} ;
                            {2}           {2*[10 10]}         Lp2Match(:,2)     {ScalingCoef} ; 
                            {3}           {2*[8 8]}           Lp2Match(:,3)     {ScalingCoef} ; ];

Objectives.FitnessFct = @(LP) RMSE_MaxAE_LP(LP,Objectives);
Objectives.UserFct    = false;

           
%% === Design Guidelines 
%                        [Symmetry,  Balanced,  Damtol,   Rule10percent,  Disorientation,  Contiguity,  InternalContinuity,  Covering];
Constraints.Vector     = [false   ,    false ,  false ,      false     ,      false     ,     false  ,      false         ,     false];
Constraints.DeltaAngle = 45;
Constraints.Implicit10PercentRule = false; 


%% === Options 
GAoptions.Npop    = 100;                     % Population size
GAoptions.Ngen    = 500;                    % Number of generations
GAoptions.NgenMin = 500;                    % Minimum number of generation calculated
GAoptions.Elitism = 0.01;                   % Percentage of elite passing to the next Gen. (from 0 to 1)
GAoptions.PC      = 0.75;                   % Percentage of crossover (from 0.1 to 1)
GAoptions.IniPopFEASIBLE = 1;               % Either (1 or 2), Ensure the initial population 1:respect all design guidelines, 2:and addition respect user function constraints

GAoptions.FitnessLimit = 1e-5;              % Value at which the GA will stop if a fitness is found below this threshold
GAoptions.PlotInterval = [10];              % Refresh plot every X itterations         
GAoptions.SaveInterval = [2];               % Save Data every    X itterations (in Results.txt)   
GAoptions.PlotFct      = @gaplotbestf;      % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;   % Custom ouput function (can be changed to output any information regarding the evolution of the population)



%% == Run
[Output]  = OptiBLESS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)


%% Plot
plotSS(Output,1) 
