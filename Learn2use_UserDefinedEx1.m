% =====                                                              ====== 
%
%                    User defined fitness function Example 
%
% This file will help you setup your own fitness function.
% Once your run this file, it will take you to an empty template for
% fitness functions
% =====                                                              ====== 

clear all; clc; format short g; format compact; close all;

addpath ./FitnessFcts
addpath ./GUI
addpath ./src


%% === Objective

% --- Corresponding Staking Sequence
% --- Bottom ply [ 45   -45    90     0    45    90     0    45] Top ply
Lp2Match = [
            0 % V1A
         0.25 % V2A
            0 % V3A
            0 % V4A
        0.125 % V1B
       0.1875 % V2B
         0.25 % V3B
            0 % V4B
     0.046875 % V1D
       0.4375 % V2D
     -0.46875 % V3D
           0];% V4D


Objectives.Type    = 'LP'; 


ScalingCoef        = ones(12,1); 


Objectives.Table   = [{'Laminate #'}     {'Nplies'}      {'LP2Match'}     {'Scaling Coefficient'} ;
                            {1}           {[8 8]}         {Lp2Match(:,1)}  {ScalingCoef} ; ];
                    
                        
% --- CHANGE HERE --- CHANGE HERE --- CHANGE HERE                           
% --- those 2 lines are the only line that have been changed compared to Learn2use_Example1.m   
% --- you may also want to have a look at the UserDefined_RMSE_LP.m file 
Objectives.UserFct    = true;    % this specify that you are using your own function for fitness calculations                      
% Objectives.FitnessFct = @(LP) UserDefined_RMSE_LP(LP,Objectives);
Objectives.FitnessFct = @(LP) UserDefined_LPFitness(LP); % for a more general function try (you will need to complete the fitness calculation)



%% === Design Guidelines 
%                        [Symmetry,  Balanced,  Damtol,   Rule10percent,  Disorientation,  Contiguity,  InternalContinuity,  Covering];
Constraints.Vector     = [false   ,    false ,  false ,      false     ,      false     ,     false  ,      false         ,     false];
Constraints.DeltaAngle = 5;                       
Constraints.ply_t      = 0.000127;          % ply thickness    


Constraints.NContiguity   = 3;  % optional (only needed if Contiguity = true)
Constraints.NInternalCont = 3;  % optional (only needed if InternalContinuity = true)
Constraints.Implicit10PercentRule = false; % the 10% rule can be implicitely or explicitely(Constraints.Vector) handled

%% === Options 
GAoptions.Npop    = 100;                    % Population size
GAoptions.Ngen    = 2000;                    % Number of generations
GAoptions.NgenMin = 1000;                    % Minimum number of generation calculated
GAoptions.Elitism = 0.01;                   % Percentage of elite passing to the next Gen. (from 0 to 1)
GAoptions.PC      = 0.75;                   % Percentage of crossover (from 0.1 to 1)
GAoptions.IniPopFEASIBLE = 1;               % Either (1 or 2), Ensure the initial population 1:respect all design guidelines, 2:and addition respect user function constraints


GAoptions.FitnessLimit = 1e-5;               % ----- Note the value here is much higher due the high value of stiffness matrices ------


GAoptions.PlotInterval = [10];              % Refresh plot every X itterations         
GAoptions.SaveInterval = [2];               % Save Data every    X itterations (in Results.txt)   
GAoptions.PlotFct      = @gaplotbestf;      % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;   % Custom ouput function (can be changed to output any information regarding the evolution of the population)


%% === Run
[Output]  = OptiBLESS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)
