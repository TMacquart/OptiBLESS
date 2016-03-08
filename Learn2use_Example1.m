% =====                                                              ====== 
%
%       Simple Single Patch Example with Lamination Parameter Matching
%
% =====                                                              ====== 

% Reset and initialise
clear all; clc; format short g; format compact; close all;

% Add acces to all functions 
addpath ./FitnessFcts
addpath ./GUI
addpath ./src



%% === Objective

% --- Staking Sequence used as example,         - Bottom ply [45/-45/90/0/45/90/0/45] Top ply
% --- Laminate paramaters are obtained using:   - Convert_SS2LP([ 45   -45    90     0    45    90     0    45])
Lp2Match = [
            0   % V1A
         0.25   % V2A
            0   % V3A
            0   % V4A
        0.125   % V1B
       0.1875   % V2B
         0.25   % V3B
            0   % V4B
     0.046875   % V1D
       0.4375   % V2D
     -0.46875   % V3D
           0];  % V4D

       
Objectives.UserFct = false;        % set to false to used one of the pre-defined fitness functions
Objectives.Type    = 'LP';         % Set the objective function to lamination parmeter based

ScalingCoef  = ones(12,1);  % Set the relative importance to all lamination parameters (equal in this example)

% --- construct the objective Table
%                        Patch ID            Lower and upper 
%                                       number of plies boundary
Objectives.Table   = [{'Laminate #'}          {'Nplies'}            {'LP2Match'}     {'Scaling Coefficient'} ;
                            {1}                {[8 8]}             {Lp2Match(:,1)}         {ScalingCoef}    ];

Objectives.FitnessFct = @(LP) RMSE_LP(LP,Objectives);       % Attribute the function handle used to calculate fitness
                    


%% === Design Guidelines 
%                        [Symmetry,  Balanced,  Damtol,   Rule10percent,  Disorientation,  Contiguity,  InternalContinuity,  Covering];
Constraints.Vector     = [false   ,    false ,  false ,      false     ,      false     ,     false  ,      false         ,     false];
Constraints.DeltaAngle = 45;                       

Constraints.NContiguity   = 3;  % optional (only needed if Contiguity = true)
Constraints.NInternalCont = 3;  % optional (only needed if InternalContinuity = true)



%% === Options 
GAoptions.Npop    = 20;                     % Population size
GAoptions.Ngen    = 100;                    % Number of generations
GAoptions.NgenMin = 100;                    % Minimum number of generation calculated
GAoptions.Elitism = 0.01;                   % Percentage of elite passing to the next Gen. (from 0 to 1)
GAoptions.PC      = 0.75;                   % Percentage of crossover (from 0.1 to 1)
GAoptions.IniPopFEASIBLE = 1;               % Either (1 or 2), Ensure the initial population 1:respect all design guidelines, 2:and addition respect user function constraints

GAoptions.FitnessLimit = 1e-5;              % Value at which the GA will stop if a fitness is found below this threshold
GAoptions.PlotInterval = [10];              % Refresh plot every X itterations         
GAoptions.SaveInterval = [2];               % Save Data every    X itterations (in Results.txt)   
GAoptions.PlotFct      = @gaplotbestf;      % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;   % Custom ouput function (can be changed to output any information regarding the evolution of the population)



%% === Run
[Output]  = OptiBLESS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)


%% === Checking output results are correctly matching (should not return an error!!!)
ScalingCoef = reshape(cell2mat(Objectives.Table(2:end,4)),12,size(Objectives.Table,1)-1);
for i = 2:size(Objectives.Table,1)
    LP2Match = Objectives.Table{i,3};
    if sum(abs(LP2Match-Output.Table{i,4}))>1e-10
        error('non matching LP2Match')
    end
    
    if abs( rms ( (Output.Table{i,5}-LP2Match).*ScalingCoef(:,i-1) )-Output.Table{i,7})>1e-10
        error('non matching SS and LPOpt')
    end
    
    LP = Convert_SS2LP(Output.Table{i,3});
    if sum(abs(LP-Output.Table{i,5}))>1e-10
        error('non matching SS and LPOpt')
    end
    
    if abs( rms ( (LP-LP2Match).*ScalingCoef(:,i-1) )-Output.Table{i,7})>1e-10
        error('non matching RSM')
    end
    
    if abs( norm ((LP-LP2Match).*ScalingCoef(:,i-1))-Output.Table{i,6})>1e-10
        error('non matching norm')
    end
end