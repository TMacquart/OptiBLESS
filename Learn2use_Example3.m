% =====                                                              ====== 
%
%       Simple Three Patch Example with Lamination Parameter Matching
%
% =====                                                              ====== 

clear all; clc; format short g; format compact; close all;

addpath ./src
addpath ./FitnessFcts
addpath ./GUI

% GuideLamDv = [-45 0 45 90 0  -45  45  90  -45  45];    
% Lam1       = [    0    90 0  -45  45  90  -45  45];   
% Lam2       = [    0    90 0           90  -45  45];   

Lp2Match = [
    % Guide     % Lam 1     % Lam 2
        0          0          0             % V1A
        0          0          0             % V2A
       -0.2        0          0.33333       % V3A
        0          0          0             % V4A
       -0.2       -0.25      -0.22222       % V1B
        0.16       0.125      0.11111       % V2B
       -0.24      -0.75      -0.88889       % V3B
        0          0          0             % V4B
        0.048      0.14062    0.22222       % V1D
       -0.048      0.14062    0.22222       % V2D
       -0.488      0.09375    0.037037      % V3D
        0          0          0];           % V4D


Objectives.Type   = 'LP';
ScalingCoef       = [1 1 1 1, 1 1 1 1, 1 1 1 1]';
Objectives.Table  = [{'Laminate #'}  {'Nplies [LB UB]'}     {'LP2Match'}     {'Scaling Coefficient'} ;
                            {1}           {[10 10]}         Lp2Match(:,1)     {ScalingCoef} ;
                            {2}           {[8 8]}           Lp2Match(:,2)     {ScalingCoef*0.7} ; 
                            {3}           {[6 6]}           Lp2Match(:,3)     {ScalingCoef*0.7} ; ];

Objectives.UserFct    = false;                        
Objectives.FitnessFct = @(LP) RMSE_LP(LP,Objectives);

% you can also try other fitness functions
% Objectives.FitnessFct = @(LP) MeanNormE_LP(LP,Objectives); 
% Objectives.FitnessFct = @(LP) MaxAE_LP(LP,Objectives);

           

%% === Design Guidelines 
%                        [Symmetry,  Balanced,  Damtol,   Rule10percent,  Disorientation,  Contiguity,  InternalContinuity,  Covering];
Constraints.Vector     = [false   ,    false ,  false ,      false     ,      false     ,     false  ,      false         ,     false];
Constraints.DeltaAngle = 45;


%% === Options 
GAoptions.Npop    = 50;                     % Population size
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



%% == Run
[Output]  = OptiBLESS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)



%% Plot
plotSS(Output,1) % plot(Output Structure, Window Sacling Factor)



%% Ouput Check (Optional Validation, should not return an error!!)
ScalingCoef = reshape(cell2mat(Objectives.Table(2:end,4)),12,size(Objectives.Table,1)-1);

for i = 2:size(Objectives.Table,1)
    LP2Match = Objectives.Table{i,3};
    if sum(abs(LP2Match-Output.Table{i,4}))>1e-10
        error('non matching LP2Match')
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
