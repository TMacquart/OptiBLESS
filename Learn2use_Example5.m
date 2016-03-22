% =====                                                              ====== 
%
%                            Generic Example 
%
% This file has been created to help you generate the objective structure
% for any stacking sequence. You can thouroughly test and get a grip on the
% toolbox by yourself this way.
% =====                                                              ====== 

clear all; clc; format short g; format compact; close all Force; fclose all;

addpath ./src
addpath ./FitnessFcts


% Enter any Guide laminates and drops to test the toolbox 
GuideLamDv  = [45 15 -5 0 -45 -60 45 90 -75 80 50 20 -10 0 90 45 -45 -20 -40 0];   % Guide laminate fiber angles                       
Drops       = [{[2 5]} {[12 18]} {[10]} {[20 19 1 4 6 7]}];                        % [{1st group of droped plies} {2nd group} ... ]
ScalingCoef = [1 1 1 1, 1 1 1 1, 1 1 1 1]';                                        % relative importance given to matching the guide laminate LPs integer [1,N], the higher the integer = the more impact on the fit. fct.



% --- Automated calculation (Do Not Change)
% --- this parts prepares the objective Table for you, you can check the
% --- stacking sequence of each laminate in Saved_SS
NUniqueLam         = length(Drops)+1;
Lp2Match           = zeros(12,NUniqueLam);
Objectives.Table   = [{'Id'} {'Nplies'} {'LP2Match'} {'Scaling'}];

for i = 1:NUniqueLam
    Lam = GuideLamDv; 
    DropsLoc = cell2mat(Drops(1:i-1));
    if i ~= 1,  
        Lam(DropsLoc) = [];   
    end
    
    Saved_SS{i}      = Lam; %#ok<SAGROW>
    Lp2Match(:,i)    = Convert_SS2LP(Lam); 
    Objectives.Table = [Objectives.Table; [{i} {[1 1]*length(Lam)} {Lp2Match(:,i)} {ScalingCoef}]];
end    
% ---

Objectives.Type       = 'LP';
Objectives.FitnessFct = @(LP) RMSE_LP(LP,Objectives);
% Objectives.FitnessFct = @(LP) RMSE_MaxAE_LP(LP,Objectives);    % you may want to try other fitness function
Objectives.UserFct = false;

%% === Design Guidelines 
%                        [Symmetry,  Balanced,  Damtol,   Rule10percent,  Disorientation,  Contiguity,  InternalContinuity,  Covering];
Constraints.Vector     = [false   ,    false ,  false ,      false     ,      false     ,     false  ,      false         ,     false];
Constraints.DeltaAngle = 5;                       
Constraints.ply_t      = 0.000127;          % ply thickness    

Constraints.NContiguity   = 3;  % optional (only needed if Contiguity = true)
Constraints.NInternalCont = 3;  % optional (only needed if InternalContinuity = true)
Constraints.Implicit10PercentRule = false; 


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


%% Plot
plotSS(Output,1) 
