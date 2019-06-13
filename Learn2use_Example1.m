% =====                                                              ====== 
%
%           This file contains the first tutorial of OptiBLESS      
% 
%   It also contains detailed comments explaning how the input file works. 
%   I know this may be a lot to read at first, but you only have to do it
%   once. Input files remain quite similar regardless of the problem you
%   are trying to solve. Please do not hesistate to tell me if you 
%   have difficulties following or suggestions on how to enhance these 
%   tutorials. You can contact me at: 
%          terence.macquart@bristol.ac.uk 
%     (or) terence.macquart@gmail.com
%  
%  In this Example:
%  A single laminate with stacking sequence [45/-45/90/0/45/90/0/45] is
%  used. The lamination parameters corresponding to the given stacking 
%  sequence are pre-calculated (i.e. Lp2Match). The optimiser tries to 
%  match the values of 'Lp2Match' by evaluating various stacking sequence
%  and comparing their lamination parameters against 'Lp2Match'. 
%  A minimisation approach is used (i.e. lower fitness = better match). 
%
% =====                                                              ====== 

% clear previous data and command window
clear all; clc; format short g; format compact; close all;

% Add acces to internal functions 
addpath ([cd '/FitnessFcts'])   % contain the fitness functions
addpath ([cd '/GUI'])           % contains the output interface
addpath ([cd '/src'])           % contains the core source code


%% === Objective 
%  The problem to optimise is setup in this section with the help of the 
%  Objectives structure.
% 
% The objective structure contains 4 fields:
% 
% 1 - Objectives.UserFct is a boolean value (0,1 or false,true)
%       *(false) if set to false, one the default fitness function includes in the
%      'FitnessFcts' folder will be used. The name of the function used is
%       specified by Objectives.FitnessFct
% 
%       *(true) if set to true, the user can specify its own fitness
%       function. Example of user-defined fitness function are also
%       included in the 'FitnessFcts' folder.
% 
% 2 - Objectives.FitnessFct contains the fitness function handle.
% 
% 3 - Objectives.Type can be set to 'LP', 'ABD', or 'SS'
%     The Type string chosen reflects which values will be used by the 
%     fitness function.
%     *(ABD) The stiffness matrices are used as input to the fitness evaluation
%     function and Objectives.FitnessFct must be an handle with @(A,B,D).
% 
%     *(LP) Lamination parameters are used as input to the fitness evaluation
%     function and Objectives.FitnessFct must be an handle with @(LP).
% 
%     *(SS) Stacking Sequence are used as input to the fitness evaluation
%     function and Objectives.FitnessFct must be an handle with @(SS). 
% 
% 4 - Objectives.Table is a cell array that regroups and summarises the
%     optimisation problem. It contains:
%     * The laminate patch number (1st colums)
%     * The Upper and Lower number of plies allowed for this laminate (2nd colums)
%     * The lamination parameter to match (if Type = 'LP') (3rd colums)
%     * Scaling Coefficients prioritise lamination parameter matching (4th column)
%       during the fitness evaluation function.
%       
%     


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
%  the various composite design guidelines are set in this section with
%  the help of the Constraints structures
%  Please see the user manual for details about each design guidelines.
% 
%                        [Symmetry,  Balanced,  Damtol,   Rule10percent,  Disorientation,  Contiguity,  InternalContinuity,  Covering];
Constraints.Vector     = [false   ,    false ,  false ,      false     ,      false     ,     false  ,      false         ,     false];
Constraints.DeltaAngle = 45;    % Available Ply laminate angle are given as [-90:DeltaAngle:90]                      

Constraints.NContiguity   = 3;  % optional (only needed if Contiguity = true)
Constraints.NInternalCont = 3;  % optional (only needed if InternalContinuity = true)

% Employs a penalty based approach to enforce the 10% rule 
% (set the Constraints.Vector (Rule10percent) to false if you set this one to true)
Constraints.Implicit10PercentRule = true; 

%% === Options 
% The optimisation option for the genetic algorithms are set in this
% section.
% 
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
% The main code OptiBLESS is called and the optimisation is run.
% The call returns the Output structure
[Output]  = OptiBLESS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)


% The Output structure contains various fields of interest, in particular:
%  - Output.Table is a cell array summarising the optimisation results
%  - Output.FEASIBLE is a boolean value indicating if the final solution is
%    feasible (true = 1) according to the active design guidelines
%  - Output.fval is the best fitness value
%  - Output.SS_Patch contains the optimised stacking sequence



%% === Checking output results are correctly matching (should not return an error!!!)
%  if this part return an error please contact me as somewthing is
%  definitly wrong.

ScalingCoef = reshape(cell2mat(Objectives.Table(2:end,4)),12,size(Objectives.Table,1)-1);
for i = 2:size(Objectives.Table,1)
    LP2Match = Objectives.Table{i,3};
    if sum(abs(LP2Match-Output.Table{i,4}))>1e-10
        error('non matching LP2Match')
    end
    
    if abs( MYrms ( (Output.Table{i,5}-LP2Match).*ScalingCoef(:,i-1) )-Output.Table{i,7})>1e-10
        error('non matching SS and LPOpt')
    end


    LP = Convert_SS2LP(Output.Table{i,3});
    if sum(abs(LP-Output.Table{i,5}))>1e-10
        error('non matching SS and LPOpt')
    end
    
    if abs( MYrms ( (LP-LP2Match).*ScalingCoef(:,i-1) )-Output.Table{i,7})>1e-10
        error('non matching RSM')
    end

    

    if abs( norm ((LP-LP2Match).*ScalingCoef(:,i-1))-Output.Table{i,6})>1e-10
        error('non matching norm')
    end
end