% =====                                                              ==== 
%  This file contains the typical default inputs that can be set. 
%  Copy and paste the part that you need into another file to help you start. 
%  
%         Note that this File will not successfully run by default.
% =====                                                              ==== 

clear all; clc; format short g; format compact; close all;

addpath ./FitnessFcts
addpath ./GUI
addpath ./src
addpath ./src/StiffnessOpt



%% === Objective
% The default fitness function pre-coded in SSORT require the objective structure. 
% The objective structure mostly contains data about fitness calculation.
% There are 3 different objective structures depending on the type of
% fitness functions that you will be using. 
%  It can either be a fitness function based on:
%           - Lamination parameters ('LP')
%           - Stiffness Matrices ('ABD')
%       or  - Stacking Sequences ('SS') 


% =================== Objectives.Table for Lamination Parameters (LP)

% Laminate #           : Laminate ID defining a particular patch (should be unique for each row) 
% Nplies [LB, UB]      : Lower and Upper number of plies allowed for this patch 
% LP2Match             : Lamination parameters that are desired for this patch
% Scaling Coefficients : Coefficient use to quantify the relative matching importance 
%                        of each lamination parameter (if no preference, use a vector of ones)

Objectives.Table   = [{'Laminate #'}      {'Nplies [LB,UB]'}        {'LP2Match'}         {'Scaling Coefficients'} ;
%                      {Integer}        {[Integer  Integer]}     {(12x1 LP Vector)}     {12x1 Real Coefficient Vector};  
                        { }                  { }                           { }                      { }           ]; % you may have as many line as desired


                    
                        
% =================== Objectives.Table for Stiffness Matrices (ABD)
% A similar one is used for Stiffness matching optimisation

Objectives.Table   = [{'Laminate #'}     {'Nplies'}          {'A2Match'}     {'B2Match'}     {'D2Match'}     {'A Scaling'}  {'B Scaling'}   {'D Scaling'} ;
%                       {Integer}    {[Integer  Integer]}  {[3x3 Matrix]}   {[3x3 Matrix]}  {[3x3 Matrix]}   {[3x3 Matrix]} {[3x3 Matrix]}  {[3x3 Matrix]}; 
                           { }                { }               { }              { }            { }              { }            { }             { }];

% A2Match   : In-plane stiffness matrice parameters that are desired for this patch                         
% A Scaling : Coefficient use to quantify the relative importance of each of the in-plane
%             stiffness matrice parameters that are desired for this patch   




% =================== Objectives.Table for Stacking Sequence (SS) and
% =================== minimum required for User based fitness functions

Objectives.Table   = [{'Laminate #'}      {'Nplies [LB,UB]'}    ;
%                      {Integer}          {[Integer  Integer]}  ;  
                         { }                      { }            ];
                     
                   
Objectives.Type = 'LP';          % The objective Type defines what kind of fitness function you are using.  Either ('LP'), ('ABD') or ('SS')      

Objectives.UserFct = false;      % set to true if your are not using one of the default fitness function provided with the toolbox 

Objectives.FitnessFct = @(LP) RMSE_LP(LP,Objectives);    % Enter your fitness function (stored in the FitnessFcts folder)




%% === Constraints 

% Constraints.Vector is a Vector of Boolean use to activate / deactivate
% composite design guidelines.

%                       [Damtol  Rule10percent  Disorientation  Contiguity   BalancedIndirect  InernalContinuity  Covering];
Constraints.Vector     = [false       false          false          false         false            false            false];

% {Damtol}.             Damage Tolerance,  +-45 plies are used for the upper and lower laminate plies.
% {Rule10percent}.      A minimum of 10% of plies in each of the 0, ±45 and 90 is enforced.
% {Disorientation}.     The change of angles between two consecutive plies should not exceed 45
% {Contiguity}.         The change of angles between two consecutive plies should not be below 5
% {DiscreteAngle}.      Discrete fibre angles are used (possible values are set with Constraints.DeltaAngle (below)
% {InernalContinuity}.  One ply must be kept spanning the entire structure every three plies.
% {Covering}.           Covering plies on the lower and upper surfaces of the laminate cannot be dropped. 


Constraints.DeltaAngle = 45;       % Possible ply angles values are defined as -90:Constraints.DeltaAngle:90

Constraints.ORDERED    = false;    % If laminates are given in descending order in Objective table, this order will be kept if set to true (keep false in most cases)

Constraints.Balanced   = false;    % All fibre angles, including 0 and 90 (current limitation), occur in pairs if set to true

Constraints.Sym        = false;    % Stacking sequence is mirrored about the mid-plane if set to true



%% === Options 
GAoptions.Npop    = 100; 	   % Population size
GAoptions.Ngen    = 100; 	   % Number of generations
GAoptions.NgenMin = 100; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.01; 	   % Percentage of elite passing to the next Gen.
GAoptions.PC      = 0.75; 	   % Percentage of crossover

GAoptions.PlotInterval = [10];                  % Refresh plot every X itterations         
GAoptions.SaveInterval = [2];                   % Save Data every X itterations   
GAoptions.PlotFct      = @gaplotbestf;          % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;



%% === Run
[Output]  = RetrieveSS(Objectives,Constraints,GAoptions); % This is the function you need to call to start the optimisation (always in this format)

display(Output)             % Display general output structure
display(Output.Table)       % Display summary of the results



%% === Plot
plotSS(Output)              % Call the GUI
