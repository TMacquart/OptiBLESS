% =====                                                              ==== 
%  This file contains the typical default inputs that can be set. 
%  Copy and paste the part that you need into another file to help you start. 
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
% fitness function that you will be using. 
%  It can either be a fitness function based on:
%           - Lamination parameters ('LP')
%           - Stiffness Matrices ('ABD')
%       or  - Stacking Sequences ('SS') 


% =================== Objectives.Table for Lamination Parameters (LP)

% --- The Table regroups the most important parameters

% Laminate #           : Laminate ID defining a particular patch (should be unique for each row) 
% Nplies [LB, UB]      : Lower and Upper number of plies allowed for this patch 
% LP2Match             : Lamination parameters that are desired for this patch
% Scaling Coefficients : Coefficient use to quantify the relative matching importance 
%                        of each lamination parameter (if no preference use a vector of ones)

Objectives.Table   = [{'Laminate #'}      {'Nplies [LB,UB]'}        {'LP2Match'}         {'Scaling Coefficients'} ;
%                      {Integer}        {[Integer  Integer]}     {(12x1 LP Vector)}     {12x1 Real Coefficient Vector};  
                        { }                  { }                           { }                      { }           ]; % you may have as many line as desired


                        
% =================== Objectives.Table for Stiffness Matrices (ABD)

Objectives.Table   = [{'Laminate #'}     {'Nplies'}          {'A2Match'}     {'B2Match'}     {'D2Match'}     {'A Scaling'}  {'B Scaling'}   {'D Scaling'} ;
%                       {Integer}    {[Integer  Integer]}  {[3x3 Matrix]}   {[3x3 Matrix]}  {[3x3 Matrix]}   {[3x3 Matrix]} {[3x3 Matrix]}  {[3x3 Matrix]}; 
                           { }                { }               { }              { }            { }              { }            { }             { }];


% =================== Objectives.Table for Stacking Sequence (SS) and
% =================== minimum required for User based fitness functions

Objectives.Table   = [{'Laminate #'}      {'Nplies [LB,UB]'}    ;
%                      {Integer}          {[Integer  Integer]}  ;  
                         { }                      { }            ];
                     
                   
% --- Objectives Type                        
%       The objective Type defines what kind of fitness function you are using. 
Objectives.Type = 'LP';  % Either ('LP'), ('ABD') or ('SS')      


Objectives.UserFct = false;   % set to true if your are not using one of the default fitness function provided with the toolbox 

Objectives.FitnessFct = @(LP) RMSE_LP(LP,Objectives); % Enter your fitness function (stored in the FitnessFcts folder)


%% === Constraints 

% Constraints.Vector is a Vector of Boolean use to activate / deactivate constraints

%                       [Damtol  Rule10percent  Disorientation  Contiguity   BalancedIndirect  InernalContinuity  Covering];
Constraints.Vector     = [false       false          false          false         false            false            false];

% {Symmetry}.           Stacking sequence is mirrored about the mid-plane.
% {Balance}.            All fibre angles, except 0$^{\circ}$ and 90$^{\circ}$, occur in $\pm$ pairs. 
% {Damtol}.             Damage Tolerance,  $\pm$45$^{\circ}$ plies are used for the upper and lower laminate plies.
% {Rule10percent}.      A minimum of 10\% of plies in each of the 0$^{\circ}$, ±45$^{\circ}$ and 90$^{\circ}$ is enforced.
% {Disorientation}.     The change of angles between two consecutive plies should not exceed 45$^{\circ}$.
% {Contiguity}.         The change of angles between two consecutive plies should not be below 5$^{\circ}$.
% {DiscreteAngle}.      Discrete fibre angles are used (possible values are set with  \textit{$\Delta$Angle}).  
% {InernalContinuity}.  One ply must be kept spanning the entire structure every three plies.
% {Covering}.           Covering plies on the lower and upper surfaces of the laminate cannot be dropped. 



Constraints.DeltaAngle = 45;       % Possible ply angles are defined as -90:Constraints.DeltaAngle:90

Constraints.ORDERED    = false;    % If laminate are given in descending order in Objective table, this order will be kept if set to true 

Constraints.Balanced   = false;    % All fibre angles, including 0 and 90 (current limitation), occur in pairs if se to true

Constraints.Sym        = false;    % Stacking sequence is mirrored about the mid-plane. if se to true



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
[Output]  = RetrieveSS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)



%% === Plot
plotSS(Output)
