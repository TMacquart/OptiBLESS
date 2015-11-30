% =====                                                              ==== 
%  This file is a typical user input file that is used to run the code.
% =====                                                              ==== 

% ----------------------------------------------------------------------- %
% Copyright (c) <2015>, <Terence Macquart>
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% The views and conclusions contained in the software and documentation are those
% of the authors and should not be interpreted as representing official policies,
% either expressed or implied, of the FreeBSD Project.
% ----------------------------------------------------------------------- %
clear all; clc; format short g; format compact; close all;

addpath ./StiffnessOpt
addpath ./FitnessFcts

% GuideLamDv = [-45 0 45 90 0  -45  45  90  -45  45]s;    
% Lam1       = [    0    90 0  -45  45  90  -45  45]s;   
% Lam2       = [    0    90 0           90  -45  45]s;   

Lp2Match = [
    % Guide     % Lam 1     % Lam 2
        0          0          0             % V1A
        0          0          0             % V2A
       -0.2        0          0.33333       % V3A
        0          0          0             % V4A
        0          0          0             % V1B
        0          0          0             % V2B
        0          0          0             % V3B
        0          0          0             % V4B
        0.162      0.22266    0.22222       % V1D
       -0.132    -0.058594    -0.027778     % V2D
       -0.092      0.58594    0.92593       % V3D  
        0          0          0];           % V4D
    

Objectives.Type   = 'LP';
ScalingCoef       = [1 1 1 1, 1 1 1 1, 1 1 1 1]';
Objectives.Table  = [{'Laminate #'}  {'Nplies [LB UB]'}     {'LP2Match'}     {'Scaling Coefficient'} ;
                            {1}           {2*[10 10]}         Lp2Match(:,1)     {ScalingCoef} ;
                            {2}           {2*[8 8]}           Lp2Match(:,2)     {ScalingCoef} ; 
                            {3}           {2*[6 6]}           Lp2Match(:,3)     {ScalingCoef} ; ];

Objectives.FitnessFct = @(LP) RMSE_MaxAE_LP(LP,Objectives);


           
% =========================== Default Options =========================== %

%                        [Damtol  Rule10percent  Disorientation  Contiguity   BalancedIndirect   InernalContinuity  Covering];
Constraints.Vector     = [false       false          false          false             false            false            false];
Constraints.DeltaAngle = 5;
Constraints.ORDERED    = false;                           
Constraints.Balanced   = true; 
Constraints.Sym        = true; 


% ---
GAoptions.Npop    = 100; 	   % Population size
GAoptions.Ngen    = 500; 	   % Number of generations
GAoptions.NgenMin = 500; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.075; 	   % Percentage of elite passing to the next Gen.
GAoptions.Plot    = true; 	   % Plot Boolean
GAoptions.PC      = 0.75; 	   % Plot Boolean



% ---
[Output]  = RetrieveSS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)
