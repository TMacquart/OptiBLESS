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
                            {2}           {[8 8]}           Lp2Match(:,2)     {ScalingCoef} ; 
                            {3}           {[6 6]}           Lp2Match(:,3)     {ScalingCoef} ; ];

Objectives.UserFct    = false;                        
Objectives.FitnessFct = @(LP) RMSE_LP(LP,Objectives);

% you can also try other fitness functions
% Objectives.FitnessFct = @(LP) MeanNormE_LP(LP,Objectives); 
% Objectives.FitnessFct = @(LP) MaxAE_LP(LP,Objectives);

           
% =========================== Default Options =========================== %

%                        [Damtol  Rule10percent  Disorientation  Contiguity   BalancedIndirect   InernalContinuity  Covering];
Constraints.Vector     = [false       false          false          false             false            false           false];
Constraints.DeltaAngle = 5;
Constraints.ORDERED    = false;                           
Constraints.Balanced   = true; 
Constraints.Sym        = false; 


% ---
GAoptions.Npop    = 100; 	   % Population size
GAoptions.Ngen    = 500; 	   % Number of generations
GAoptions.NgenMin = 500; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.075; 	   % Percentage of elite passing to the next Gen.
GAoptions.PC      = 0.75; 	   % Plot Boolean

GAoptions.PlotInterval = [10];                  % Refresh plot every X itterations         
GAoptions.SaveInterval = [];                    % Save Data every X itterations   
GAoptions.PlotFct      = @gaplotbestf;          % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;


% ---
[Output]  = RetrieveSS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)


%% Plot
plotSS(Output)


%% Ouput Check (Optional Validation)
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
