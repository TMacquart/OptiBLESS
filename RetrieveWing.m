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
%
%
% ===                                                                   === 
%         Retrive Stacking Sequence from an optimised wing based on 
%                       lamination parameters                       
% ===                                                                   === 


clear all; close all force; clc; format short g; format compact;

addpath ./FitnessFcts
addpath ./src
addpath ./GUI
addpath ./src/StiffnessOpt


%% Objective Structure
Objectives.Type   = 'LP';
ScalingCoef       = [1 1 1 1, 0 0 0 0, 0 0 0 0]';
Objectives.Table  = [{'Laminate #'}  {'Nplies [LB UB]'}     {'LP2Match'}     {'Scaling Coefficient'} ];

% load ('HighFid_Opt_14LC')
% constant = AnalysisInputs.constant;

load ('BlendedWing')
dvFull_s = dvFull;

TopId = cell2mat(constant.lam.TopID);
TopId = TopId(:);
Nlam  = length(TopId);

Lp2Match = zeros(12,Nlam);

tply = constant.mat.tply * 1;
for j=1:Nlam
    Lp2Match([1 2 3 4  9 10 11 12],j) = dvFull_s ((j-1)*9 + [1:8],end);
    thickness(j) = dvFull_s(j*9,end);
    Nply(j)      = round(thickness(j)/tply);  
    
    if Nply(j)<25,
        Nply(j)=25;
    end
    Objectives.Table = [Objectives.Table; {j} {round([Nply(j)*1 Nply(j)*1])}  Lp2Match(:,j) {ScalingCoef}];  
end


Objectives.UserFct    = false;                        
Objectives.FitnessFct = @(LP) RMSE_LP(LP,Objectives);

%% === Design Guidelines 
%                        [Symmetry,  Balanced,  Damtol,   Rule10percent,  Disorientation,  Contiguity,  InternalContinuity,  Covering];
% Constraints.Vector     = [true   ,    false ,  false ,      false     ,      false     ,     false  ,      false         ,     false];
Constraints.Vector     =   [true   ,    false ,  true ,      false     ,      true     ,     true  ,      true         ,     true];
Constraints.DeltaAngle = 15;       

Constraints.NContiguity   = 3;  % optional (only needed if Contiguity = true)
Constraints.NInternalCont = 3;  % optional (only needed if InternalContinuity = true)


%% === Options 
GAoptions.Npop    = 200;                     % Population size
GAoptions.Ngen    = 10000;                    % Number of generations
GAoptions.NgenMin = 5000;                    % Minimum number of generation calculated
GAoptions.Elitism = 0.01;                   % Percentage of elite passing to the next Gen. (from 0 to 1)
GAoptions.PC      = 0.75;                   % Percentage of crossover (from 0.1 to 1)
GAoptions.IniPopFEASIBLE = 1;               % Either (1 or 2), Ensure the initial population 1:respect all design guidelines, 2:and addition respect user function constraints

GAoptions.FitnessLimit = 1e-5;              % Value at which the GA will stop if a fitness is found below this threshold
GAoptions.PlotInterval = [];              % Refresh plot every X itterations         
GAoptions.SaveInterval = [1];               % Save Data every    X itterations (in Results.txt)   
GAoptions.PlotFct      = @gaplotbestf;      % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;   % Custom ouput function (can be changed to output any information regarding the evolution of the population)



%% == Run
[Output]  = OptiBLESS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)


%% Plot
Coord = constant.lam.coord;

for j=1:Nlam
     PatchXYZ{j}.X = Coord{j}(:,1);
     PatchXYZ{j}.Y = Coord{j}(:,2);
     PatchXYZ{j}.Z = Coord{j}(:,3);
end

plotSS(Output,1,PatchXYZ)


%% write text file output with LP
fr

rms(Output.LP(1:12,1)-Lp2Match(1:12,1))
% LPText = Lp2Match;
LPText = Output.LP;
LPText([5 6 7 8],:) =[];
fid = fopen('LPResults.txt', 'wt'); % Open for writing
for i=1:size(LPText,1) % rows
        fprintf(fid, '%1.0d ', i );
    for j=1:size(LPText,2) % col
        fprintf(fid, '%1.4f ', LPText(i,j) );
    end
    fprintf(fid, '\n' );
end
fclose(fid);


