% =====                                                              ==== 
%               This file is an user input file example of SSR
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

clear all; clc; format short g; format compact; close all Force;

global Pop

if 0
    for i=1:length(Pop)
        MeanNply(i) = mean(Pop{i}(:,1));
        StdNply(i) = std(Pop{i}(:,1));
        
        MeanTheta1(i) = mean(Pop{i}(:,2));
        StdTheta1(i) = std(Pop{i}(:,2));
    end
    figure
    hold all
%     plot(MeanTheta1)
%     plot(StdTheta1)
    plot(MeanNply)
    plot(StdNply)
end

addpath ./FitnessFcts

% --- Creating a real SS to match
Objectives.Table   = [{'Laminate Index'} {'Nplies'} {'LP2Match'} {'Importance'}];
                        
ply_t = 0.000127;



% GuideLamDv = [+45 0 -45 0]; % test

%% Balanced Symmetric
% GuideLamDv = [+45 0 -45 90];                                                   % theta 1
% GuideLamDv = [+45 0 -45 90 45 0 -45 0];                                        % theta 2
GuideLamDv = [-45 0 45  90 0  -45  45  90  -45  45];                           % theta 3
% GuideLamDv = [0 -45 45 45 -45   0 45 90  90 0 -45 45 0  -45 45 90 90 -45 0   0]; % theta 4
% GuideLam   = [GuideLamDv, fliplr(GuideLamDv)];

%% Symmetric 
% GuideLamDv = [45   -15    30   -55    40   -80   -40   -80   -70    60];                                                              % theta 5
% GuideLamDv = [40   -30    85   -80   -10   -20    50    55   -55     0  -5    30    40    50   -40    35    30   -60   -65     0];    % theta 6
% GuideLam   = [GuideLamDv, fliplr(GuideLamDv)];

%% Balanced
% GuideLamDv = [0 30 -40 -30 85 -85 40 -70 70 0];                                                                                       % theta 7
% GuideLamDv = [-50 85 -40  -25  20 25  -45 -85 50 -20 -40  5 -5 -75 40  75  45 -85 40 85];                                             % theta 8

%% Generic 
% GuideLamDv = [10   -65   -60   -40    65   -40    60   -45    80   -25];                                                              % theta 9
% GuideLamDv = [-50   -40    25     0   -25    60    20    10    80   -35 50    50   -20    15   -75   -80    10    55    80   -65];    % theta 10

% GuideLamDv = randi([1 36],1,100)*5-90;
% Drops = num2cell(randperm(100,50));

%%
% GuideLamDv = [-90 -45 -45];
Drops      = [{[4 6 7 8]}] %[{[4 5 6 11 15 17 20]}]%[{[2 3 8]}]; %[2 4 6];
% GuideLam   = [GuideLamDv, fliplr(GuideLamDv)];
% GuideLam   = [GuideLamDv, -GuideLamDv];
% GuideLam   = [GuideLamDv, -GuideLamDv, fliplr([GuideLamDv, -GuideLamDv])]'; % balanced/symetric

ScalingCoef = [1 1 1 1, 1 1 1 1, 1 1 1 1]';      % relative importance given to matching the guide laminate LPs integer [1,N], the higher the integer = the more impact on the fit. fct.
NUniqueLam  = length(Drops)+1;
Lp2Match    = zeros(12,NUniqueLam);
for i = 1:NUniqueLam
    Lam = GuideLamDv;
    DropsLoc = cell2mat(Drops(1:i-1));
    if i ~= 1,  
        Lam(DropsLoc) = [];   
    end
    
       Lam   = [Lam, fliplr(Lam)];
       
%     Lam = [Lam, -Lam]'; % balanced/symetric 
%     Lam = [Lam, -Lam, fliplr([Lam, -Lam])]'; % balanced/symetric 
%     keyboard
    Lp2Match(:,i)    = Convert_SS2LP(Lam);
    Objectives.Table = [Objectives.Table; [{i} {[1 1]*length(Lam)} {Lp2Match(:,i)} {ScalingCoef}]];
end

Objectives.Type        = 'LP';
Objectives.FitnessFct = @(LP) RMSE_MaxAE_LP(LP,Objectives);

% =========================== Default Options =========================== %

%                        [Damtol  Rule10percent  Disorientation  Contiguity   BalancedIndirect  InernalContinuity  Covering  ];
Constraints.Vector     = [false       false          false          false         false            false            false       ];
Constraints.DeltaAngle = 45;
Constraints.ply_t      = ply_t;      % ply thickness
Constraints.Balanced   = true;      % Direct Constraint Handling
Constraints.Sym        = true; 
Constraints.ORDERED    = false;           



% ---
GAoptions.Npop    = 200; 	   % Population size
GAoptions.Ngen    = 500; 	   % Number of generations
GAoptions.NgenMin = 500; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.075; 	   % Percentage of elite passing to the next Gen.
GAoptions.Plot    = false; 	   % Plot Boolean
GAoptions.PC      = 0.55;

% ---
Nrun = 1
output_Match = cell(1,Nrun);
feasible     = zeros(1,Nrun);
fval         = zeros(1,Nrun);
MeanRMS      = zeros(1,Nrun);
Ngens        = zeros(1,Nrun);
NFctEval     = zeros(1,Nrun);
MeanNorm     = zeros(1,Nrun);


for i = 1:Nrun
    display(i)
    [output_Match{i}] = RetrieveSS(Objectives,Constraints,GAoptions);
    feasible(i)       = output_Match{i}.FEASIBLE;
    fval(i)           = output_Match{i}.fval;
    
    MeanNormE(i)      = mean(cell2mat(output_Match{i}.Table(2:end,end-3)));
    MeanRMSE(i)       = mean(cell2mat(output_Match{i}.Table(2:end,end-2)));
    MeanMAE(i)        = mean(cell2mat(output_Match{i}.Table(2:end,end-1)));
    MeanMaxAE(i)      = mean(cell2mat(output_Match{i}.Table(2:end,end)));
    
    NFctEval(i)       = output_Match{i}.NfctEval;  
    Ngens(i)          = output_Match{i}.NGen; 
    display(fval(i))
end

fr
if 0
    NormMAE = [[-0.01:0.01:0.1] [0.12 0.14 0.16 0.18 0.2 0.3 0.4 0.5]]
    figure(2)
    hist(MeanMAE,NormMAE)
    [N1,X] = hist(MeanMAE,NormMAE);
    [N2,X] = hist(MeanMaxAE,NormMAE);
    [N1' N2'  X']

    
    [[1:Nrun]' fval' NFctEval' Ngens' MeanNorm']
end

% output_Match = output_Match{1}

%% Checking output results are correct
ScalingCoef = reshape(cell2mat(Objectives.Table(2:end,4)),12,size(Objectives.Table,1)-1);

for i = 2:size(Objectives.Table,1)
    LP2Match = Objectives.Table{i,3};
    if sum(abs(LP2Match-output_Match.Table{i,4}))>1e-10
        error('non matching LP2Match')
    end
    
    LP = Convert_SS2LP(output_Match.Table{i,3});
    if sum(abs(LP-output_Match.Table{i,5}))>1e-10
        error('non matching SS and LPOpt')
    end
    
    if abs( rms ( (LP-LP2Match).*ScalingCoef(:,i-1) )-output_Match.Table{i,7})>1e-10
        error('non matching RSM')
    end
    
    if abs( norm ((LP-LP2Match).*ScalingCoef(:,i-1))-output_Match.Table{i,6})>1e-10
        error('non matching norm')
    end
end