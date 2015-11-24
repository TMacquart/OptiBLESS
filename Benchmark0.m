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
                        
ply_t      = 0.000127;
GuideLamDv = [+45 -45 90 0];
% GuideLamDv = [-90 -45 -45];
Drops      = []; %[2 4 6];
% GuideLam   = [GuideLamDv, -GuideLamDv];
GuideLam   = [GuideLamDv, -GuideLamDv, fliplr([GuideLamDv, -GuideLamDv])]'; % balanced/symetric

ScalingCoef = [1 0 1 0, 1 1 1 1, 1 1 1 1]';      % relative importance given to matching the guide laminate LPs integer [1,N], the higher the integer = the more impact on the fit. fct.
NUniqueLam = length(Drops)+1;
Lp2Match   = zeros(12,NUniqueLam);
for i = 1:NUniqueLam
    Lam = GuideLamDv;
    if i ~= 1,        
        Lam(Drops(1:i-1)) = [];   
    end
%     Lam = [Lam, -Lam]'; % balanced/symetric 
    Lam = [Lam, -Lam, fliplr([Lam, -Lam])]'; % balanced/symetric 
    
    Lp2Match(:,i) = Convert_SS2LP(Lam);
    Objectives.Table = [Objectives.Table; [{i} {[1 1]*length(Lam)} {Lp2Match(:,i)} {ScalingCoef}]];
end

Objectives.Type        = 'LP';
Objectives.FitnessFct = @(LP) SumRMSLP(LP,Objectives);

                        
% =========================== Default Options =========================== %

%                        [Damtol  Rule10percent  Disorientation  Contiguity   DiscreteAngle  InernalContinuity  Covering  Balanced];
Constraints.Vector     = [false       false          false          false         true            false            false  true];
Constraints.DeltaAngle = 45;
Constraints.ply_t      = ply_t;          % ply thickness
Constraints.Balanced   = false; 
Constraints.Sym        = false; 
Constraints.ORDERED    = false;           


% ---
GAoptions.Npop    = 300; 	   % Population size
GAoptions.Ngen    = 1000; 	   % Number of generations
GAoptions.NgenMin = 500; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.075; 	   % Percentage of elite passing to the next Gen.
GAoptions.Plot    = true; 	   % Plot Boolean
GAoptions.PC      = 0.5;

% ---
Nrun = 1;
output_Match = cell(1,Nrun);
feasible     = zeros(1,Nrun);
fval         = zeros(1,Nrun);
MeanRMS      = zeros(1,Nrun);
Ngens        = zeros(1,Nrun);
NFctEval     = zeros(1,Nrun);
MeanNorm     = zeros(1,Nrun);
for i = 1:Nrun
    [output_Match{i}] = RetrieveSS(Objectives,Constraints,GAoptions);
    feasible(i)       = output_Match{i}.FEASIBLE;
    fval(i)           = output_Match{i}.fval;
    MeanRMS(i)        = mean(cell2mat(output_Match{i}.Table(2:end,end)));
    MeanNorm(i)       = mean(cell2mat(output_Match{i}.Table(2:end,end-1)));
    NFctEval(i)       = output_Match{i}.NfctEval;  
    Ngens(i)          = output_Match{i}.NGen; 

end

if 0
    figure(2)
    hist(fval,20)
    hist(NFctEval,20)
    hist(Ngens)
    
end
fr
% output_Match = output_Match{i}

%% Checking output results are correct
ScalingCoef = reshape(cell2mat(Objectives.Table(2:end,4)),12,size(Objectives.Table,1)-1);

for i = 2:size(Objectives.Table,1)
    LP2Match = Objectives.Table{i,3};
    if sum(abs(LP2Match-output_Match.Table{i,4}))>1e-10
        error('non matching LP2Match')
    end
    
    LP = Convert_SS2LP(output_Match.Table{i,3});
    if sum(abs(LP-output_Match.Table{i,6}))>1e-10
        error('non matching SS and LPOpt')
    end
    
    if abs( rms ( (LP-LP2Match).*ScalingCoef(:,i-1) )-output_Match.Table{i,8})>1e-10
        error('non matching RSM')
    end
    
    if abs( norm ((LP-LP2Match).*ScalingCoef(:,i-1))-output_Match.Table{i,7})>1e-10
        error('non matching norm')
    end
    
    if abs( 100*sum(abs(  ((LP-LP2Match)./LP2Match).*ScalingCoef(:,i-1) )) - output_Match.Table{i,6})>1e-10
        error('non matching error percent')
    end
end