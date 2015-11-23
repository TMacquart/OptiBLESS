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

addpath ./FitnessFcts

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

% --- bottom [ 45   -45    90     0    45    90     0    45] Top
Lp2Match = [
            0 % V1A
         0.25 % V2A
            0 % V3A
            0 % V4A
        0.125 % V1B
       0.1875 % V2B
         0.25 % V3B
            0 % V4B
     0.046875 % V1D
       0.4375 % V2D
     -0.46875 % V3D
           0];% V4D


ScalingCoef = ones(12,1); 
Objectives.Table   = [{'Laminate #'}     {'Nplies'}      {'LP2Match'}     {'Scaling Coefficient'} ;
                            {1}          {[6 10]}         {Lp2Match(:,1)}  {ScalingCoef} ; ];

Objectives.Type       = 'LP'; 
Objectives.FitnessFct = @(LP) SumRMSLP(LP,Objectives);


% =========================== Default Options =========================== %

%                        [Damtol  Rule10percent  Disorientation  Contiguity   DiscreteAngle  InernalContinuity  Covering];
Constraints.Vector     = [false       false          false          false         true            false            false];
Constraints.DeltaAngle = 45;
Constraints.ply_t      = 0.000127;          % ply thickness
Constraints.ORDERED    = true;                         
Constraints.Balanced   = true; 
Constraints.Sym        = true; 



% ---
GAoptions.Npop    = 100; 	   % Population size
GAoptions.Ngen    = 500; 	   % Number of generations
GAoptions.NgenMin = 250; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.01; 	   % Percentage of elite passing to the next Gen.
GAoptions.PC      = 0.5; 	   % Percentage of crossover
GAoptions.Plot    = true; 	   % Plot Boolean



% ---
[output_Match]  = RetrieveSS(Objectives,Constraints,GAoptions);

display(output_Match)
display(output_Match.Table)

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