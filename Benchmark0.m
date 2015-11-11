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


% --- Creating a real SS to match
Objectives.Table   = [{'Laminate Index'} {'Nplies'} {'LP2Match'} {'Importance'}];
                        
ply_t      = 0.000127;
GuideLamDv = [+45 -45 90 0 45 90 0 45];
Drops      = [2 4];
GuideLam   = [GuideLamDv, -GuideLamDv, fliplr([GuideLamDv, -GuideLamDv])]'; % balanced/symetric

ImportanceFactor = [1 1 1];     % relative importance given to matching the guide laminate LPs integer [1,N], the higher the integer = the more impact on the fit. fct.
NUniqueLam = length(Drops)+1;
Lp2Match   = zeros(12,NUniqueLam);
for i = 1:NUniqueLam
    Lam = GuideLamDv;
    if i ~= 1,        
        Lam(Drops(1:i-1)) = [];   
    end
    Lam = [Lam, -Lam, fliplr([Lam, -Lam])]'; % balanced/symetric 
    
    Lp2Match(:,i) = SS2LP(Lam);
    Objectives.Table = [Objectives.Table; [{i} {length(Lam)} {Lp2Match(:,i)} {ImportanceFactor(i)}]];
end

Objectives.IndexLP = [1 3 9 10 11 12];



                        
% =========================== Default Options =========================== %

%                        [Damtol  Rule10percent  Disorientation  Contiguity   DiscreteAngle  InernalContinuity  Covering];
Constraints.Vector     = [false       false          false          false         true            false            false];
Constraints.DeltaAngle = 45;
Constraints.ply_t      = ply_t;          % ply thickness
Constraints.Balanced   = true; 
Constraints.Sym        = true; 
Constraints.NRange     = 1;
Constraints.ORDERED    = true;           
Constraints.alpha      = 0;               

% ---
GAoptions.Npop    = 100; 	   % Population size
GAoptions.Ngen    = 250; 	   % Number of generations
GAoptions.NgenMin = 250; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.05; 	   % Percentage of elite passing to the next Gen.
GAoptions.Plot    = true; 	   % Plot Boolean

% ---
Nrun = 10;
output_Match = cell(1,Nrun);
feasible     = zeros(1,Nrun);
fval         = zeros(1,Nrun);
for i = 1:Nrun
    [output_Match{i}]  = RetrieveSS(Objectives,Constraints,GAoptions);
    feasible(i)        = output_Match{i}.FEASIBLE;
    fval(i)            = output_Match{i}.fval;
    fitRMS(i)          = output_Match{i}.fitRMS;
end
display(feasible)
display(fval)

