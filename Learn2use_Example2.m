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

LpIni = [
    0.1821	 0.2102	 0.3000   % V1A
   -0.3643	-0.2871	-0.1732   % V2A
    0.0667	 0.1539	 0.1000   % V3A
   -0.1155	-0.1332	-0.1732   % V4A
    0.0000	 0.0000	 0.0000   % V1B
    0.0000	 0.0000	 0.0000   % V2B
    0.0000	 0.0000	 0.0000   % V3B
    0.0000	 0.0000	 0.0000   % V4B
    0.1699	 0.1227	-0.0741   % V1D
   -0.2584	-0.1579	 0.1131   % V2D
    0.2261	 0.3518	 0.4120   % V3D
   -0.3177	-0.2444	-0.3811]; % V4D
E1  = 13.0e9;
E2  = 72.0e9;
G12 = 26.9e9;
v12 = 0.33;
h   = 0.000127 * 10;

Objectives.mat = [E1 E2 G12 v12 h];

% convert lamination parameters into ABD matrices
for i = 1:size(LpIni,2)
    [A2Match{i},B2Match{i},D2Match{i}]    = Convert_LP2ABD (E1,E2,v12,G12,h,LpIni(:,i),true);                              %#ok<SAGROW>
    [Lp2Match{i},Abar{i},Bbar{i},Dbar{i}] = Convert_ABD2LP (E1,E2,v12,G12,h,A2Match{i},B2Match{i},D2Match{i},true);     %#ok<SAGROW>
end


Objectives.A = A2Match;
Objectives.B = B2Match;
Objectives.D = D2Match;

NPliesIni          = [20 22 12];
ImportanceFactor   = [1 1 1]; 
Objectives.IndexLP = [1:12];
Objectives.Table   = [{'Laminate Index'}      {'Nplies'}     {'LP2Match'}     {'Importance'} ;
                            {1}            {NPliesIni(1)}    Lp2Match{1}   {ImportanceFactor(1)} ;
                            {2}            {NPliesIni(2)}    Lp2Match{2}   {ImportanceFactor(2)} ;
                            {3}            {NPliesIni(3)}    Lp2Match{3}   {ImportanceFactor(3)} ; ];

                        
                        
% =========================== Default Options =========================== %

%                        [Damtol  Rule10percent  Disorientation  Contiguity   DiscreteAngle  InernalContinuity  Covering];
Constraints.Vector     = [true       false          false          false         true            false            true];
Constraints.DeltaAngle = 5;
Constraints.ply_t      = 0.000127;          % ply thickness
Constraints.ORDERED    = true;                           
Constraints.Balanced   = false; 
Constraints.Sym        = true; 
Constraints.NRange     = 1.0;


% ---
GAoptions.Npop    = 100; 	   % Population size
GAoptions.Ngen    = 200; 	   % Number of generations
GAoptions.NgenMin = 200; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.1; 	   % Percentage of elite passing to the next Gen.
GAoptions.Plot    = true; 	   % Plot Boolean


Objectives.Type        = 'ABD'
Objectives.FitnessFct = @(A,B,D) SumRMSABD(A,B,D,Objectives);


% ---
[output_Match]  = RetrieveSS(Objectives,Constraints,GAoptions);

display(output_Match)
display(output_Match.Table)

% --- Back from SS 2 ABD
for i = 1:length(output_Match.SS)
    [AOpt{i},BOpt{i},DOpt{i}] = Convert_SS2ABD(E1,E2,v12,G12,0.000127,output_Match.SS{i},true);                       %#ok<SAGROW>
end

A2Match{i}
AOpt{i}

B2Match{i}
BOpt{i}

D2Match{i}
DOpt{i}
