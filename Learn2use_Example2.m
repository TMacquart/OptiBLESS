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
addpath ./src
addpath ./src/StiffnessOpt
addpath ./GUI

% --- bottom [ 45   -45    90     0    45    90     0    45] Top   
E1   = 13.0e9;
E2   = 72.0e9;
G12  = 26.9e9;
v12  = 0.33;

tply = 0.000127;  % ply thickness
h    = 8*tply;

% --- Stacking sequence corresponding stiffness matrices

A2Match ={[
   1.0874e+11   5.8225e+10  -9.2917e+09
   5.8225e+10   1.0874e+11  -9.2917e+09
  -9.2917e+09  -9.2917e+09   2.5255e+10]};

B2Match ={[
  -9.7029e+09   4.1122e+08  -6.9687e+09
   4.1122e+08   8.8804e+09  -6.9687e+09
  -6.9687e+09  -6.9687e+09   4.1122e+08]};

D2Match ={[
   1.0602e+11   5.7454e+10   -1.626e+10
   5.7454e+10   1.1299e+11   -1.626e+10
   -1.626e+10   -1.626e+10   2.4484e+10]};

% --- ( you can used [A,B,D]=Convert_SS2ABD(E1,E2,v12,G12,tply,[ 45  -45    90     0    45    90     0    45],true) 
% --- to convert a stacking sequence into lamination parameters)

Objectives.mat = [E1 E2 G12 v12 h];
 
IndexAStiff = ones(3,3);
IndexBStiff = ones(3,3);
IndexDStiff = ones(3,3);

Objectives.Table   = [{'Laminate #'}     {'Nplies'}   {'A2Match'}  {'B2Match'} {'D2Match'}  {'A Scaling'} {'B Scaling'} {'D Scaling'} ;
                            {1}          {[8 8]}       A2Match       B2Match     D2Match   {IndexAStiff} {IndexBStiff} {IndexDStiff}];

                        
Objectives.Type       = 'ABD';
Objectives.UserFct    = false;
Objectives.FitnessFct = @(A,B,D) RMSE_ABD(A,B,D,Objectives);


%% =========================== Default Options =========================== %

%                        [Damtol  Rule10percent  Disorientation  Contiguity   BalancedIndirect  InernalContinuity  Covering];
Constraints.Vector     = [false       false          false          false         false            false            false];
Constraints.DeltaAngle = 45;
Constraints.ply_t      = 0.000127;          % ply thickness
Constraints.ORDERED    = true;                         
Constraints.Balanced   = false; 
Constraints.Sym        = false; 





% ---
GAoptions.Npop    = 100; 	   % Population size
GAoptions.Ngen    = 500; 	   % Number of generations
GAoptions.NgenMin = 500; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.01; 	   % Percentage of elite passing to the next Gen.
GAoptions.PC      = 0.75; 	   % Percentage of crossover

GAoptions.PlotInterval = [10];                  % Refresh plot every X itterations         
GAoptions.SaveInterval = [1];                  % Save Data every X itterations   
GAoptions.PlotFct      = @gaplotbestf;          % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;


% ---
[Output]  = RetrieveSS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)

%% Plot
plotSS(Output)

%% --- Back from SS 2 ABD
for i = 1:length(Output.SS)
    [AOpt{i},BOpt{i},DOpt{i}] = Convert_SS2ABD(E1,E2,v12,G12,tply,Output.SS{i},true);                       %#ok<SAGROW>
end

A2Match{i}
AOpt{i}

B2Match{i}
BOpt{i}

D2Match{i}
DOpt{i}
