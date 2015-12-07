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
clear all; clc; format short g; format compact; close all Force; fclose all;

addpath ./src
addpath ./FitnessFcts

% Enter any Guide laminates and drops to test the code
GuideLamDv  = [45 15 -5 0 -45 -60 45 90 -75 80 50 20 -10 0 90 45 -45 -20 -40 0];                          
Drops       = [{[2 5]} {[12 18]} {[10]} {[20 19 1 4 6 7]}]; 
ScalingCoef = [1 1 1 1, 1 1 1 1, 1 1 1 1]';      % relative importance given to matching the guide laminate LPs integer [1,N], the higher the integer = the more impact on the fit. fct.


% --- Automated LP2Match calculation (Do Not Change)
NUniqueLam         = length(Drops)+1;
Lp2Match           = zeros(12,NUniqueLam);
Objectives.Table   = [{'Laminate Index'} {'Nplies'} {'LP2Match'} {'Importance'}];

for i = 1:NUniqueLam
    Lam = GuideLamDv;
    DropsLoc = cell2mat(Drops(1:i-1));
    if i ~= 1,  
        Lam(DropsLoc) = [];   
    end
    
    Lp2Match(:,i)    = Convert_SS2LP(Lam);
    Objectives.Table = [Objectives.Table; [{i} {[1 1]*length(Lam)} {Lp2Match(:,i)} {ScalingCoef}]];
end    
% ---

Objectives.Type       = 'LP';
% Objectives.FitnessFct = @(LP) RMSE_MaxAE_LP(LP,Objectives);
Objectives.FitnessFct = @(LP) RMSE_LP(LP,Objectives);

% =========================== Default Options =========================== %

%                        [Damtol  Rule10percent  Disorientation  Contiguity   BalancedIndirect   InernalContinuity  Covering];
Constraints.Vector     = [false       false          false          false             false            false            false];
Constraints.DeltaAngle = 5;
Constraints.ORDERED    = false;                           
Constraints.Balanced   = false; 
Constraints.Sym        = false; 


% ---
GAoptions.Npop     = 150;      % Population size
GAoptions.Ngen     = 2000;     % Number of generations
GAoptions.NgenMin  = 2000;     % Minimum number of generation calculated
GAoptions.Elitism  = 0.075;   % Percentage of elite passing to the next Gen.
GAoptions.PC       = 0.75;    % Plot Boolean

GAoptions.PlotInterval = [10];                  % Refresh plot every X itterations         
GAoptions.SaveInterval = [];                  % Save Data every X itterations   
GAoptions.PlotFct      = @gaplotbestf;          % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;



% ---
[Output]  = RetrieveSS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)
