% ===                                                                   === 
%                        HorseShoe Blending Example                       
% ===                                                                   === 


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


clear all; close all force; clc; format short g; format compact;

addpath ./FitnessFcts
addpath ./src
addpath ./GUI
% addpath ./StiffnessOpt

Objectives.Table   = [{'Laminate #'}  {'Nplies [LB UB]'}];
for i=1:18
    Objectives.Table = [Objectives.Table; {i}  {[10 40]}];
end
                        
% ---
if 1    % Problem definition
    Parameters.Dim{1}      = [18 24];                                % Panel dimension (inch)
    Parameters.Dim{2}      = [20 12];                                % Panel dimension (inch)
    Parameters.PanelDim    = [1 1 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2]';
    Parameters.Nx          = [700 375 270 250 210 305 290 600 1100 900 375 400 330 190 300 815  320 300]';      % Paper Problem
    Parameters.Ny          = [400 360 325 200 100 360 195 480 600  400 525 320 330 205 610 1000 180 410]';      % Paper Problem
    Parameters.E1          = 20.5e6;
    Parameters.E2          = 1.31e6;
    Parameters.G12         = 0.62e6;
    Parameters.v12         = 0.32;
    Parameters.ply_t       = 0.0075;
    Parameters.ply_tMax    = 0.0075 *100;
    Parameters.rho         = 0.056999; % density (lb/inch^3) == 1577.727kg/m^3
    
    NormD = [12 18 20 24]/12;
    if 1
        PatchXYZ{1}.X = [0  1.5 1.5  0];
        PatchXYZ{1}.Y = [-2 -2  0    0];
        
        PatchXYZ{2}.X = 1.5+[0  1.5 1.5  0];
        PatchXYZ{2}.Y = [-2 -2  0    0];
        
        PatchXYZ{3}.X = 3+[0  1.6667 1.6667  0];
        PatchXYZ{3}.Y = [-1 -1  0    0];
        
        PatchXYZ{4}.X = 4.6667+[0  1.6667 1.6667  0];
        PatchXYZ{4}.Y = [-1 -1  0    0];
        
        PatchXYZ{5}.X = 6.3334+[0  1.6667 1.6667  0];
        PatchXYZ{5}.Y = [-1 -1  0    0];
        
        PatchXYZ{6}.X = 3+[0  1.6667 1.6667  0];
        PatchXYZ{6}.Y = [-2 -2  -1    -1];
        
        PatchXYZ{7}.X = 4.6667+[0  1.6667 1.6667  0];
        PatchXYZ{7}.Y = [-2 -2  -1    -1];
        
        PatchXYZ{8}.X = 6.3334+[0  1.6667 1.6667  0];
        PatchXYZ{8}.Y = [-2 -2  -1    -1];
        
        PatchXYZ{9}.X = [0  1.5 1.5  0];
        PatchXYZ{9}.Y = [-4 -4  -2    -2];
        
        PatchXYZ{10}.X = 1.5+[0  1.5 1.5  0];
        PatchXYZ{10}.Y = [-4 -4  -2    -2];
        
        PatchXYZ{11}.X = [0  1.5 1.5  0];
        PatchXYZ{11}.Y = [-6 -6  -4    -4];
        
        PatchXYZ{12}.X = 1.5+[0  1.5 1.5  0];
        PatchXYZ{12}.Y = [-6 -6  -4    -4];
        
        
        PatchXYZ{13}.X = 3+[0  1.6667 1.6667  0];
        PatchXYZ{13}.Y = -4+[-1 -1  0    0];
        
        PatchXYZ{14}.X = 4.6667+[0  1.6667 1.6667  0];
        PatchXYZ{14}.Y = -4+[-1 -1  0    0];
        
        PatchXYZ{15}.X = 6.3334+[0  1.6667 1.6667  0];
        PatchXYZ{15}.Y = -4+[-1 -1  0    0];
        
        PatchXYZ{16}.X = 3+[0  1.6667 1.6667  0];
        PatchXYZ{16}.Y = -4+[-2 -2  -1    -1];
        
        PatchXYZ{17}.X = 4.6667+[0  1.6667 1.6667  0];
        PatchXYZ{17}.Y = -4+[-2 -2  -1    -1];
        
        PatchXYZ{18}.X = 6.3334+[0  1.6667 1.6667  0];
        PatchXYZ{18}.Y = -4+[-2 -2  -1    -1];
        
    end
    
    for i=1:18
        PatchXYZ{i}.Z = [0 0 0 0];
    end
    
    Parameters.mMax = 1;
    Parameters.nMax = 1;
    
    Parameters.connectivity = [
        0 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
        1 0 1 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0
        0 1 0 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0
        0 0 1 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0
        0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0
        0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0
        0 0 1 1 1 1 0 1 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0
        1 1 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0
        1 1 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 1 1 1 0 1 0 0 1 0 0
        0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 1 0
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 1 1
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1
        0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 1 0
        0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 1
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 0];
    
end


Objectives.Type       = 'SS'; 
Objectives.FitnessFct = @(SS)  HS_EvaluationFct(SS,Parameters);


% ---
Constraints.Vector     = [false       false          false          false         false            false            false];
Constraints.DeltaAngle = 15;
Constraints.ORDERED    = false;                         
Constraints.Balanced   = false; 
Constraints.Sym        = true; 


% ---
GAoptions.Npop    = 100; 	   % Population size
GAoptions.Ngen    = 500; 	   % Number of generations
GAoptions.NgenMin = 500; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.05; 	   % Percentage of elite passing to the next Gen.
GAoptions.PC      = 0.75; 	   

GAoptions.PlotInterval = [10];                  % Refresh plot every X itterations         
GAoptions.SaveInterval = [];                  % Save Data every X itterations   
GAoptions.PlotFct      = @gaplotbestf;          % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;


% ---
[Output] = RetrieveSS(Objectives,Constraints,GAoptions);

% --- 
plotSS(Output,3,PatchXYZ)


% ---
Parameters.mMax = 2;
Parameters.nMax = 2;
[Fitness2,output2] = HS_EvaluationFct(Output.SS,Parameters)