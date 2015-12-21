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
addpath ./src/StiffnessOpt

Optimisation = 0; % 0=direct, 1=indirect


if Optimisation ==0
    Objectives.Table   = [{'Laminate #'}  {'Nplies [LB UB]'}];
    for i=1:18
        Objectives.Table = [Objectives.Table; {i}  {[18 48]}];
    end
else
    Table =[ % Nply  V1D     V3D
                32	0.208	-0.843 
                28	0.092	-0.714
                20	-0.722	0.054
                18	-0.582	-0.228
                16	-0.477	-0.235
                22	-0.469	-0.335
                18	-0.582	-0.228
                24	-0.597	-0.252
                38	0.192	-0.657
                34	0.308	-0.776
                30	-0.241	-0.816
                28	0.092	-0.714
                22	-0.469	-0.335
                18	-0.582	-0.228
                24	-0.597	-0.252
                30	-0.241	-0.816
                18	-0.582	-0.228
                22	-0.469	-0.335];
    
            ScalingCoef = [0 0 0 0, 0 0 0 0, 1 0 1 0]';
            Objectives.Table   = [{'Laminate Index'} {'Nplies [LB UB]'} {'LP2Match'} {'Importance'}];
            for i=1:18
                Lp2Match(:,i) = [zeros(8,1); Table(i,2); 0; Table(i,3);0 ];
                Objectives.Table = [Objectives.Table; {i}  {[1 1]*Table(i,1)} {Lp2Match(:,i)}  ScalingCoef];
            end
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


Objectives.UserFct    = true;
if Optimisation ==0
    Objectives.Type       = 'SS';
    Objectives.FitnessFct = @(SS)  HS_EvaluationFct(SS,Parameters);
else
    Objectives.Type       = 'LP';
    Objectives.FitnessFct = @(LP)  RMSE_LP(LP,Objectives);
end

% ---
%                        [Damtol  Rule10percent  Disorientation  Contiguity   BalancedIndirect  InernalContinuity  Covering];
Constraints.Vector     = [false       false          false          false         false            false            false];
Constraints.DeltaAngle = 15;
Constraints.ORDERED    = false;                         
Constraints.Balanced   = true; 
Constraints.Sym        = true; 
Constraints.PatchXYZ   = PatchXYZ;


% --- Format Geometric Input
if 1
    if isfield(Constraints,'PatchXYZ') && ~isfield(Constraints,'PatchConnectivity')
        PatchConnectivity = Format_GeometricInput(Constraints.PatchXYZ);
        display(PatchConnectivity)
        UserInput = input(' Please check the Patch Connectivity matrix. Do you want to continue? [Y/N]: ','s');
        if ~strcmp(UserInput,'Y')
            error('Stopped by user. If the automatic Patch Connectivity matrix is incorrect you can input it directly as Constraints.PatchConnectivity')
        end
        Constraints.PatchConnectivity = PatchConnectivity;
    end
end


% ---
GAoptions.Npop    = 25; 	   % Population size
GAoptions.Ngen    = 250; 	   % Number of generations
GAoptions.NgenMin = 250; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.05; 	   % Percentage of elite passing to the next Gen.
GAoptions.PC      = 0.75; 	   

GAoptions.PlotInterval = [10];                  % Refresh plot every X itterations         
GAoptions.SaveInterval = [];                  % Save Data every X itterations   
GAoptions.PlotFct      = @gaplotbestf;          % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;


% ---
[Output] = RetrieveSS(Objectives,Constraints,GAoptions);

% --- 
plotSS(Output,PatchXYZ)

NGeoConstraints = CheckContinuity(Output.SS,Constraints.PatchConnectivity)


% ---
Parameters.mMax = 2;
Parameters.nMax = 2;
[Fitness2,output2] = HS_EvaluationFct(Output.SS,Parameters)