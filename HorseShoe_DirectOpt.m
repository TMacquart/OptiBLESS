% ======================================================================= %
%                        HorseShoe Blending Example                       %
%
% fmincon(@EvaluationFct,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
% ======================================================================= %
clear all; close all force; clc; format short g; format compact;

addpath ./FitnessFcts
addpath ./src
addpath ./VisualGUI
addpath ./HorseShoeOpt
addpath ./StiffnessOpt

Objectives.Table   = [{'Laminate #'}  {'Nplies [LB UB]'}];
for i=1:18
    Objectives.Table = [Objectives.Table; {i}  {[20 40]}];
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
    
    
    % --- for gradients
    v21 = Parameters.v12*Parameters.E2/Parameters.E1;
    Q11 = Parameters.E1/(1-Parameters.v12*v21);
    Q22 = Parameters.E2/(1-Parameters.v12*v21);
    Q12 = Parameters.v12*Parameters.E2/(1-Parameters.v12*v21);
    Q66 = Parameters.G12;

end


Objectives.Type       = 'SS'; 
Objectives.FitnessFct = @(SS)  EvaluationFct(SS,Parameters);


% ---
Constraints.Vector     = [false       false          false          false         false            false            false];
Constraints.DeltaAngle = 45;
Constraints.ORDERED    = false;                         
Constraints.Balanced   = true; 
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
plotSS(Output,2)