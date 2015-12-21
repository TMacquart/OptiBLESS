% =====                                                              ==== 
%  This file is a typical user input file that is used to run the code.
% =====                                                              ==== 

clear all; clc; format short g; format compact; close all;

addpath ./src
addpath ./FitnessFcts
addpath ./GUI

% GuideLamDv = [-45 0 45 90 0  -45  45  90  -45  45];    
% Lam1       = [    0    90 0  -45  45  90  -45  45];   
% Lam2       = [    0    90 0           90  -45  45];   
% Lam3       = [    0    90 0           90         ];

Lp2Match = [
    % Guide     % Lam 1     % Lam 2          % Lam 3
        0          0          0              0          % V1A
        0          0          0              0          % V2A
       -0.2        0          0.33333        1          % V3A
        0          0          0              0          % V4A
       -0.2       -0.25      -0.22222       -0.5        % V1B
        0.16       0.125      0.11111        0          % V2B
       -0.24      -0.75      -0.88889        0          % V3B
        0          0          0              0          % V4B
        0.048      0.14062    0.22222        0          % V1D
       -0.048      0.14062    0.22222        0          % V2D
       -0.488      0.09375    0.037037       1          % V3D
        0          0          0              0 ];       % V4D


Objectives.Type   = 'LP';
ScalingCoef       = [1 1 1 1, 1 1 1 1, 1 1 1 1]';
Objectives.Table  = [{'Laminate #'}  {'Nplies [LB UB]'}     {'LP2Match'}     {'Scaling Coefficient'} ;
                            {1}           {[10 10]}         Lp2Match(:,1)     {ScalingCoef} ;
                            {2}           {[8 8]}           Lp2Match(:,2)     {ScalingCoef} ; 
                            {3}           {[6 6]}           Lp2Match(:,3)     {ScalingCoef} ;
                            {4}           {[4 4]}           Lp2Match(:,4)     {ScalingCoef} ];

Objectives.UserFct    = false;                        
Objectives.FitnessFct = @(LP) RMSE_LP(LP,Objectives);


if 1 % --- Format Geometric Input
    Patch{1}.X = [0  1.5 1.5  0];
    Patch{1}.Y = [-2 -2  0    0];
    Patch{1}.Z = [0   0  0    0];
    
    Patch{2}.X = 1.5 +[0  3 3 0];
    Patch{2}.Y = [-2 -2  -1  -1];
    Patch{2}.Z = [0   0  0    0];
    
    Patch{3}.X = 3+[0  1.5 1.5  0];
    Patch{3}.Y = [-1 -1  0    0];
    Patch{3}.Z = [0   0  0    0];
    
    Patch{4}.X = 4.5+[0   1.5 1.5  0];
    Patch{4}.Y = [-1 -1  0    0];
    Patch{4}.Z = [0   0  0    0];
    
    Constraints.Patch = Patch;
    
    if ~isfield(Constraints,'PatchConnectivity')
        PatchConnectivity = Format_GeometricInput(Constraints.Patch);
        display(PatchConnectivity)
        UserInput = input(' Please check the Patch Connectivity matrix. Do you want to continue? [Y/N]: ','s');
        if ~strcmp(UserInput,'Y')
            error('Stopped by user. If the automatic Patch Connectivity matrix is incorrect you can input it directly as Constraints.PatchConnectivity')
        end
        Constraints.PatchConnectivity = PatchConnectivity;
    end
end



% =========================== Default Options =========================== %

%                        [Damtol  Rule10percent  Disorientation  Contiguity   BalancedIndirect   InernalContinuity  Covering];
Constraints.Vector     = [false       false          false          false             false            false            false];
Constraints.DeltaAngle = 15;
Constraints.ORDERED    = false;
Constraints.Balanced   = true;
Constraints.Sym        = false;


% ---
GAoptions.Npop    = 150; 	   % Population size
GAoptions.Ngen    = 500; 	   % Number of generations
GAoptions.NgenMin = 500; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.075; 	   % Percentage of elite passing to the next Gen.
GAoptions.PC      = 0.75; 	   % Plot Boolean

GAoptions.PlotInterval = [10];                  % Refresh plot every X itterations
GAoptions.SaveInterval = [];                    % Save Data every X itterations
GAoptions.PlotFct      = @gaplotbestf;          % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;


% ---
[Output]  = RetrieveSS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)


%% Plot
% plotSS(Output)
plotSS(Output,Patch)


%% Ouput Check (Optional Validation)
ScalingCoef = reshape(cell2mat(Objectives.Table(2:end,4)),12,size(Objectives.Table,1)-1);

for i = 2:size(Objectives.Table,1)
    LP2Match = Objectives.Table{i,3};
    if sum(abs(LP2Match-Output.Table{i,4}))>1e-10
        error('non matching LP2Match')
    end
    
    LP = Convert_SS2LP(Output.Table{i,3});
    if sum(abs(LP-Output.Table{i,5}))>1e-10
        error('non matching SS and LPOpt')
    end
    
    if abs( rms ( (LP-LP2Match).*ScalingCoef(:,i-1) )-Output.Table{i,7})>1e-10
        error('non matching RSM')
    end
    
    if abs( norm ((LP-LP2Match).*ScalingCoef(:,i-1))-Output.Table{i,6})>1e-10
        error('non matching norm')
    end
    
end
