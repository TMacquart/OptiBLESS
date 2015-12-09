% =====                                                              ==== 
%  This file is a typical user input file that is used to run the code.
% =====                                                              ==== 

clear all; clc; format short g; format compact; close all;

addpath ./FitnessFcts
addpath ./GUI
addpath ./src


%% === Objective

% --- Corresponding Staking Sequence
% --- Bottom ply [ 45   -45    90     0    45    90     0    45] Top ply
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


Objectives.Type    = 'LP'; 


ScalingCoef        = ones(12,1); 


Objectives.Table   = [{'Laminate #'}     {'Nplies'}      {'LP2Match'}     {'Scaling Coefficient'} ;
                            {1}           {[8 8]}         {Lp2Match(:,1)}  {ScalingCoef} ; ];
                        
                        
Objectives.FitnessFct = @(LP) RMSE_LP(LP,Objectives);



%% === Constraints 

%                        [Damtol  Rule10percent  Disorientation  Contiguity   BalancedIndirect  InernalContinuity  Covering];
Constraints.Vector     = [false       false          false          false         false            false            false];
Constraints.DeltaAngle = 45;
Constraints.ORDERED    = false;                         
Constraints.Balanced   = false; 
Constraints.Sym        = false; 
Constraints.UserFct    = false;


%% === Options 
GAoptions.Npop    = 100; 	   % Population size
GAoptions.Ngen    = 100; 	   % Number of generations
GAoptions.NgenMin = 100; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.01; 	   % Percentage of elite passing to the next Gen.
GAoptions.PC      = 0.75; 	   % Percentage of crossover

GAoptions.PlotInterval = [10];                  % Refresh plot every X itterations         
GAoptions.SaveInterval = [2];                   % Save Data every X itterations   
GAoptions.PlotFct      = @gaplotbestf;          % Refresh plot every X itterations
GAoptions.OutputFct    = @GACustomOutput;



%% === Run
[Output]  = RetrieveSS(Objectives,Constraints,GAoptions);

display(Output)
display(Output.Table)



%% === Plot
plotSS(Output)



%% === Checking output results are correct
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