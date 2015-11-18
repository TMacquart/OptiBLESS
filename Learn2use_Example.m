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

% ---
Lp2Match = [
% LP2Match1 LP2Match2  LP2Match3
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

NPliesIni = [28 20 16];
ScalingCoef = [1 1 1]; 
Objectives.IndexLP = [1 3];
Objectives.Table   = [{'Laminate #'}     {'Nplies'}      {'LP2Match'}     {'Scaling Coefficient'} ;
                            {1}            {NPliesIni(1)}    {Lp2Match(:,1)}   {ScalingCoef(1)} ;
                            {2}            {NPliesIni(2)}    {Lp2Match(:,2)}   {ScalingCoef(2)} ;
                            {3}            {NPliesIni(3)}    {Lp2Match(:,3)}   {ScalingCoef(3)} ; ];

                        
                        
% =========================== Default Options =========================== %

%                        [Damtol  Rule10percent  Disorientation  Contiguity   DiscreteAngle  InernalContinuity  Covering];
Constraints.Vector     = [false       true          true          true         true            true            false];
Constraints.DeltaAngle = 5;
Constraints.ply_t      = 0.000127;          % ply thickness
Constraints.ORDERED    = false;                         
Constraints.Balanced   = false; 
Constraints.Sym        = false; 
Constraints.NRange     = 1.2;

Objectives.Type        = 'LP'; % 'ABD' 'SS' 'LP'
Objectives.FitnessFct = @(LP) SumRMSLP(LP,Objectives);

% ---
GAoptions.Npop    = 100; 	   % Population size
GAoptions.Ngen    = 250; 	   % Number of generations
GAoptions.NgenMin = 250; 	   % Minimum number of generation calculated
GAoptions.Elitism = 0.05; 	   % Percentage of elite passing to the next Gen.
GAoptions.Plot    = true; 	   % Plot Boolean
GAoptions.PC      = 0.5; 	   


% ---
[output_Match]  = RetrieveSS(Objectives,Constraints,GAoptions);

display(output_Match)
display(output_Match.Table)

%% Checking output results are correct
IndexLp = Objectives.IndexLP;
for i = 2:length(NPliesIni)+1
    LP2Match = Objectives.Table{i,3};
    if sum(abs(LP2Match-output_Match.Table{i,5}))>1e-10
        error('non matching LP2Match')
    end
    
    LP = Convert_SS2LP(output_Match.Table{i,4});
    if sum(abs(LP-output_Match.Table{i,6}))>1e-10
        error('non matching SS and LPOpt')
    end
    
    if abs( rms (LP(IndexLp)-LP2Match(IndexLp))-output_Match.Table{i,9})>1e-10
        error('non matching RSM')
    end
    
    if abs( norm (LP(IndexLp)-LP2Match(IndexLp))-output_Match.Table{i,8})>1e-10
        error('non matching norm')
    end
    
    if abs( 100*sum(abs(  (LP(IndexLp) - LP2Match(IndexLp))./LP2Match(IndexLp) )) - output_Match.Table{i,7})>1e-10
        error('non matching error percent')
    end
end