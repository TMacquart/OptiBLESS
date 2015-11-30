% =====                                                              ==== 
%    Use GA to find stacking sequence ply angles matching Lam. Param.     
%          and varying thickness employing Stacking Sequence Table        
%
% [output] = RetrieveSS(Objectives,Constraints,GAoptions)
%
% 
%  The individual of GA is composed of 3 parts:
%  --------------------------------------------
%  [ [Nply(1) ... Nply(Npatch)]                                  -- the Number of plies
%  [ Theta(1) ... Theta(N) ]                                     -- N is the guide Nply
%  [ (Location of -Theta(1)) ... (Location of -Theta(N)) ]       -- location of balanced angle pairs
%  [ Drop(1)  ... Drop(M)  ]                                     -- M is the Delta Nply
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


function [output] = RetrieveSS(Objectives,Constraints,GAoptions)


%% Treat Inputs  
[Nvar,NpatchVar,NthetaVar,NdropVar,LamType,LB,UB,AllowedNplies] = FormatInput(Objectives,Constraints);
% Nvar          - Number of design variables 
% NpatchVar     - Number of patches with variable number of plies
% NthetaVar     - Number of fibre angles used to describe the guide laminate
% NdropVar      - Number of drops from the guide laminate
% LamType       - Type of laminates
% LB            - Lower bound of design variables
% UB            - Upper bound of design variables
% AllowedNplies - Number of plies allowed for each patches





%% Set GA, see --- doc gaoptimset --- for more option
options  = gaoptimset('PopulationSize',GAoptions.Npop,...
                      'Generation',GAoptions.Ngen, ...
                      'StallGenLimit',GAoptions.NgenMin, ...
                      'EliteCount',ceil(GAoptions.Elitism*GAoptions.Npop),...
                      'FitnessLimit' ,1e-5,...
                      'TolFun' ,1e-10,...
                      'PlotInterval',1,...
                      'CrossoverFraction',GAoptions.PC);
                  
if GAoptions.Plot
    options  = gaoptimset(options,'PlotFcns' ,{@gaplotbestf});
end


% Handle of the fitness function 
fct_handle = @(x)Eval_Fitness(x,Objectives,Constraints,NpatchVar,NthetaVar,AllowedNplies,LamType);  



%% Generate Ini. Pop.
for i = 1:5
    try
        [IniPop] = Generate_IniPop (Nvar,GAoptions.Npop,NpatchVar,NthetaVar,NdropVar,Constraints,AllowedNplies,LamType);
        break; 
    catch
        fprintf('Inipop Failed. Retrying ...\n');
        if i == 5
            error(' Did not manage to generate a feasible initial population. Try relaxing some constraints.')
        end
    end
end
%  IniPop(1,:) = [(90+[-45 0 45 90 0  -45  45  90  -45  45])/Constraints.DeltaAngle [1 3 6 7]]
% keyboard
%  IniPop(1,:) = [(90+[-45 0 45 90 0 -45  90 45 0 -45 0 45])/Constraints.DeltaAngle [2 11 6 12]]
options = gaoptimset(options,'InitialPopulation' ,IniPop);

% keyboard
% [Fval,output]=fct_handle([(90+[-45 0 45 90 0 -45  90 45 0 -45 0 45])/Constraints.DeltaAngle [2 11 6 12]])

%% run GA
fprintf(strcat('Running GA \n'))
[xOpt,fval,~,StandardOutputGA] = ga(fct_handle,Nvar,[],[],[],[],LB,UB,[],1:Nvar,options);
display('GA(s) Terminated Successfully')



%% Format Results
[~,output] = fct_handle(xOpt);

SS = output.SS;

if strcmp(Objectives.Type,'LP')
    LPMatched = output.LP;
    Table = [{'Lam #'} {'Nplies SST'} {'Ply Angles'} {'LP2Match'} {'LPOpt'} {'NormE'} {'RMSE'} {'MAE'} {'MaxAE'}];
    for j = 1:length(AllowedNplies)
        LP2Match    = Objectives.Table{j+1,3};
        ScalingCoef = Objectives.Table{j+1,4};
        
        
        QualIndex1 = norm( (LPMatched(:,j) - LP2Match(:)).*ScalingCoef );
        QualIndex2 = rms( (LPMatched(:,j) - LP2Match(:)).*ScalingCoef );
        QualIndex3 = mae( (LPMatched(:,j) - LP2Match(:)).*ScalingCoef );
        QualIndex4 = max( abs((LPMatched(:,j) - LP2Match(:)).*ScalingCoef) );
        Table      = [Table ;  {j} {length(SS{j})} SS(j) {LP2Match} {LPMatched(:,j)} {QualIndex1} {QualIndex2} {QualIndex3} {QualIndex4}]; %#ok<AGROW>
    end
end


if strcmp(Objectives.Type,'ABD')

    Table = [{'Lam #'} {'Nplies SST'} {'Ply Angles'} {'A2Match'} {'AOpt'} {'Error % A'} {'Error Norm A'} {'Error RMS A'} ...
            {'B2Match'} {'BOpt'} {'Error % B'} {'Error Norm B'} {'Error RMS B'} ...
            {'D2Match'} {'DOpt'} {'Error % D'} {'Error Norm D'} {'Error RMS D'}];
    for j = 1:length(AllowedNplies)
        A_Matched = output.A{j};
        B_Matched = output.B{j};
        D_Matched = output.D{j};
        A2Match   = Objectives.Table{j+1,3};
        B2Match   = Objectives.Table{j+1,4};
        D2Match   = Objectives.Table{j+1,5};
        
        AScaling = Objectives.Table{j+1,6};
        BScaling = Objectives.Table{j+1,7};
        DScaling = Objectives.Table{j+1,8};
    
        QualIndex1A = 100*sum(abs(  AScaling(:).*((A_Matched(:) - A2Match(:))./A2Match(:)) ));
        QualIndex2A = norm( AScaling(:).*(A_Matched(:) - A2Match(:)) );
        QualIndex3A = rms(  AScaling(:).*(A_Matched(:) - A2Match(:)) );
        
        QualIndex1B = 100*sum(abs(  BScaling(:).*((B_Matched(:) - B2Match(:))./B2Match(:)) ));
        QualIndex2B = norm( BScaling(:).*(B_Matched(:) - B2Match(:)) );
        QualIndex3B = rms(  BScaling(:).*(B_Matched(:) - B2Match(:)) );
        
        QualIndex1D = 100*sum(abs(  DScaling(:).*((D_Matched(:) - D2Match(:))./D2Match(:)) ));
        QualIndex2D = norm( DScaling(:).*(D_Matched(:) - D2Match(:)) );
        QualIndex3D = rms(  DScaling(:).*(D_Matched(:) - D2Match(:)) );
        
        Table = [Table ;  {j} {length(SS{j})} SS(j) ...                                      
                    {A2Match} {A_Matched} {QualIndex1A} {QualIndex2A} {QualIndex3A}...
                    {B2Match} {B_Matched} {QualIndex1B} {QualIndex2B} {QualIndex3B} ...
                    {D2Match} {D_Matched} {QualIndex1D} {QualIndex2D} {QualIndex3D}];               %#ok<AGROW>
    end

end


output.NfctEval  = StandardOutputGA.funccount;
output.NGen      = StandardOutputGA.generations;
output.Table     = Table;
output.xOpt      = xOpt;
output.fval      = fval;


if ~output.FEASIBLE,    warning('The SST GA could not find a single feasible solution!'); end
end
