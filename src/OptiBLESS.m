% =====                                                              ==== 
%   Use GA to optimise stacking sequence ply angles, Drops and Number of      
%     plies. Direct and Indirect Optimisation Framework can be used        
%
% [output] = OptiBLESS(Objectives,Constraints,GAoptions)
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


function [output] = OptiBLESS(Objectives,Constraints,GAoptions)


%% Format Inputs  

[NStruct,LamType,LB,UB,BCs,AllowedNplies] = FormatInput(Objectives,Constraints);
% Nvar          - Number of design variables 
% NpatchVar     - Patches with variable number of plies (set to 1 if variable)
% NthetaVar     - Number of fibre angles used to describe the guide laminate
% NdropVar      - Number of drops from the guide laminate
% LamType       - Type of laminates
% LB            - Lower bound of design variables
% UB            - Upper bound of design variables
% AllowedNplies - Number of plies allowed for each patches

%% Set GA, see --- doc gaoptimset --- for more option
 options  = gaoptimset('PopulationSize',GAoptions.Npop,...
                      'Generation',GAoptions.Ngen, ...
                      'StallGenLimit',GAoptions.NgenMin, ...                    % Minimum Number of Generation computed
                      'EliteCount',ceil(GAoptions.Elitism*GAoptions.Npop),...   % Elitism
                      'FitnessLimit' ,1e-5,...                                  % Stoping fitness criterion
                      'TolFun' ,1e-10,...                                       % Stoping change in fitness criterion
                      'CrossoverFraction',GAoptions.PC,...                      % crossover fraction
                      'PlotInterval',GAoptions.PlotInterval);
                  
if ~isempty(GAoptions.SaveInterval)
    options  = gaoptimset(options,'OutputFcns',{@(options,state,flag)GAoptions.OutputFct(options,state,flag,GAoptions.SaveInterval,Objectives.Type)});  % Output function used  
end                                                  
if ~isempty(GAoptions.PlotInterval)
    options  = gaoptimset(options,'PlotFcns',{GAoptions.PlotFct});              % Plot function used  
end

    

% Handle of the fitness function 
fct_handle = @(x)Eval_Fitness(x,Objectives,Constraints,NStruct,AllowedNplies,LamType);  


keyboard

for i=1:length(UB)
    Individual(i,1) = randi([LB(i) UB(i)],1,1);
end
fct_handle(Individual)


keyboard
%% Generate Initial Population
for i = 1:5
    try
        [IniPop] = Generate_IniPop (NStruct,GAoptions.Npop,Constraints,Objectives,AllowedNplies,LamType,BCs,fct_handle);
        break; 
    catch
        fprintf('Inipop Failed. Retrying ...\n');
        if i == 5
            error('Did not manage to generate a feasible initial population. There might be something wrong.')
        end
    end
end   
% IniPop(1,:) = [(90+[10   -65   -60   -40    65   -40    60   -45    80   -25])/Constraints.DeltaAngle [] []] % not balanced
% IniPop(1,:) = [40*ones(1,18) zeros(1,70)]
% keyboard
options = gaoptimset(options,'InitialPopulation' ,IniPop);


%% run GA
display('Running GA')

% fct_handle(IniPop(1,:))
[xOpt,fval,~,OutputGA] = ga(fct_handle,Nvar,[],[],[],[],LB,UB,[],1:Nvar,options);

display('GA(s) Terminated Successfully')


[~,output]       = fct_handle(xOpt);                                        % Evaluate the best individual found during GA, returns the output structure
output.NfctEval  = OutputGA.funccount;                                      % Number of function evaluation that have been computed
output.NGen      = OutputGA.generations;                                    % Number of generation computed
output.xOpt      = xOpt;                                                    % Genotype of the best found individual
output.fval      = fval;                                                    % Fintess value of the best found individual
output.LamType   = LamType;

if ~output.FEASIBLE,  
    warning('The optimal solution found is not feasible!');
end



%% Add output Lamination Parameter Results
if strcmp(Objectives.Type,'LP')
    Table     = [{'Lam #'} {'Nplies'} {'Ply Angles'} {'LP2Match'} {'LP Retrieved'} {'NormE'} {'RMSE'} {'MAE'} {'MaxAE'}];
    LPMatched = output.LP;                                                          % Lamination parameters retrieved by the GA
    
    for j = 1:length(AllowedNplies)
        LP2Match    = Objectives.Table{j+1,3};                                      % Lamination parameters given as objectives 
        ScalingCoef = Objectives.Table{j+1,4};                                      % Scaling coefficients given as objectives 
        
        QualIndex1 = norm( (LPMatched(:,j) - LP2Match(:)).*ScalingCoef );           % Norm Error
        QualIndex2 = rms( (LPMatched(:,j) - LP2Match(:)).*ScalingCoef );            % Root mean square error
        QualIndex3 = mae( (LPMatched(:,j) - LP2Match(:)).*ScalingCoef );            % Mean absolute error
        QualIndex4 = max( abs((LPMatched(:,j) - LP2Match(:)).*ScalingCoef) );       % Maximum absolute error
        
        Table = [Table ;  {j} {length(output.SS{j})} output.SS(j) {LP2Match} {LPMatched(:,j)} {QualIndex1} {QualIndex2} {QualIndex3} {QualIndex4}]; %#ok<AGROW>
    end
    
    output.Table     = Table;                                                   % Table sumarising results
end



%% Add output Stiffness Results
if strcmp(Objectives.Type,'ABD')

    
    Table = [{'Lam #'} {'Nplies SST'} {'Ply Angles'} {'A2Match'} {'AOpt'} {'Error % A'} {'Error Norm A'} {'Error RMS A'} ...
            {'B2Match'} {'BOpt'} {'Error % B'} {'Error Norm B'} {'Error RMS B'} ...
            {'D2Match'} {'DOpt'} {'Error % D'} {'Error Norm D'} {'Error RMS D'}];
    for j = 1:length(AllowedNplies)
        A_Matched = output.A{j};                                            % In-plane stiffness matrix retrieved by the GA
        B_Matched = output.B{j};                                            % Coupling stiffness matrix gretrieved by the GA
        D_Matched = output.D{j};                                            % Out-of-plane stiffness matrix retrieved by the GA
        A2Match   = Objectives.Table{j+1,3};                                % In-plane stiffness matrix given as objectives
        B2Match   = Objectives.Table{j+1,4};                                % Coupling stiffness matrix given as objectives
        D2Match   = Objectives.Table{j+1,5};                                % Out-of-plane stiffness matrix given as objectives
        
        AScaling = Objectives.Table{j+1,6};                                 % In-plane scaling coefficients
        BScaling = Objectives.Table{j+1,7};                                 % Coupling scaling coefficients
        DScaling = Objectives.Table{j+1,8};                                 % Out-of-plane scaling coefficients
    
        QualIndex1A = 100*sum(abs(  AScaling(:).*((A_Matched(:) - A2Match(:))./A2Match(:)) ));
        QualIndex2A = norm( AScaling(:).*(A_Matched(:) - A2Match(:)) );
        QualIndex3A = rms(  AScaling(:).*(A_Matched(:) - A2Match(:)) );
        
        QualIndex1B = 100*sum(abs(  BScaling(:).*((B_Matched(:) - B2Match(:))./B2Match(:)) ));
        QualIndex2B = norm( BScaling(:).*(B_Matched(:) - B2Match(:)) );
        QualIndex3B = rms(  BScaling(:).*(B_Matched(:) - B2Match(:)) );
        
        QualIndex1D = 100*sum(abs(  DScaling(:).*((D_Matched(:) - D2Match(:))./D2Match(:)) ));
        QualIndex2D = norm( DScaling(:).*(D_Matched(:) - D2Match(:)) );
        QualIndex3D = rms(  DScaling(:).*(D_Matched(:) - D2Match(:)) );
        
        Table = [Table ;  {j} {length(output.SS{j})} output.SS(j) ...                                      
                    {A2Match} {A_Matched} {QualIndex1A} {QualIndex2A} {QualIndex3A}...
                    {B2Match} {B_Matched} {QualIndex1B} {QualIndex2B} {QualIndex3B} ...
                    {D2Match} {D_Matched} {QualIndex1D} {QualIndex2D} {QualIndex3D}];               %#ok<AGROW>
    end
    
    output.Table     = Table;                                                   % Table sumarising results
end




