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


% =====                                                              ==== %
%    Use GA to find stacking sequence ply angles matching Lam. Param.     %
%          and varying thickness employing Stacking Sequence Table        %
%                       Terence Macquart (09/11/2015)                     %
%
%  Blended laminate
%  ----------------
%  objective:   match other LP by removing plies while minimising thickness 
%  constraints implemented: covering, max stoping, continuity, alternation
%
% 
%  The individual of GA is composed of 3 parts:
%  --------------------------------------------
%  [ [Nply(1) ... Nply(Nsec)]                   -- the Number of plies
%  [ Theta(1) ... Theta(N)   nan ... nan]       -- N is the guide Nply
%  [ Drop(1)  ... Drop(M)    nan ... nan]       -- M is the Delta Nply
% =====                                                              ==== %

function [output] = RetrieveSS(Objectives,Constraints,GAoptions)


% ---
if 1    % Define objectives (i.e. LP to match)
    Obj.LamInd    = cell2mat(Objectives.Table(2:end,1));
    Obj.Nplies    = cell2mat(Objectives.Table(2:end,2));
    Obj.Lps2Match = Objectives.Table(2:end,3);
    IndexLp       = Objectives.IndexLP;
    Nsec          = length(Obj.Nplies);
    Obj.ImportanceFactor = cell2mat(Objectives.Table(2:end,4));
    
    if Constraints.Sym   && Constraints.Balanced,
        LamType = 'Balanced_Sym';
        NpliesUpBound   = round(Obj.Nplies*Constraints.NRange/4)*4;
        NpliesLowBound  = round(Obj.Nplies/4)*4;
        Nmax = max(NpliesUpBound)/4 ;
        Nmin = min(NpliesLowBound)/4;
        for i=1:length(NpliesUpBound)
            AllowedNplies{i} = (NpliesLowBound(i):4:NpliesUpBound(i))/4;
        end
    end
    if Constraints.Sym   && ~Constraints.Balanced,
        LamType = 'Sym';
        NpliesUpBound  = round(Obj.Nplies*Constraints.NRange/2)*2;
        NpliesLowBound  = round(Obj.Nplies/2)*2;
        Nmax = max(NpliesUpBound)/2 ;
        Nmin = min(NpliesLowBound)/2;
        for i=1:length(NpliesUpBound)
            AllowedNplies{i} = (NpliesLowBound(i):2:NpliesUpBound(i))/2;
        end
    end
    if ~Constraints.Sym  && Constraints.Balanced,
        LamType = 'Balanced';
        NpliesUpBound   = round(Obj.Nplies*Constraints.NRange/2)*2;
        NpliesLowBound  = round(Obj.Nplies/2)*2;
        Nmax = max(NpliesUpBound)/2 ;
        Nmin = min(NpliesLowBound)/2;
        for i=1:length(NpliesUpBound)
            AllowedNplies{i} = (NpliesLowBound(i):2:NpliesUpBound(i))/2;
        end
    end
    if  ~Constraints.Sym  && ~Constraints.Balanced,
        LamType = 'Generic';
        NpliesUpBound   = round(Obj.Nplies*Constraints.NRange);
        NpliesLowBound  = round(Obj.Nplies);
        Nmax = max(NpliesUpBound) ;
        Nmin = min(NpliesLowBound);
        for i=1:length(NpliesUpBound)
            AllowedNplies{i} = NpliesLowBound(i):NpliesUpBound(i);
        end
    end

    
    LP_obj = zeros(12,Nsec);
    for i = 1 : Nsec
        LP_obj(:,i) = cell2mat(Obj.Lps2Match(i));
    end
    
    [SortedObj.Nplies,sortIndex] = sort(Obj.Nplies,'descend');
    SortedObj.LP_obj             = LP_obj(:,sortIndex);
    SortedObj.ImportanceFactor   = Obj.ImportanceFactor(sortIndex);
    
    
end


if Constraints.NRange == 1 
    Nsec = 0;
end

% ---
if 1    % Set GA options
    options  = gaoptimset('PopulationSize' ,GAoptions.Npop,'Generation',GAoptions.Ngen, ...
        'StallGenLimit',GAoptions.NgenMin,'EliteCount',ceil(GAoptions.Elitism*GAoptions.Npop),...
        'FitnessLimit' ,1e-5,'TolFun' ,1e-10);
    if GAoptions.Plot
        options  = gaoptimset(options,'PlotFcns' ,{@gaplotbestf});
    end
    
    ConstraintVector = Constraints.Vector;                                  % [Damtol  Rule10percent  Disorientation  Contiguity  DiscreteAngle  InernalContinuity  Covering];
    DeltaAngle       = Constraints.DeltaAngle;
    
    Nvar       = Nsec + Nmax +(Nmax-Nmin) ; % Nsec (thickness) + Nmax Guide plies + Nmax potential drops 

    fct_handle = @(x)FitnessEval(x,Constraints,IndexLp,SortedObj,Nsec,AllowedNplies,LamType);  % handle of the fitness function becomes constraint
    
    % ---
    if 1  % design variable boundaries
        if ConstraintVector(5) % if DiscreteAngle
            IntegerDV = 1:Nvar;                                             % discrete
            Nd_state  = length(-90:DeltaAngle:90);                           % number of discrete state
            LB = [0*ones(Nmax,1)            ; ones(Nmax-Nmin,1)];
            UB = [(Nd_state-1)*ones(Nmax,1) ; Nmax*ones(Nmax-Nmin,1)];
        else
            keyboard
            IntegerDV = [[1:Nsec] [(Nsec+Nmax):(Nsec+2*Nmax)]];             % continuous
            LB = [-90*ones(Nmax,1); ones(Nmax-Nmin,1)];
            UB = [90*ones(Nmax,1);  Nmax*ones(Nmax-Nmin,1)];
        end
        for iply=1:Nsec
            LB = [AllowedNplies{iply}(1)   ;LB];
            UB = [AllowedNplies{iply}(end) ;UB];
        end
        
        if ConstraintVector(1) % if Damtol, make the first ply +- 45
            if ~ConstraintVector(5),
                IntegerDV = [Nsec+1 IntegerDV];
            end
            LB(Nsec+1) = 1;
            UB(Nsec+1) = 2;
        end
    end

end


% --- Generate Ini. Pop.
for i = 1:5
    try
        [IniPop] = CreateBlendedIniPop (Nvar,GAoptions.Npop,Nmax,Nmin,Constraints,AllowedNplies,LamType);
        break; 
    catch
        fprintf('Retrying IniPop Gen. ...\n');
        if i == 5
             UserResponse = input('Non-feasible Elmts in the IniPop, Do you want to continue (y|n)','s');
             if strcmp(UserResponse,'y')
                 temp = Constraints;
                 temp.Vector = 0* temp.Vector;
                 [IniPop] = SST_CreateBlendedIniPop (Nvar,GAoptions.Npop,Nsec,Nmax,Nmin,temp,LamType);
             else 
                 error('Code stopped by the user')
             end
        end
    end
    
end
options = gaoptimset(options,'InitialPopulation' ,IniPop);


% --- run GA
fprintf(strcat('Running GA \n'))
[xOpt,fval] = ga(fct_handle,Nvar,[],[],[],[],LB ,UB,[],IntegerDV,options);
display('GA(s) Terminated Successfully')

% --- Results
[~,LPMatched,output,fitRMS] = fct_handle(xOpt);
SS    = output.SS;
Table = [{'Lam #'} {'Nplies Ori'} {'Nplies SST'} {'Ply Angles'} {'Lam. Param.'} {'Error %'} {'Error Norm'} {'Error RMS'}];
for j = 1:length(sortIndex)
    LP2Match   = SortedObj.LP_obj(:,j);
    QualIndex1 = 100*sum(abs(  (LPMatched(IndexLp,j) - LP2Match(IndexLp))./LP2Match(IndexLp) ));
    QualIndex2 = norm(LPMatched(IndexLp,j) - LP2Match(IndexLp));
    QualIndex3 = rms(LPMatched(IndexLp,j) - LP2Match(IndexLp));
    Table      = [Table ;  {sortIndex(j)} {SortedObj.Nplies(j)} {length(SS{j})} SS(j) {LPMatched(:,j)} {QualIndex1} {QualIndex2} {QualIndex3}];
end

output.Table     = Table;
output.xOpt      = xOpt;
output.fval      = fval;
output.fitRMS    = fitRMS;

if ~output.FEASIBLE
    warning('The SST GA could not find a single feasible solution!')
end
end
