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
%                       Terence Macquart (09/11/2015)                     %
%                                                                         %
%
%  The individual of GA is composed of 2 parts:
%  --------------------------------------------
%  [ Theta(1) ... Theta(N)]  -- the guide laminate ply angles
%  [ Drop(1)  ... Drop(M) ]  -- the guide ply drops 
% =====                                                              ==== %

function [output] = RetrieveSS_MatchLP(Objectives,Constraints,GAoptions)


if 1 % Read and Format objectives (i.e. LPs to match) -- do not change

    Obj.LamInd       = cell2mat(Objectives.Table(2:end,1));
    Obj.Nplies       = cell2mat(Objectives.Table(2:end,2));
    Obj.Lps2Match    = Objectives.Table(2:end,3);
    IndexLp          = Objectives.IndexLP; 
    ImportanceFactor = cell2mat(Objectives.Table(2:end,4));
    
    % --- Re-organise Nplies and LPs to match from thickest to thinnest lam. 
    [NpliesUnique,UniqueIndex] = unique(Obj.Nplies);
    NonUniqueIndex = find(~ismember([1:length(Obj.Nplies)],UniqueIndex));

    NplyVector        = cell(length(NpliesUnique),3);
    NplyVector(:,1:2) = num2cell([NpliesUnique UniqueIndex]);
    for i = 1 : length(NonUniqueIndex)
        irow = find(Obj.Nplies(NonUniqueIndex(i))==NpliesUnique,1);
        NplyVector(irow,3) = {[cell2mat(NplyVector(irow,3)) NonUniqueIndex(i)]};
    end

    NpliesSortedArray = sortrows(NplyVector,-1);
    clear NplyVector NpliesUnique UniqueIndex

    Nplies     = (cell2mat(NpliesSortedArray(:,1)));
    NDropPlies = (cell2mat(NpliesSortedArray(1:end-1,1)) - cell2mat(NpliesSortedArray(2:end,1)));      
    
    LP_obj = zeros(12,size(NpliesSortedArray,1));
    for i = 1 : size(NpliesSortedArray,1)
        if isempty(NpliesSortedArray{i,3})
            LP_obj(:,i) = cell2mat(Obj.Lps2Match(NpliesSortedArray{i,2}));
        else
            % if multiple laminates with the same number of plies take the mean of LPs 
            LPs = cell2mat(Obj.Lps2Match(NpliesSortedArray{i,2}));
            for j = 1 : (length(NpliesSortedArray{i,3}))
                LPs = LPs + cell2mat( Obj.Lps2Match(NpliesSortedArray{i,3}(j)) );
            end
            LP_obj(:,i) = LPs/(length(NpliesSortedArray{i,3})+1);
        end
    end
end




% --- 
if Constraints.Sym   && Constraints.Balanced,     LamType = 'Balanced_Sym';  NguidePlies = max(Nplies)/4;  NguideDrops = NDropPlies/4;   end  
if Constraints.Sym   && ~Constraints.Balanced,    LamType = 'Sym';           NguidePlies = max(Nplies)/2;  NguideDrops = NDropPlies/2;   end  
if ~Constraints.Sym  && Constraints.Balanced,     LamType = 'Balanced';      NguidePlies = max(Nplies)/2;  NguideDrops = NDropPlies/2;   end  
if ~Constraints.Sym  && ~Constraints.Balanced,    LamType = 'Generic';       NguidePlies = max(Nplies);    NguideDrops = NDropPlies;     end 
Nvar = NguidePlies + sum(NguideDrops);



% ---
if 1  % Set GA options
    options  = gaoptimset('PopulationSize' ,GAoptions.Npop,'Generation',GAoptions.Ngen, ...
        'StallGenLimit',GAoptions.NgenMin,'EliteCount',ceil(GAoptions.Elitism*GAoptions.Npop),...
        'FitnessLimit' ,1e-5,'TolFun' ,1e-10);
    if GAoptions.Plot
        options  = gaoptimset(options,'PlotFcns' ,{@gaplotbestf});
    end
    ConstraintVector = Constraints.Vector;
    DeltaAngle       = Constraints.DeltaAngle;
    
    fct_handle = @(x)LPMatch_FitnessEval(x,Constraints,NguidePlies,NguideDrops,IndexLp,LP_obj,ImportanceFactor,LamType);   % handle of the fitness function
                                        
    if ConstraintVector(5) % if DiscreteAngle
        IntegerDV = 1:Nvar;                                                              % set all dv discrete
        Nd_state  = length(-90:DeltaAngle:90);                                           % number of discrete possible states
        LB = [0*ones(NguidePlies,1)            ; ones(sum(NguideDrops),1)];              % dv lower bound (discrete)
        UB = [(Nd_state-1)*ones(NguidePlies,1) ; NguidePlies*ones(sum(NguideDrops),1)];  % dv upper bound (discrete)
    else
        IntegerDV = NguidePlies + [1:sum(NguideDrops)];                                 % only ply drops are discrete
        LB = [-90*ones(NguidePlies,1) ; ones(sum(NguideDrops),1)];
        UB = [+90*ones(NguidePlies,1) ; NguidePlies*ones(sum(NguideDrops),1)];
    end
    
    if ConstraintVector(1) % if Damtol, make the first ply a integer dv. (+- 45)
        if ~ConstraintVector(5),
            IntegerDV = [1 IntegerDV];
        end
        LB(1) = 1; 
        UB(1) = 2;
    end
    
end


% --- Generate Initial Population
for i = 1:5
    try
        [IniPop] = LPMatch_CreateBlendedIniPop (GAoptions.Npop,NguidePlies,NguideDrops,Constraints,LamType);
        break;  % Break out of the i-loop on success
    catch
        fprintf('Retrying IniPop Gen. ...\n');
    end
end
options = gaoptimset(options,'InitialPopulation' ,IniPop);


% --- run GA
display('Running GA')
[xOpt,fval] = ga(fct_handle,Nvar,[],[],[],[],LB,UB,[],IntegerDV,options);
display('GA(s) Terminated Successfully')


% --- Results
[~,LPMatched,output] = fct_handle(xOpt);
DropsIndexes = output.DropsIndexes;
SS           = output.SS;

Table = [{'Lam #'} {'Nplies'} {'Ply Angles'} {'Lam. Param.'} {'Error %'} {'Error Norm'}];
for iDrop = 1 : length(DropsIndexes)+1
    LP2Match = LP_obj(:,iDrop); % sorted objectives
    
    QualityMatch = QualityEval(LP2Match,LPMatched(:,iDrop),IndexLp);
    Table        = [Table ;  NpliesSortedArray(iDrop,2) NpliesSortedArray(iDrop,1) SS(iDrop) {LPMatched(:,iDrop)} {QualityMatch(1)} {QualityMatch(2)}];
    
    if ~isempty(NpliesSortedArray{iDrop,3})
        for j = 1:length(NpliesSortedArray{iDrop,3})
            LP2Match     = cell2mat(Obj.Lps2Match(NpliesSortedArray{iDrop,3}(j)));   
            QualityMatch = QualityEval(LP2Match,LPMatched(:,iDrop),IndexLp);
            Table        = [Table ;  NpliesSortedArray{iDrop,3}(j) NpliesSortedArray(iDrop,1) SS(iDrop)  {LPMatched(:,iDrop)} {QualityMatch(1)} {QualityMatch(2)}];
        end
    end
end



output.Table = Table;
output.xOpt  = xOpt;
output.fval  = fval;
if ~output.FEASIBLE
    warning('The GA could not find a single feasible solution!')
end

end

function QualityMatch = QualityEval(LP2Match,LPMatched,IndexLp) 
    QualityMatch(1) = 100*sum(abs( (LPMatched(IndexLp) - LP2Match(IndexLp))./LP2Match(IndexLp) )); % sum of the absolute percentage of matching error
    QualityMatch(2) = norm(LP2Match(IndexLp) - LPMatched(IndexLp));                                % norm of the matching error
end
