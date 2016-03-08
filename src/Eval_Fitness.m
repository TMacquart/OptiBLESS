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
%
%
% =====                                                              ====== 
%                 Fitness evaluation function used by the GA              
%                                                                                                        
% [fitness,output] = Eval_Fitness (Individual,Objectives,Constraints,NStruct,AllowedNplies,LamType)
%
% =====                                                              ====== 

function [fitness,output] = Eval_Fitness (Individual,Objectives,Constraints,NStruct,NStructMin,AllowedNplies,LamType)

% Extract laminate numbers
Nlam      = size(Objectives.Table,1)-1;
LamNumber = int8(zeros(Nlam,1));
for j=1:Nlam
    LamNumber(j) = Objectives.Table{j+1,1}; 
end


% Convert individual genotype into stacking sequence table

[NpliesPerLam,SSTable,NplySS] = Convert_Genotype(Individual,Constraints,NStruct,NStructMin,AllowedNplies,LamType);

[NUniquePliesPerLam] = unique(NpliesPerLam(:,1));
NUniquePliesPerLam   = flipud(NUniquePliesPerLam);
NUnique              = length(NUniquePliesPerLam);

FEASIBLE = Check_Feasibility(Constraints,SSTable,NplySS);



%% Compute Fitness

% --- pre-allocation
    SS_Unique = cell(size(SSTable));
    SS_Patch  = cell(Nlam,size(SSTable,2));

    if strcmp(Objectives.Type,'LP')
        LP = zeros(12,NUnique);
        LP_Patch = zeros(12,Nlam);
    end

    if strcmp(Objectives.Type,'ABD')
        A = cell(NUnique,1);
        B = cell(NUnique,1);
        D = cell(NUnique,1);

        A_Patch = cell(NUnique,1);
        B_Patch = cell(NUnique,1);
        D_Patch = cell(NUnique,1);
    end
% ---


% Only keep unique stacking sequences 
for iSS = 1 : NUnique
    RowIndex = find(NplySS==NUniquePliesPerLam(iSS),1);
    
    SS_Unique(iSS,:)   = SSTable(RowIndex,:);    
    NplySS_Unique(iSS) = length(find(~cellfun(@isempty,SS_Unique(iSS,:)))); %#ok<AGROW>
    FiberAngles        = cell2mat(SS_Unique(iSS,:));
    
    if strcmp(Objectives.Type,'LP')
        LP(:,iSS) = Convert_SS2LP(FiberAngles);          % evaluate lamination parameters for the droped laminates
    end
    if strcmp(Objectives.Type,'ABD')
        [A{iSS},B{iSS},D{iSS}] = Convert_SS2ABD (Objectives.mat(1),Objectives.mat(2),Objectives.mat(4),Objectives.mat(3),Constraints.ply_t,FiberAngles,true);
    end
end


% Assign each patch with its corresponding stacking sequence, LP and ABD
for iPatch = 1:Nlam
    RowIndex = find(NplySS_Unique == NpliesPerLam(iPatch),1);
    SS_Patch(iPatch,:) = SS_Unique(RowIndex,:);
    
    if strcmp(Objectives.Type,'LP')
        LP_Patch(:,iPatch) = LP(:,RowIndex);
    end

    if strcmp(Objectives.Type,'ABD')
        A_Patch(iPatch) = A(RowIndex);
        B_Patch(iPatch) = B(RowIndex);
        D_Patch(iPatch) = D(RowIndex);
    end
end

% Remove Empty Columns (ply not used)
Temp  = sum(cellfun(@isempty,SS_Patch));
EmptyCol_Index = find(Nlam==Temp);
if ~isempty(EmptyCol_Index)
    SS_Patch(:,EmptyCol_Index) = [];
end



% Compute Lamination parameter based fitness
if strcmp(Objectives.Type,'LP')
    if ~Objectives.UserFct
        fitness          = Objectives.FitnessFct(LP_Patch);                       % Default Fitness Function (Do not Change)
    else
        [fitness,output] = Objectives.FitnessFct(LP_Patch);                       % User Fitness Function Calls
    end
    output.LP = LP_Patch;
end


% Compute Stiffness based fitness
if strcmp(Objectives.Type,'ABD')
    if ~Objectives.UserFct
        fitness = Objectives.FitnessFct(A_Patch,B_Patch,D_Patch);                             % Default Fitness Function (Do not Change)
    else
        [fitness,output] = Objectives.FitnessFct(A_Patch,B_Patch,D_Patch);                    % User Fitness Function Calls
    end 
    output.A  = A_Patch;
    output.B  = B_Patch;
    output.D  = D_Patch;
end


% Compute Stacking sequence based fitness
if strcmp(Objectives.Type,'SS')
    [fitness,output] = Objectives.FitnessFct(SS_Patch);                           % User Fitness Function Calls
end



% --- check individual ply continuity (only if structure geometry is given)
    if isfield(Constraints,'PatchConnectivity') 
        NGeoConstraints = CheckContinuity(SS_Patch,Constraints.PatchConnectivity);
    else
        NGeoConstraints = 0;
    end
% ---

% keyboard
% Penalise infeasible solutions
MaxNGeoConstraints = 50;            % abritraty for now 
fitness = fitness * (1 + NGeoConstraints/MaxNGeoConstraints);

if ~FEASIBLE  
    if isnan(fitness) || isinf(fitness) || ~isreal(fitness) || size(fitness,1)~=1 ||  size(fitness,2)~=1 
        error('Non appropriate Fitness (i.e. Non-Scalar, NaN, inf or complex) has been detected')
    end
    fitness = fitness * 2;
end


% return
output.SS_Patch = SS_Patch;
output.FEASIBLE = FEASIBLE;
output.NGeoConstraints = NGeoConstraints;


