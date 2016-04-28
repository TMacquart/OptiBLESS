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
%           Convert the coded genotype into a Stacking sequence table, 
%                       and the number of plies per patch 
%
% [NpliesPerLam,SSTable,NplySS] = Convert_Genotype(Individual,Constraints,NStruct,AllowedNplies,LamType)
% =====                                                              ====== 

function [NpliesPerLam,SSTable,NplySS] = Convert_Genotype(Individual,Constraints,NStruct,NStructMin,AllowedNplies,LamType)

% Separate design variables from the individual
NVarPatch       = sum(NStruct.NpatchVar); % number of patches with variable thickness 
Thetas          = Individual(NVarPatch + [1:NStruct.NthetaVar]);                
BalancedLoc     = Individual(NVarPatch + NStruct.NthetaVar + [1:NStruct.NbalVar]);
TenPercentLoc   = Individual(NVarPatch + NStruct.NthetaVar + NStruct.NbalVar + [1:NStruct.N10percentVar]);

if NStruct.NMidPlane>=NStruct.NDV_NMidPlane
    % only happens for Sym. Balanced (NStruct.NMidPlane can be = 3 at max), in this case need to ensure symmetry is preserved
    Thetas_Mid      = Individual(NVarPatch + NStruct.NthetaVar + NStruct.NbalVar + NStruct.N10percentVar + [1:NStruct.NDV_NMidPlane]);
else
    % NStruct.NMidPlane = 1
    Thetas_Mid      = Individual(NVarPatch + NStruct.NthetaVar + NStruct.NbalVar + NStruct.N10percentVar + [1:NStruct.NMidPlane]);
end

InsertIndexes  = Individual(NVarPatch + NStruct.NthetaVar + NStruct.NbalVar + NStruct.N10percentVar + NStruct.NDV_NMidPlane + [1:NStruct.NInsertVar]);



% compute stacking sequence table
[SSTable,NplySS] = ComputeSSTable(Thetas,InsertIndexes,BalancedLoc,TenPercentLoc,Thetas_Mid,LamType,Constraints,NStruct,NStructMin,false);   



% --- Extract Number of ply per patches and repair to match SS_Table options
NpliesPerLam = nan*ones(length(NStruct.NpatchVar),1);
for iPly = 1:length(NStruct.NpatchVar)
    if NStruct.NpatchVar(iPly)==0
        % if not a design variable (only 1 choice)
        NpliesPerLam(iPly) = AllowedNplies{iPly};
        
    else
        % if a design variable de-code the information
        NpliesPerLam(iPly) = AllowedNplies{iPly}(Individual(iPly));
        
    end
end


% replace Nply by the closest one available
for j = 1:length(NpliesPerLam(:,1))
    if isempty(find(NpliesPerLam(j,1)==NplySS,1))
        [~,NIndex] = min(abs(NplySS-NpliesPerLam(j,1) ));
        NpliesPerLam(j,1) = NplySS(NIndex);
    end
end

end
