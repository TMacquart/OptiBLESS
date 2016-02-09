% =====                                                              ==== 
%               Convert the coded genotype into a guide laminate, 
%           balanced angles shuffling locations and ply drop locations
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

function [NpliesPerLam,SSTable] = Convert_Genotype(Individual,Constraints,NStruct,AllowedNplies,LamType)


% Extract Number of ply per patches
NpliesPerLam = nan*ones(length(NStruct.NpatchVar),1);
for iPly = 1:length(NStruct.NpatchVar)
    if NStruct.NpatchVar(iPly)==0
        NpliesPerLam(iPly) = AllowedNplies{iPly};
    else
        NpliesPerLam(iPly) = AllowedNplies{iPly}(Individual(iPly));
    end
end

NpliesPerLam = [NpliesPerLam [1:length(NpliesPerLam)]'];
NpliesPerLam = flipud(sortrows(NpliesPerLam,1));

NStruct_0 = AttributeNply(NpliesPerLam(:,1),Constraints,LamType);


% Extract All design variable from the individual
Thetas          = Individual(sum(NStruct.NpatchVar) + [1:NStruct_0.NthetaVar]);                
BalancedLoc     = Individual(sum(NStruct.NpatchVar) + NStruct.NthetaVar + [1:NStruct_0.NbalVar]);
TenPercentLoc   = Individual(sum(NStruct.NpatchVar) + NStruct.NthetaVar + NStruct.NbalVar + [1:NStruct_0.N10percentVar]);
Thetas_Mid      = Individual(sum(NStruct.NpatchVar) + NStruct.NthetaVar + NStruct.NbalVar + NStruct.N10percentVar + [1:NStruct_0.NMidPlane]);
PlyDrops        = Individual(sum(NStruct.NpatchVar) + NStruct.NthetaVar + NStruct.NbalVar + NStruct.N10percentVar + NStruct.NDV_NMidPlane + [1:NStruct_0.NdropVar]);


 Individual0 = [NpliesPerLam(:,1); 
                          Thetas;
                          BalancedLoc;
                          TenPercentLoc;
                          Thetas_Mid;
                          PlyDrops ];
                      
% compute SSTAble
SSTable = NewComputeSSTable(Thetas,PlyDrops,BalancedLoc,TenPercentLoc,Thetas_Mid,LamType,Constraints);   
end