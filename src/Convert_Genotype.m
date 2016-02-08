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

function [SortedLamNumber,GuideAngles,ShuffleLoc,DropIndexes] = Convert_Genotype(Individual,LamNumber,Constraints,NpatchVar,NthetaVar,NbalVar,N10percentVar,AllowedNplies)

keyboard

%% Extract Number of ply per patches
IndexPly = 1;
NpliesperLam = nan*ones(length(NpatchVar),1);
for iPly = 1:length(NpatchVar)
    if NpatchVar(iPly)==0
        NpliesperLam(iPly) = AllowedNplies{iPly};
    else
        NpliesperLam(iPly) = Individual(IndexPly);
        IndexPly = IndexPly+1;
    end
end



%% sort ply order to have Guide First 
if Constraints.ORDERED
    NpliesperLam = sort(NpliesperLam,'descend');                            % repair to ensure thickness is ordered
    SortIndex    = 1:length(NpliesperLam);
    if find(diff(NpliesperLam)>0,1)
        keyboard % should never happen
    end
else
    [NpliesperLam,SortIndex] = sort(NpliesperLam,'descend');
end

SortedLamNumber = LamNumber(SortIndex);                                        % after 2nd sorting
NGuidePlies     = max(NpliesperLam);                                           % number of plies in the guide laminate (take half for Sym.)
NDropPlies      = abs(diff(NpliesperLam));                                     % number of ply drops between each laminates
GuideAngles     = Individual(sum(NpatchVar) + [1:NGuidePlies]);                % Extract variable fibre angles of the guide


[NthetaVar0,NbalVar0,N10percentVar0] = AttributeDesignVariable(NGuidePlies*16,Constraints);

% if Constraints.Vector(2) == 1 % add the 10% rule framework
%     keyboard
%     Fixed
% end


%% --- Shuffle Location (i.e. location of angles pairs) for balanced Lam.
if Constraints.Balanced
    ShuffleLoc = Individual(sum(NpatchVar) + NthetaVar + [1:NGuidePlies]);
%     StartIndex = sum(NpatchVar) + NthetaVar*2;
else
    ShuffleLoc = [];
%     StartIndex = sum(NpatchVar) + NthetaVar;
end

%% --- Extract 10% rule data

if Constraints(2)
    
   
    TenPercentData = Individual(sum(NpatchVar) + NthetaVar + NbalVar + [1:NGuidePlies]);
else
    TenPercentData = [];
end


%% --- organise ply drop sequences
Ndrop = length(NDropPlies);
DropIndexes = cell(1,Ndrop);                                                   % matrix of drops plies bundled together

for iDrop = 1 : Ndrop
    DropIndexesTEMP     = Individual(StartIndex + [1:NDropPlies(iDrop)]);
    DropIndexes{iDrop}  = DropIndexesTEMP(~isnan(DropIndexesTEMP));
    StartIndex          = StartIndex + NDropPlies(iDrop);
end

end