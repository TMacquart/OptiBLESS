% =====                                                              ==== %
%                Used to format input structures for the GA               %
%
% [Nvar,NpatchVar,NthetaVar,NdropVar,LamType,LB,UB,AllowedNplies]  = FormatInput(Objectives,Constraints)
%
% 
% Nvar          - Number of design variables 
% NpatchVar     - Number of patches with variable number of plies
% NthetaVar     - Number of fibre angles used to describe the guide laminate
% NdropVar      - Number of drops from the guide laminate
% LamType       - Type of laminates
% LB            - Lower bound of design variables
% UB            - Upper bound of design variables
% AllowedNplies - Number of plies allowed for each patches
% =====                                                              ==== %

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


function [Nvar,NpatchVar,NthetaVar,NdropVar,LamType,LB,UB,AllowedNplies]  = FormatInput(Objectives,Constraints)

%% NthetaVar and LamType

NplyIni = cell2mat(Objectives.Table(2:end,2));

if Constraints.Sym  && Constraints.Balanced,
    DeltaNPly = 4;                  % The allowed ply numbers are multiple of DeltaNPly
    LamType   = 'Balanced_Sym';     % Type of Laminate
end

if Constraints.Sym && ~Constraints.Balanced,
    DeltaNPly = 2;
    LamType   = 'Sym';
end

if ~Constraints.Sym && Constraints.Balanced,
    DeltaNPly = 2;
    LamType   = 'Balanced';
end

if  ~Constraints.Sym && ~Constraints.Balanced,
    LamType   = 'Generic';
    DeltaNPly = 1;
end

Nplies    = round(NplyIni/DeltaNPly)*DeltaNPly;                    % Upper and lower Number of plies allowed for each patch
NthetaVar = max(Nplies(:))/DeltaNPly;                              % Max number of variable angles

%% NdropVar
[MaxNplies,rowIndMax] = max(Nplies(:,2)); 
NpliesTemp = Nplies;
NpliesTemp(rowIndMax,:) = [];
if ~isempty(NpliesTemp), 
    NdropVar = (MaxNplies-min(NpliesTemp(:)))/DeltaNPly;           % Max number of guide drops variable
    clear NpliesTemp
else
    NdropVar = 0;
end

%%
if Constraints.Balanced
    NbalVar = NthetaVar;                                           % Max number of balanced angles
else
    NbalVar = 0;
end

%% NpatchVar and AllowedNplies
AllowedNplies = cell(size(NplyIni,1),1);
for i=1:size(NplyIni,1)
    AllowedNplies{i} = (Nplies(i,1):DeltaNPly:Nplies(i,2))/DeltaNPly;       % Number of Ply Range
end


Nrange = cellfun(@max,AllowedNplies,'UniformOutput', true) - cellfun(@min,AllowedNplies,'UniformOutput', true); % max ply - min ply per lam.
NpatchVar = zeros(length(Nrange),1);
NpatchVar(Nrange>0)=1; % number of variable thickness lam.



%% Genotype Upper and Lower Bounds (LB and UB)

Nvar      = sum(NpatchVar) + NthetaVar + NbalVar + NdropVar;
Nd_state  = length(-90:Constraints.DeltaAngle:90);              % number of discrete state for variable angles

LB = [];
UB = [];

LB = [LB; cellfun(@min,AllowedNplies(find(NpatchVar)),'UniformOutput', true)];  % Variable plies
UB = [UB; cellfun(@max,AllowedNplies(find(NpatchVar)),'UniformOutput', true)];

LB = [LB; 0*ones(NthetaVar,1)];                                                 % Guide angles
UB = [UB; (Nd_state-1)*ones(NthetaVar,1)];

if ~Constraints.Vector(1) % Damtol
    LB = [LB; ones(NbalVar,1)];                                                 % Balanced angles
    UB = [UB; 2*NbalVar*ones(NbalVar,1)];
else
    LB = [LB; 2*ones(NbalVar,1)];                                               % Ensure the First ply is from the guide (+- 45)
    if ~strcmp(LamType,'Generic')
        UB = [UB; 2*NbalVar*ones(NbalVar,1)];
    else
        UB = [UB; (2*NbalVar-1)*ones(NbalVar,1)];                               % Ensure the last ply is also from the guide (+- 45)
    end
end
    
if Constraints.Vector(7) % covering
    if strcmp(LamType,'Generic')
        LB = [LB; 2*ones(NdropVar,1)];                  % remove the first
        UB = [UB; (NthetaVar-1)*ones(NdropVar,1)];      % remove the last
    else
        LB = [LB; 2*ones(NdropVar,1)];                  % remove the first
        UB = [UB; NthetaVar*ones(NdropVar,1)];
    end
else
    LB = [LB; ones(NdropVar,1)];
    UB = [UB; NthetaVar*ones(NdropVar,1)];
end

if Constraints.Vector(1) % if Damtol, make the first ply +- 45
    LB(sum(NpatchVar)+1) = 1;
    UB(sum(NpatchVar)+1) = 2;
end
end