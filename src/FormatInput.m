% =====                                                              ==== %
%                Used to format input structures for the GA               %
%                       and output usefull quantities
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


function [NStruct,LamType,LB,UB,BCs,AllowedNplies]  = FormatInput(Objectives,Constraints)

%% NthetaVar and LamType

if Constraints.Sym  && Constraints.Balanced,    LamType = 'Balanced_Sym';     end
if Constraints.Sym && ~Constraints.Balanced,    LamType = 'Sym';              end
if ~Constraints.Sym && Constraints.Balanced,    LamType = 'Balanced';         end
if ~Constraints.Sym && ~Constraints.Balanced,   LamType = 'Generic';          end


Nplies = round(cell2mat(Objectives.Table(2:end,2)));    % Upper and lower Number of plies allowed for each patch

NStruct =  AttributeNply(Nplies,Constraints,LamType);

%% NpatchVar and AllowedNplies
AllowedNplies = cell(size(Nplies,1),1);            % Number of Ply Range allowed
for i=1:size(Nplies,1)
    if strcmp(LamType,'Balanced_Sym')
        AllowedNplies{i} = Nplies(i,1):2:Nplies(i,2);      
    else
        AllowedNplies{i} = Nplies(i,1):Nplies(i,2);       
    end
end


Nrange = cellfun(@max,AllowedNplies,'UniformOutput', true) - cellfun(@min,AllowedNplies,'UniformOutput', true); % max ply - min ply per lam.
NStruct.NpatchVar = boolean(zeros(length(Nrange),1));
NStruct.NpatchVar(Nrange>0)=1; % number of variable thickness lam.



%% Genotype Upper and Lower Bounds (LB and UB)

%[Nply, Theta's, -Theta's Location , 10% rule location, MidPlaneAngle , ply drops ]



NStruct.Nvar = sum(NStruct.NpatchVar) + NStruct.NDvs + NStruct.NdropVar;                                   % total number of design variables
Nd_state     = length(-90:Constraints.DeltaAngle:90);                              % number of discrete state for variable angles

LB = [];
UB = [];

% Bounds for Variable number of plies
LB = [LB; ones(sum(NStruct.NpatchVar),1)];  
UB = [UB; cellfun(@length,AllowedNplies(NStruct.NpatchVar),'UniformOutput', true)];

BCs.LB.Nply = ones(sum(NStruct.NpatchVar),1);
BCs.UB.Nply = cellfun(@length,AllowedNplies(NStruct.NpatchVar),'UniformOutput', true);


 % Bounds for theta's
LB = [LB; ones(NStruct.NthetaVar,1)];  
UB = [UB; (Nd_state-1)*ones(NStruct.NthetaVar,1)];

BCs.LB.Thetas = ones(NStruct.NthetaVar,1);
BCs.UB.Thetas = (Nd_state-1)*ones(NStruct.NthetaVar,1);


% Bounds for Balanced angles and 10% rule (first and last ply are Theta's)
LB = [LB; 2*ones(NStruct.NbalVar + NStruct.N10percentVar,1)];                                             
UB = [UB; (NStruct.NDvs-1)*ones(NStruct.NbalVar + NStruct.N10percentVar,1)];                               

BCs.LB.BalancedLoc = 2*ones(NStruct.NbalVar + NStruct.N10percentVar,1);
BCs.UB.BalancedLoc = (NStruct.NDvs-1)*ones(NStruct.NbalVar + NStruct.N10percentVar,1);


% Bounds for Midplane Angles (max of 2) 
% (either 0 or 90)

NStruct.NDV_NMidPlane = 2;
if Constraints.Balanced
    LB = [LB; ones(2,1)];
    UB = [UB; 2*ones(2,1)];

    BCs.LB.MidPlane = ones(2,1); 
    BCs.UB.MidPlane = 2*ones(2,1); 
else
    LB = [LB; ones(2,1)];
    UB = [UB; (Nd_state-1)*ones(2,1)];

    BCs.LB.MidPlane = ones(2,1); 
    BCs.UB.MidPlane = (Nd_state-1)*ones(2,1); 
end


% Bounds for Drop ply locations
if Constraints.Vector(7) % covering
    LB = [LB; 2*ones(NStruct.NdropVar,1)];                                              % remove the first
    BCs.LB.PlyDrop = 2*ones(NStruct.NdropVar,1);
    
    if strcmp(LamType,'Generic')
       UB = [UB; (NStruct.N10percentVar + NStruct.NthetaVar + NStruct.NDV_NMidPlane-1)*ones(NStruct.NdropVar,1)];
       BCs.UB.PlyDrop = (NStruct.N10percentVar + NStruct.NthetaVar +NStruct.NDV_NMidPlane -1)*ones(NStruct.NdropVar,1);
       
    else
       UB = [UB; (NStruct.N10percentVar + NStruct.NthetaVar + NStruct.NDV_NMidPlane)*ones(NStruct.NdropVar,1)]; 
       BCs.UB.PlyDrop = (NStruct.N10percentVar + NStruct.NthetaVar +NStruct.NDV_NMidPlane)*ones(NStruct.NdropVar,1);
    end
else
    LB = [LB; ones(NStruct.NdropVar,1)];
    BCs.LB.PlyDrop = ones(NStruct.NdropVar,1); 
    
    UB = [UB; (NStruct.N10percentVar + NStruct.NthetaVar + NStruct.NDV_NMidPlane)*ones(NStruct.NdropVar,1)];
    BCs.UB.PlyDrop = (NStruct.N10percentVar + NStruct.NthetaVar + NStruct.NDV_NMidPlane)*ones(NStruct.NdropVar,1);
end

if Constraints.Vector(1) % if Damtol, make the first ply +- 45
    LB(sum(NStruct.NpatchVar)+1) = 1;
    UB(sum(NStruct.NpatchVar)+1) = 2;
    
    BCs.LB.Thetas(1) = 1; 
    BCs.UB.Thetas(1) = 2; 
end
end