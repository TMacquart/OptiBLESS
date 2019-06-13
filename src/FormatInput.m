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
%                         Formats inputs for the GA
%
% - Returns the number of design variable for each genotype fields and
%   design variables bounds.
%
%  [NStruct,LamType,LB,UB,BCs,AllowedNplies] = FormatInput(Objectives,Constraints)
%
%
% NStruct.NthetaVar;        % number of theta's  design variables
% NStruct.NbalVar;          % number of balanced design variables (will always be equal to NthetaVar, but kept for sake of clarity)
% NStruct.N10percentVar;    % number of 10% rule design variables
% NStruct.NMidPlane;        % number of MidPlane design variables
% NStruct.NInsertVar;         % number of Drop Location design variables
%
% LamType       - Type of laminates
% LB            - Lower bound of design variables (GA Format)
% UB            - Upper bound of design variables (GA Format)
% BCs           - Design variable bounds stored in a more readable format
% AllowedNplies - Number of plies allowed for each patch
% =====                                                              ======


function [NStruct,NStructMin,LamType,LB,UB,BCs,AllowedNplies] = FormatInput(Objectives,Constraints)

if Constraints.Vector(1)  && Constraints.Vector(2),    LamType = 'Balanced_Sym';     end
if Constraints.Vector(1)  && ~Constraints.Vector(2),   LamType = 'Sym';              end
if ~Constraints.Vector(1) && Constraints.Vector(2),    LamType = 'Balanced';         end
if ~Constraints.Vector(1) && ~Constraints.Vector(2),   LamType = 'Generic';          end


Nplies  = round(cell2mat(Objectives.Table(2:end,2)));   % Upper and lower Number of plies allowed for each patch (as given in the Objective structure)

NStruct    = Attribute_NDvs(Nplies,Constraints,LamType);   % return the number of each design variable field (structure of the genotype)
NStructMin = Attribute_NDvs([min(Nplies(:)) min(Nplies(:))],Constraints,LamType);   % return the number of each design variable field (structure of the genotype)


%% Allowed number of plies for each patch and identify the ones with variable thickness
AllowedNplies = cell(size(Nplies,1),1);
for i=1:size(Nplies,1)
    AllowedNplies{i} = Nplies(i,1):Nplies(i,2);
end

Nrange = cellfun(@max,AllowedNplies,'UniformOutput', true) - cellfun(@min,AllowedNplies,'UniformOutput', true); % max ply - min ply per lam.
NStruct.NpatchVar = boolean(zeros(length(Nrange),1));
NStruct.NpatchVar(Nrange>0)=1; % variable thickness lam. are set to 1 (others are set to zero and are not considered as design variables)



%% Genotype design variables coded upper and lower bounds (BCs and LB, UB)
% Genotype Fields = [ [Nply], [Theta's], [-Theta's Balanced Location] , [10% rule angle location], [MidPlaneAngle] , [ply drops] ]

NStruct.NDV_NMidPlane = 2;                                                          % forced number of design variable for midplance (sometines unused - i.e. non coding genes)
NStruct.Nvar = sum(NStruct.NpatchVar) + NStruct.NthetaVar + NStruct.NbalVar ...     % total number of design variables
    + NStruct.N10percentVar + NStruct.NInsertVar + NStruct.NDV_NMidPlane ;
Nd_state     = length(-90:Constraints.DeltaAngle:90);                               % number of discrete state for variable theta's angles



% Bounds for Variable number of plies (coded from 1 to Nrange+1)
if ~isempty(sum(NStruct.NpatchVar))
    BCs.LB.Nply = ones(sum(NStruct.NpatchVar),1);
    BCs.UB.Nply = cellfun(@length,AllowedNplies(NStruct.NpatchVar),'UniformOutput', true);
else
    BCs.LB.Nply = [];
    BCs.UB.Nply = [];
end


% Bounds for Theta's (coded from 1 to (Nd_state-1))
BCs.LB.Thetas = ones(NStruct.NthetaVar,1);
BCs.UB.Thetas = (Nd_state-1)*ones(NStruct.NthetaVar,1);     % (Nd_state-1) to avoid repeating -90 and 90


% Bounds for Balanced angles location within the laminate
NDv_Location = NStruct.NthetaVar + NStruct.NbalVar + NStruct.N10percentVar + NStruct.NMidPlane;
if NStruct.NbalVar~=0
    if Constraints.Vector(3) || Constraints.Vector(8)
        % if Damtol or covering are active the first and last plies of
        % theta's cannot be substituted (and later, cannot be dropped)
        BCs.LB.BalancedLoc = 2*ones(NStruct.NbalVar,1);
        
        if Constraints.Vector(1) % symmetry
            BCs.UB.BalancedLoc = NDv_Location     * ones(NStruct.NbalVar,1);
        else
            BCs.UB.BalancedLoc = (NDv_Location-1) * ones(NStruct.NbalVar,1);
        end
        
    else
        % Any ply can be substituted
        BCs.LB.BalancedLoc = ones(NStruct.NbalVar,1);
        BCs.UB.BalancedLoc = NDv_Location*ones(NStruct.NbalVar,1);
    end
else
    BCs.LB.BalancedLoc = [];
    BCs.UB.BalancedLoc = [];
end


% Bounds for 10% rule angle location within the laminate
if NStruct.N10percentVar~=0
    
    if 0 % Option 1
        % 10% design variable are uniformly disitributed throughout the thickness
        % (helps generating feasible solution complying with disorientation and possibily increase robustness)
        LinSpacing = round(linspace(2,NDv_Location,(NStruct.N10percentVar+1)));
        BCs.LB.TenPercentLoc = LinSpacing(1:end-1)' +[0; ones(length(LinSpacing(1:end-2)),1)];
        BCs.UB.TenPercentLoc = LinSpacing(2:end)';
        
    else % Option 2
        BCs.LB.TenPercentLoc = 2*ones(NStruct.N10percentVar,1);
        BCs.UB.TenPercentLoc = NDv_Location*ones(NStruct.N10percentVar,1);
        
    end

    if ~isempty(find( (BCs.UB.TenPercentLoc-BCs.LB.TenPercentLoc)<0,1))
        display([BCs.LB.TenPercentLoc BCs.UB.TenPercentLoc])
        keyboard % internal check, should never happen
    end
else
    BCs.LB.TenPercentLoc = [];
    BCs.UB.TenPercentLoc = [];
end



% Bounds for Midplane Angles (max of 2 angles)
if Constraints.Vector(2) % Balanced
    % either 0 or 90 deg are possible (coded as 1 and 2)
    BCs.LB.MidPlane = ones(NStruct.NDV_NMidPlane,1);
    BCs.UB.MidPlane = 2*ones(NStruct.NDV_NMidPlane,1);
else
    % any angles used for theta's
    BCs.LB.MidPlane = ones(NStruct.NDV_NMidPlane,1);
    BCs.UB.MidPlane = (Nd_state-1)*ones(NStruct.NDV_NMidPlane,1);
end



% keyboard
% Bounds for insert ply locations
if NStruct.NInsertVar~=0
    
    BCs.LB.InsertIndex = 0*ones(NStruct.NInsertVar,1);
    
    if 1
        BCs.UB.InsertIndex = NDv_Location*ones(NStruct.NInsertVar,1);
    else
        LastDV_INdex       = (NStruct.NthetaVar + NStruct.N10percentVar + NStruct.NMidPlane);
        LinSpaceIndex      = round(linspace(2,LastDV_INdex,NStruct.NInsertVar));
        
        for i=1:NStruct.NInsertVar
            if i ==1
                BCs.UB.InsertIndex(i,1) = LinSpaceIndex(i+2);
            elseif i == NStruct.NInsertVar
                BCs.UB.InsertIndex(i,1) = LinSpaceIndex(i);
            else
                BCs.UB.InsertIndex(i,1) = LinSpaceIndex(i+1);
            end
        end
    end
else
    BCs.LB.InsertIndex = [];
    BCs.UB.InsertIndex = [];
end


% concatenation format GA
LB = [BCs.LB.Nply;   BCs.LB.Thetas;   BCs.LB.BalancedLoc;   BCs.LB.TenPercentLoc;  BCs.LB.MidPlane;   BCs.LB.InsertIndex];     % design variable lower bounds
UB = [BCs.UB.Nply;   BCs.UB.Thetas;   BCs.UB.BalancedLoc;   BCs.UB.TenPercentLoc;  BCs.UB.MidPlane;   BCs.UB.InsertIndex];     % design variable upper bounds



% if Damtol, limit the first ply angle to +- 45 (and for the last ply as well py for generic laminates) - coded using 1 and 2
if Constraints.Vector(3)
    LB(sum(NStruct.NpatchVar)+1) = 1; % first theta (-45)
    UB(sum(NStruct.NpatchVar)+1) = 2; % first theta (+45)
    
    BCs.LB.Thetas(1) = 1;
    BCs.UB.Thetas(1) = 2;
    
    if strcmp(LamType,'Generic')
        LB( sum(NStruct.NpatchVar) + NStructMin.NthetaVar) =  1; % last theta
        UB( sum(NStruct.NpatchVar) + NStructMin.NthetaVar) =  2; % last theta
        
        BCs.LB.Thetas(NStructMin.NthetaVar) = 1;
        BCs.UB.Thetas(NStructMin.NthetaVar) = 2;
    end
end

end