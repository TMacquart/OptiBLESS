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
%      Computes the number of design variables for each genotype Field
%
% NStruct = Attribute_NDvs(Nplies,Constraints,LamType)
%
% NStruct.NthetaVar;        % number of theta's  design variables
% NStruct.NbalVar;          % number of balanced design variables (will always be equal to NthetaVar, but kept for sake of clarity)
% NStruct.N10percentVar;    % number of 10% rule design variables
% NStruct.NMidPlane;        % number of MidPlane design variables
% NStruct.NInsertVar;         % number of Drop Location design variables
% =====                                                              ======



function NStruct = Attribute_NDvs(Nplies,Constraints,LamType)

MaxNplies = max(Nplies(:));


%% Generic Laminates
if strcmp(LamType,'Generic')
    NbalVar   = 0;  % no balanced design variables
    NMidPlane = 0;  % no mid plane design variables
    
    if Constraints.Vector(4)==0
        % no 10% rule, no 10% design variables
        N10percentVar = 0;
        NthetaVar     = MaxNplies;
    else
        % 10% rule, about 40% (or more) of design variable are allocated to the 10% rule
        N10percentVar = (round( ceil(0.4*MaxNplies) /4) *4); % is a multiple of 4
        NthetaVar     = MaxNplies-N10percentVar;
    end
    
    NtotalPly = NthetaVar + N10percentVar; % used for internal check ( NtotalPly must be = to MaxNplies)
end



%% Symmetric Laminates
if strcmp(LamType,'Sym')
    NbalVar   = 0; % no balanced design variables
    
    if Constraints.Vector(4)==0
        % no 10% rule, no 10% design variables
        N10percentVar = 0;
        if rem(MaxNplies,2) == 0
            % even number of plies
            NthetaVar = MaxNplies/2;
            NMidPlane = 0;
        else
            % odd number of plies
            NthetaVar = floor(MaxNplies/2);
            NMidPlane = 1;                   % add a midplane ply
        end
        
    else
        % 10% rule, about 40% (or more) of design variable are allocated to the 10% rule
        N10percentVar = (round( ceil(0.4*MaxNplies/2) /4) *4); % is a multiple of 4
        if rem( (MaxNplies-N10percentVar*2) ,2) == 0
            % even number of plies left
            NthetaVar = (MaxNplies-N10percentVar*2)/2;
            NMidPlane = 0;
        else
            % odd number of plies left
            NthetaVar = floor((MaxNplies-N10percentVar*2)/2);
            NMidPlane = 1;
        end
        
    end
    
    NtotalPly = 2*(NthetaVar + NbalVar + N10percentVar) + NMidPlane;   % used for internal check ( NtotalPly must be = to MaxNplies)
end



%% Balanced Laminates
if strcmp(LamType,'Balanced')
    
    if Constraints.Vector(4)==0
        % no 10% rule, no 10% design variables
        N10percentVar = 0;
        if rem(MaxNplies,2) == 0
            % even number of plies
            NthetaVar = MaxNplies/2;
            NbalVar   = MaxNplies/2;
            NMidPlane = 0;
        else
            % odd number of plies
            NthetaVar = floor(MaxNplies/2);
            NbalVar   = floor(MaxNplies/2);
            NMidPlane = 1;
        end
        
    else
        % 10% rule, about 40% (or more) of design variable are allocated to the 10% rule
        N10percentVar = (round( ceil(0.4*MaxNplies) /4) *4); % is a multiple of 4
        if rem(MaxNplies-N10percentVar,2) == 0
            % even number of plies left
            NthetaVar = (MaxNplies-N10percentVar)/2;
            NbalVar   = (MaxNplies-N10percentVar)/2;
            NMidPlane = 0;
        else
            % odd number of plies left
            NthetaVar = floor((MaxNplies-N10percentVar)/2);
            NbalVar   = floor((MaxNplies-N10percentVar)/2);
            NMidPlane = 1;
        end
        
    end
    
    NtotalPly = NthetaVar + NbalVar + N10percentVar + NMidPlane;
end



%% Symmetric and Balanced Laminate
if strcmp(LamType,'Balanced_Sym')
    
    if Constraints.Vector(4)==0
        % no 10% rule
        N10percentVar = 0;
        if rem(MaxNplies,4) == 0
            % is a multiple of 4
            NthetaVar = MaxNplies/4;
            NbalVar   = MaxNplies/4;
            NMidPlane = 0;
        else
            NthetaVar = floor(MaxNplies/4);
            NbalVar   = floor(MaxNplies/4);
            NMidPlane = MaxNplies-2*(NthetaVar+NbalVar);
        end
        
    else
        
        N10percentVar = (round( ceil(0.4*MaxNplies) /4) *4) /2;
        if rem(N10percentVar,4)~=0
            N10percentVar = N10percentVar + rem(N10percentVar,4);
        end
        
        if rem(MaxNplies - N10percentVar*2,4) == 0
            % left number of plies is a multiple of 4
            NthetaVar = (MaxNplies - N10percentVar*2)/4;
            NbalVar   = (MaxNplies - N10percentVar*2)/4;
            NMidPlane = 0;
        else
            NthetaVar = floor((MaxNplies - N10percentVar*2)/4);
            NbalVar   = floor((MaxNplies - N10percentVar*2)/4);
            NMidPlane = MaxNplies -2*(NthetaVar+NbalVar+N10percentVar);
            
            if NMidPlane>=2
                N10percentVar = N10percentVar + floor(NMidPlane/2);
                NMidPlane = MaxNplies -2*(NthetaVar+NbalVar+N10percentVar);
            end
        end
    end
    
    NtotalPly = 2*(NthetaVar + NbalVar + N10percentVar)+ NMidPlane;
end



%% Internal Checks

if NthetaVar<=0
    error('The Laminate does not contain enough plies (N10percentVar<=0). Deactivate the constraint or increase the number of plies.')
end
if NtotalPly~=MaxNplies
    keyboard % should never happen (used for check)
end
if NMidPlane<0
    keyboard  % should never happen (used for check)
end




%% Compute the number of drop ply design variables: NInsertVar
% note that NInsertVar does not have to be the exact number, it can be
% greater than the actual number of drops

NInsert = MaxNplies - min(Nplies(:)); % total insert in the Stacking sequence table

if NInsert~=0
    if strcmp(LamType,'Generic'),
        NInsertVar = NInsert;
    end
    
    
    if strcmp(LamType,'Sym') || strcmp(LamType,'Balanced')
        if rem(NInsert,2) == 0
            NInsertVar = NInsert/2;
        else
            NInsertVar = floor(NInsert/2) + 1;
        end
    end
    
    if strcmp(LamType,'Balanced_Sym')
        if rem(NInsert,4) == 0
            NInsertVar = NInsert/4;
        else
            NInsertVar = floor(NInsert/4) + 1;
        end
    end
    
    if NMidPlane~=0
        NInsertVar = NInsertVar + 1;  % account for possible mid-plane drops
    end
    
    if Constraints.Vector(2)
        NInsertVar = NInsertVar + round(N10percentVar/2);
    end
    
    NInsertVar = round(1.25*NInsertVar);
else
    % no drops
    NInsertVar = 0 ;
end


NStruct.MinNplies       = min(Nplies(:));
NStruct.MaxNplies       = MaxNplies;        % Maximum number of plies
NStruct.NthetaVar       = NthetaVar;        % number of theta's  design variables
NStruct.NbalVar         = NbalVar;          % number of balanced design variables
NStruct.N10percentVar   = N10percentVar;    % number of 10% rule design variables
NStruct.NMidPlane       = NMidPlane;        % number of MidPlane design variables
NStruct.NInsertVar      = NInsertVar;         % number of Drop Location design variables

