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
%            Creates Initial Population of Ply angles for SST             %
%                                                                         %
% Recommended for heavily constrained problems                            %
%
% Requires 4 inputs
% nvar       : # a design variables = half (the # of plies + # of drops)
% Npop       : Size of the initial population
% NDropPlies : Vector containing the # of plies to drop per section
% EDoutput   : true or false, enable/disable EuclideanDist calculation
%
% Returns up to 3 outputs 
% IniPop        : Initial population matrix ( Stacking Sequence )
% IniPopLP      : Initial population LP     
% EuclideanDist : Euclidean Distance between the LP of each Section
% =====                                                              ==== %

function [IniPop,IniPopLP,EuclideanDist] = SST_CreateBlendedIniPop (nvar,Npop,Nmax,Nmin,Constraints,AllowedNplies,LamType)


% NDropPlies is the number of plies to drop in the symmetric part only

if sum(arrayfun(@(i) length(AllowedNplies{i}),1:length(AllowedNplies) )) == length(AllowedNplies)
    ConstantThickness = true;
    IniPop = zeros(Npop,nvar + length(AllowedNplies));
else
    ConstantThickness = false;
    IniPop = zeros(Npop,nvar);
end

fprintf(' Creating IniPop ... ' )
ConstraintVector = Constraints.Vector;
DeltaAngle       = Constraints.DeltaAngle;

EuclideanDist = cell(Npop,1);
IniPopLP      = cell(Npop,1);
ipop          = 1 ;
NTried        = 0;

while ipop < Npop + 1

    NpliesPerLam = zeros(1,length(AllowedNplies));
    for iply = 1:length(AllowedNplies)
        NpliesPerLam(iply) = AllowedNplies{iply}(randi([1 length(AllowedNplies{iply})],1,1));
    end
    NpliesPerLam  = sort(NpliesPerLam,'descend'); 
    NPliesGuide   = max(NpliesPerLam);
    
    FEASIBLE    = true;
    GuideAngles = zeros(1,NPliesGuide);
%     GuideAngles = nan*ones(1,Nmax);
    
    % ---
    if 1 % enforce constraints on IniPop - Build the angles step by step
        

        if ConstraintVector(1)                                              % Damtol
            A = [-1 1];
            GuideAngles(1) = 45*A(ceil(2*rand)); % 1st ply is +- 45
        else
            GuideAngles(1) = -90 + 180*rand;
        end
        
        
        for iAngle = 2 : NPliesGuide
            if ConstraintVector(3)
                if ConstraintVector(4)
                    A = [(-45 + 40*rand) (5 + 40*rand)];
                    AddedAngle = GuideAngles(iAngle-1) + A(ceil(2*rand));
                    while abs(AddedAngle)>90
                        AddedAngle = GuideAngles(iAngle-1) + A(ceil(2*rand)); % not too close (min +5) not too far (max 45)
                    end
                else
                    AddedAngle = GuideAngles(iAngle-1) + (-45 + 90*rand);
                    while abs(AddedAngle)>90
                        AddedAngle = GuideAngles(iAngle-1) + (-45 + 90*rand); % not too far (max 45)
                    end
                end
            else
                if ConstraintVector(4)
                    A = [(-90 + 85*rand) (5 + 85*rand)];
                    AddedAngle = GuideAngles(iAngle-1) + A(ceil(2*rand));
                    while abs(AddedAngle)>90
                        AddedAngle = GuideAngles(iAngle-1) + A(ceil(2*rand)); % not too close (min +5)
                    end
                else
                    AddedAngle = GuideAngles(iAngle-1) + (-90 + 180*rand);    % no constraint - full range
                    while abs(AddedAngle)>90
                        AddedAngle = GuideAngles(iAngle-1) + (-90 + 180*rand);
                    end
                end
            end
            GuideAngles(iAngle) = AddedAngle;
        end                                     % Stack plies according to constraints activated
        
        if ConstraintVector(5)
            GuideAngles = round((GuideAngles)/DeltaAngle)*DeltaAngle; % made multiples of DeltaAngle
        end
        
        if ConstraintVector(2) % 10% rule 
             [GuideAngles] = Enforcing_10PercentRule(GuideAngles);
        end
        
    end
    % ---
    
    if FEASIBLE

        NdropPlies = (NPliesGuide)-min(NpliesPerLam);
        
        if NdropPlies>0
            DropsIndexes = Generating_DropIndexes (GuideAngles,NdropPlies,ConstraintVector);
            if isempty(DropsIndexes)
                FEASIBLE = false;
            end
        else
            DropsIndexes = [];
        end
        
        if FEASIBLE
            [FEASIBLE] = CheckFeasibility(ConstraintVector,GuideAngles,DropsIndexes,LamType);
        end

        DropsIndexes = [DropsIndexes nan*ones(1,(Nmax-Nmin)-length(DropsIndexes))];

        if ConstraintVector(5)
            GuideAngles = round((90+GuideAngles)/DeltaAngle);
        end
        
        GuideAngles = [GuideAngles nan*ones(1,Nmax-length(GuideAngles)) ]; %#ok<AGROW>
        
        if FEASIBLE
            IniPop(ipop,:)    = [NpliesPerLam GuideAngles DropsIndexes];
            ipop = ipop + 1;
        end
    end
    
    NTried = NTried + 1;
    if  NTried == 5000 && ipop<2
        error('Hard Constrained Problem. Difficulties Generating IniPop')
    end
    
end

if ConstantThickness
    IniPop(:,1:length(AllowedNplies)) = []; % remove thickness variables
end
fprintf(' IniPop Created ... ' )
end