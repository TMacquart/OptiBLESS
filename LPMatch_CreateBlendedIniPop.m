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
%            Creates Initial Population of Ply angles for GA              %
%                                                                         %
% Recommended for heavily constrained problems                            %
%
% Requires 4 inputs
% Nvar       : # a design variables = half (the # of plies + # of drops)
% Npop       : Size of the initial population
% NDropPlies : Vector containing the # of plies to drop per section
% EDoutput   : true or false, enable/disable EuclideanDist calculation
%
% Returns up to 3 outputs 
% IniPop        : Initial population matrix ( Stacking Sequence )
% IniPopLP      : Initial population LP     
% EuclideanDist : Euclidean Distance between the LP of each Section
% =====                                                              ==== %

function [IniPop,IniPopLP,EuclideanDist] = LPMatch_CreateBlendedIniPop (Npop,NguidePlies,NguideDrops,Constraints,LamType)

Nvar = NguidePlies + sum(NguideDrops);
ConstraintVector = Constraints.Vector;
DeltaAngle       = Constraints.DeltaAngle;

fprintf(' Creating IniPop ... ' )

EuclideanDist = cell(Npop,1);
IniPopLP      = cell(Npop,1);
IniPop        = zeros(Npop,Nvar);

ipop          = 1 ;
NTried        = 0;


while ipop < Npop + 1
    
    GuideAngles = zeros(1,NguidePlies);
    
    % ---
    if 1 % enforce constraints on IniPop - Build the angles step by step
        
        if ConstraintVector(1)                                              % Damtol
            A = [-1 1];
            GuideAngles(1) = 45*A(ceil(2*rand)); % 1st ply is +- 45
        else
            GuideAngles(1) = -90 + 180*rand;
        end
        
        
        for iAngle = 2 : NguidePlies
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
            else % no disorientation
                if ConstraintVector(4)
                    A = [(-90 + 85*rand) (5 + 85*rand)];
                    AddedAngle = GuideAngles(iAngle-1) + A(ceil(2*rand));
                    while abs(AddedAngle)>90
                        AddedAngle = GuideAngles(iAngle-1) + A(ceil(2*rand)); % not too close (min +5)
                    end
                else
                     AddedAngle = -90 + 180*rand;    % no constraint - full range
                end
            end
            GuideAngles(iAngle) = AddedAngle;
        end                                     % Stack plies according to constraints activated
        
        if ConstraintVector(2) % 10% rule                                       
            [GuideAngles] = Enforcing_10PercentRule(GuideAngles);
        end
        
        if ConstraintVector(5)
            GuideAngles = round((GuideAngles)/DeltaAngle)*DeltaAngle; % made multiples of DeltaAngle
        end
        
    end
    % ---
    
    FEASIBLE =true;
    if sum(NguideDrops)>0
        DropsIndexes = Generating_DropIndexes (GuideAngles,sum(NguideDrops),ConstraintVector);
        if isempty(DropsIndexes)
            FEASIBLE = false;
        end
    else
        DropsIndexes = [];
    end
    
    if FEASIBLE
        [FEASIBLE] = CheckFeasibility(ConstraintVector,GuideAngles,DropsIndexes,LamType);
    end
    
%     if ipop == 1
%         GuideAngles  = [+45 -45 90 0 45 90 0]
%         DropsIndexes = [2 4]
%     end
    
    % formating angles in design variables
    if ConstraintVector(5) % if DiscreteAngles = use integer indexes
        if ~ConstraintVector(1)
            GuideAngles = round((90+GuideAngles)/DeltaAngle);
        else
            GuideAngles(2:end) = round((90+GuideAngles(2:end))/DeltaAngle);
        end
    end
    if ConstraintVector(1) % make the first ply a discrete value with only 2 possible state
        if (GuideAngles(1) == -45), GuideAngles(1) = 1; end
        if (GuideAngles(1) == 45),  GuideAngles(1) = 2; end
    end
    
    if FEASIBLE
        IniPop(ipop,:)    = [GuideAngles DropsIndexes]; % here GuideAngles is not necessarily in degrees
        ipop = ipop + 1;
    end
    
    
    NTried = NTried + 1;
    if  NTried == 5000 && ipop<2
        error('Hard Constrained Problem. Difficulties Generating IniPop')
    end
    
end



fprintf(' IniPop Created ... ' )
end