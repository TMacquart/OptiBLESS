% =====                                                              ==== 
%           Check the feasibility of a blended stacking sequence
%            
% Non-feasible individual are attributed a penalised fitness values
% GuideAngles must be input in degrees!
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



function [FEASIBLE] = Check_Feasibility(ConstraintVector,GuideAngles,ShuffleLoc,DropsIndexes,NGuidePlies,NDropPlies,LamType,NContiguity)

FEASIBLE = true;

if (length(DropsIndexes) ~= length(unique(DropsIndexes))),
    FEASIBLE = false; % penalise GA individual with non-unique ply drops
    return
end 

if length(DropsIndexes) ~= sum(NDropPlies)
    FEASIBLE = false; 
    return
end


if max(DropsIndexes)>length(GuideAngles) 
    FEASIBLE = false; 
    return
end 

if ConstraintVector(1)
    if abs(GuideAngles(1)) ~= 45
        FEASIBLE = false;
        return
    end
    if strcmp(LamType,'Generic') && abs(GuideAngles(end)) ~= 45
        FEASIBLE = false;
        return
    end
end


if ConstraintVector(3) || ConstraintVector(4) % Disorientation and Contiguity

    for iDrop= 0:length(DropsIndexes)
        
        if iDrop == 0
            [FiberAngles] = Convert_dvAngles2FiberAngles(GuideAngles,[],ShuffleLoc,LamType);
        else
            DropsLoc = unique(DropsIndexes(1:iDrop));
            DropsLoc(DropsLoc>NGuidePlies) = [];
            FiberAngles = Convert_dvAngles2FiberAngles(GuideAngles,DropsLoc,ShuffleLoc,LamType);
        end
        
        DetlaAngle = zeros(numel(FiberAngles)-1,1);
        Contiguity = 0;
        for iply = 1:numel(FiberAngles)-1
            if FiberAngles(iply)>=-45 && FiberAngles(iply)<=45
                DetlaAngle(iply) = abs(FiberAngles(iply)-FiberAngles(iply+1));
            elseif FiberAngles(iply)>45
                if FiberAngles(iply+1)>-45
                    DetlaAngle(iply) = abs(FiberAngles(iply)-FiberAngles(iply+1));
                else
                    DetlaAngle(iply) = abs(FiberAngles(iply)-(180+FiberAngles(iply+1)));
                end
            elseif FiberAngles(iply)<-45
                if FiberAngles(iply+1)<45
                    DetlaAngle(iply) = abs(FiberAngles(iply)-FiberAngles(iply+1));
                else
                    DetlaAngle(iply) = abs(FiberAngles(iply)-(-180+FiberAngles(iply+1)));
                end
            end
            
            
            if ConstraintVector(3) && DetlaAngle(iply)>45 % disorientation
                FEASIBLE = false;
                return
            end
            
            if ConstraintVector(4) % contiguity
                if DetlaAngle(iply)==0
                    Contiguity = Contiguity + 1;
                else
                    Contiguity = 0;
                end
                if Contiguity >= NContiguity
                    FEASIBLE = false;
                    return
                end
            end
        end

    end
end


if (ConstraintVector(7) || ConstraintVector(6)) && ~isempty(DropsIndexes)    % check covering and internal continuity constraints
   
    DropsLoc = sort(DropsIndexes);
    
    if ConstraintVector(7) && ~isempty(find(DropsLoc == 1,1)), % covering check first ply is not removed
        keyboard % this constraints should be removed, it has been replaced by a direct constraints on LB , UB
        FEASIBLE = false;  
        return
    end 
    
    if ConstraintVector(7) && strcmp(LamType,'Generic') && ~isempty(find(DropsLoc == length(GuideAngles),1)), % covering check first ply is not removed
        keyboard % this constraints should be removed, it has been replaced by a direct constraints on LB , UB
        FEASIBLE = false;
        return
    end
    
    if ConstraintVector(6)
        for i = 1 : length(DropsLoc)-3
            deltaLoc = diff(DropsLoc(i:i+3));
            if sum(abs(deltaLoc))==3
                FEASIBLE = false;
                return
            end
        end
    end
end



if ConstraintVector(5) % check if balanced for indirect constraint handling
    % guide
    [FiberAngles] = Convert_dvAngles2FiberAngles(GuideAngles,[],LamType);
    UniqueAngle = unique(FiberAngles);
    UniqueAngle(abs(UniqueAngle)==90) = [];
    UniqueAngle(UniqueAngle==0)       = [];
    for j=1:length(UniqueAngle)
        if length(find(UniqueAngle(j)==FiberAngles)) ~= length(find(-UniqueAngle(j)==FiberAngles))
            FEASIBLE = false;
            return
        end
    end
    
    % guideDroped laminates
    for iDrop=1:length(DropsIndexes)
        ply_angles = GuideAngles;
        DropsLoc = unique(DropsIndexes(1:iDrop));
        if max(DropsLoc)>NGuidePlies
            DropsLoc(DropsLoc>NGuidePlies) = [];
        end
        
        ply_angles(DropsLoc(DropsLoc<=length(ply_angles))) = [];  % drop plies
        FiberAngles      = Convert_dvAngles2FiberAngles(ply_angles,[],LamType);
        UniqueAngle = unique(FiberAngles);
        UniqueAngle(abs(UniqueAngle)==90) = [];
        UniqueAngle(UniqueAngle==0)       = [];
        
        for j=1:length(UniqueAngle)
            if length(find(UniqueAngle(j)==FiberAngles)) ~= length(find(-UniqueAngle(j)==FiberAngles))
                FEASIBLE = false;
                return
            end
        end
    end
end


end