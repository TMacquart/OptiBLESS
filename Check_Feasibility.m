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



function [FEASIBLE] = CheckFeasibility(ConstraintVector,GuideAngles,DropsIndexes,NGuidePlies,NDropPlies,LamType)

FEASIBLE = true;

if length(GuideAngles) ~= NGuidePlies
    FEASIBLE = false; 
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


if (length(DropsIndexes) ~= length(unique(DropsIndexes))),
    FEASIBLE = false; % penalise GA individual with non-unique ply drops
    return
end 


if ConstraintVector(1)
    if abs(GuideAngles(1)) ~= 45
        FEASIBLE = false;
        return
    end
end


if ConstraintVector(3) % Disorientation
    [FiberAngles] = Convert_dvAngles2FiberAngles(GuideAngles,LamType);
    for iply = 1:numel(FiberAngles)-1
        if FiberAngles(iply)>=-45 && FiberAngles(iply)<=45
            DetlaAngle = abs(FiberAngles(iply)-FiberAngles(iply+1));
        elseif FiberAngles(iply)>45
            if FiberAngles(iply+1)>=-45
                DetlaAngle = abs(FiberAngles(iply)-FiberAngles(iply+1));
            else
                DetlaAngle = abs(FiberAngles(iply)-(180+FiberAngles(iply+1)));
            end
        elseif FiberAngles(iply)<-45
            if FiberAngles(iply+1)<=45
                DetlaAngle = abs(FiberAngles(iply)-FiberAngles(iply+1));
            else
                DetlaAngle = abs(FiberAngles(iply)-(-180+FiberAngles(iply+1)));
            end
        end
        if DetlaAngle>45
            FEASIBLE = false;
            return
        end
    end
end


if ConstraintVector(4) % Contiguity
    DeltaOrientation = abs(diff(GuideAngles));
    if ~isempty(find(DeltaOrientation<5,1)),
        FEASIBLE = false;
        return
    end
end


if (ConstraintVector(7) || ConstraintVector(6)) && ~isempty(DropsIndexes)    % check covering and internal continuity constraints
   
    DropsLoc = sort(DropsIndexes);
    
    if ConstraintVector(7) && ~isempty(find(DropsLoc == 1,1)), % covering check first ply is not removed
        FEASIBLE = false;  
        return
    end 
    
    if ConstraintVector(7) && strcmp(LamType,'Generic') && ~isempty(find(DropsLoc == length(GuideAngles),1)), % covering check first ply is not removed
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

end