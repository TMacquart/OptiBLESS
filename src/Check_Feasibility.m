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



function [FEASIBLE,ConstViolated] = Check_Feasibility(Constraints,NpliesPerLam,SSTable)

ConstViolated = '';
FEASIBLE = true;

NpliesPerLam = sort(unique(NpliesPerLam),'descend');
for j = 1:size(SSTable,1)
    NplySS(j) = length(find(~cellfun(@isempty,SSTable(j,:))));
end
if find(ismember(NpliesPerLam,NplySS)==0,1)
    FEASIBLE = false;
    ConstViolated = 'Drop Indexes';
    return
end


% if Constraints.Vector(1) % damtol
%     if abs(SSTable{1,1}) ~= 45
%         keyboard % should never happem and can be removed
%         FEASIBLE = false;
%         ConstViolated = 'Damtol';
%         return
%     end
%     if strcmp(LamType,'Generic') && abs(GuideAngles(end)) ~= 45
%         FEASIBLE = false;
%         ConstViolated = 'Damtol';
%         return
%     end
% end


if Constraints.Vector(2) || Constraints.Vector(3) || Constraints.Vector(4) % 10%rule, Disorientation and Contiguity
    

    for i = 1:length(NplySS) 
        FiberAngles = cell2mat(SSTable(j,:));
        DetlaAngle = ComputeDeltaAngle(FiberAngles);
        Contiguity = 0;
        for iply = 1:numel(FiberAngles)-1

            if Constraints.Vector(2)
                FEASIBLE = Check_10PercentRule(FiberAngles);
                if ~FEASIBLE
                    ConstViolated = 'TenPercentRule';
                    return
                end
            end
                
            if Constraints.Vector(3) && DetlaAngle(iply)>45 % disorientation
                FEASIBLE = false;
                ConstViolated = 'Disorientation';
                return
            end
            
            if Constraints.Vector(4) % contiguity
                if DetlaAngle(iply)==0
                    Contiguity = Contiguity + 1;
                else
                    Contiguity = 0;
                end
                if Contiguity >= Constraints.Contiguity
                    ConstViolated = 'Contiguity';
                    FEASIBLE = false;
                    return
                end
            end
        end

    end
end

if Constraints.Vector(6)
    FEASIBLE = Check_InternalContinuity(SSTable);
    if ~FEASIBLE
          ConstViolated = 'InternalContinuity';
        return
    end
end

if Constraints.Vector(7)  % check covering 
   if ~isempty(find(isempty(SSTable(:,1)),1))
        ConstViolated = 'Covering';
       FEASIBLE = false;
       return
   end
end



% if ConstraintVector(5) % check if balanced for indirect constraint handling
%     % guide
%     [FiberAngles] = Convert_dvAngles2FiberAngles(GuideAngles,[],LamType);
%     UniqueAngle = unique(FiberAngles);
%     UniqueAngle(abs(UniqueAngle)==90) = [];
%     UniqueAngle(UniqueAngle==0)       = [];
%     for j=1:length(UniqueAngle)
%         if length(find(UniqueAngle(j)==FiberAngles)) ~= length(find(-UniqueAngle(j)==FiberAngles))
%             FEASIBLE = false;
%             return
%         end
%     end
%     
%     % guideDroped laminates
%     for iDrop=1:length(DropsIndexes)
%         ply_angles = GuideAngles;
%         DropsLoc = unique(DropsIndexes(1:iDrop));
%         if max(DropsLoc)>NGuidePlies
%             DropsLoc(DropsLoc>NGuidePlies) = [];
%         end
%         
%         ply_angles(DropsLoc(DropsLoc<=length(ply_angles))) = [];  % drop plies
%         FiberAngles      = Convert_dvAngles2FiberAngles(ply_angles,[],LamType);
%         UniqueAngle = unique(FiberAngles);
%         UniqueAngle(abs(UniqueAngle)==90) = [];
%         UniqueAngle(UniqueAngle==0)       = [];
%         
%         for j=1:length(UniqueAngle)
%             if length(find(UniqueAngle(j)==FiberAngles)) ~= length(find(-UniqueAngle(j)==FiberAngles))
%                 FEASIBLE = false;
%                 return
%             end
%         end
%     end
% end


end