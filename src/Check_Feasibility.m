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
% =====                                                              ==== 
%           Check the feasibility of a blended stacking sequence
%            
% Non-feasible individual are attributed a penalised fitness values
% =====                                                              ==== 




function [FEASIBLE,output] = Check_Feasibility(Constraints,SSTable,NplySS,SSTableSymbolic,varargin)

% Only implicit constraints are check in this functions, 
% constraints not checked are explicitly coded and are therefore ensured

output.ConstViolated = '';
FEASIBLE = true;


% if NplySS not provided, compute it
if nargin==2 || isempty(NplySS)
    NplySS = nan*ones(size(SSTable,1),1);
    for j = 1:size(SSTable,1)
        NplySS(j) = length(find(~cellfun(@isempty,SSTable(j,:))));
    end
end

% Internal Continuity
if Constraints.Vector(7)
    SSDrops  = sort(find(cellfun(@isempty,SSTable(end,:)))); % check last line of SS Table
    if length(SSDrops)>Constraints.NInternalCont
        for i = 1 : length(SSDrops)-Constraints.NInternalCont
            deltaLoc = diff(SSDrops(i:i+Constraints.NInternalCont));
            if sum(abs(deltaLoc))==Constraints.NInternalCont
                FEASIBLE = false;
                output.ConstViolated = 'InternalContinuity';
                return
            end
        end
    end
end


% These Constraints are checked for all patches
% Rule10percent,  Disorientation and Contiguity
if Constraints.Vector(4) || Constraints.Vector(5) || Constraints.Vector(6)
    
    for i = 1:length(NplySS) 
        FiberAngles = cell2mat(SSTable(i,:));
        
        % --- Rule10percent 
            if 1 %Constraints.Vector(2) 
                
                TenpercentDV = round(length(FiberAngles)*0.1);
                N0Plies      = length(find(FiberAngles==0));
                N90Plies     = length(find(abs(FiberAngles)==90));
                N45Plies     = length(find(FiberAngles==45));
                NM45Plies    = length(find(FiberAngles==-45));
                
                if N0Plies<TenpercentDV || N90Plies<TenpercentDV || N45Plies<TenpercentDV || NM45Plies<TenpercentDV
                    FEASIBLE = false;
                    output.ConstViolated = 'TenPercentRule';
                    output.TenPercentVec = [N0Plies<TenpercentDV  N90Plies<TenpercentDV  N45Plies<TenpercentDV  NM45Plies<TenpercentDV];
                    return
                end

            end
        % ---
        
        
        DetlaAngle  = ComputeDeltaAngle(FiberAngles); % delta fibre angles between plies
        
        % --- Disorientation
            if Constraints.Vector(5) && ~isempty( find(DetlaAngle>45,1) ) 
                FEASIBLE = false;
                output.ConstViolated = 'Disorientation';
                
                if nargin == 4 % for inipop only
                    % symbolic for Inipop Gen. only, find where
                    % disorientation occured and return to repair
                        
                    DvIndex    = find(DetlaAngle>45,1);
                    DvIndex_SS = find(~cellfun(@isempty,SSTable(i,:)));
                    Dvs        = cell2mat(SSTableSymbolic(i,DvIndex_SS(DvIndex+[0 1])));

                    if find(Dvs>0,1)
                        Index = find(Dvs>0);
                        if length(Index)==1
                            output.FailedDv = Dvs(find(Dvs>0,1));           % failed due to a theta
                        else
                            output.FailedDv = Dvs(2);      
                        end
                    elseif find(Dvs<0,1)
                        Index = find(Dvs<0);
                        if length(Index)==1
                            output.FailedDv = abs(Dvs(find(Dvs<0,1)));           % failed due to a theta
                        else
                            output.FailedDv = abs(Dvs(2));      
                        end
                    end
                end
                return
            end
        % ---
        
        
        % --- Contiguity
            if Constraints.Vector(6) 
                Contiguity = 0;
                for iply = 1:numel(FiberAngles)-1

                    if DetlaAngle(iply)==0
                        Contiguity = Contiguity + 1;
                    else
                        Contiguity = 0;
                    end
                    if Contiguity > (Constraints.NContiguity-1)
                        output.ConstViolated = 'Contiguity';
                        FEASIBLE = false;
                        return
                    end
                end
            end
        % ---
        
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