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
%                            LP Evaluation function                       %
%                          Used as fitness fct by the GA                  %
%
%  The individual is composed of 3 parts:
%  --------------------------------------------
%  [ [Nply(1) ... Nply(Nsec)]                   -- the Number of plies
%  [ Theta(1) ... Theta(N)   nan ... nan]       -- N is the guide Nply
%  [ Drop1    ... DropM]                      
% =====                                                              ==== %

function [fitness,LP,output] = Eval_Fitness (Individual,Constraints,IndexLp,SortedObj,Nsec,AllowedNplies,LamType)

FEASIBLE  = true;  
LP_obj    = SortedObj.LP_obj;
ImportanceFactor = SortedObj.ImportanceFactor;

if Constraints.NRange == 1
    NpliesperLam = cell2mat(AllowedNplies);
else
    NpliesperLam = Individual(1:Nsec);
end

if Constraints.ORDERED
    NpliesperLam = sort(NpliesperLam,'descend');                            % repair to ensure thickness is ordered
    SortIndex    = 1:length(NpliesperLam);
    if find(diff(NpliesperLam)>0,1)
        FEASIBLE = false;
    end
else
    [NpliesperLam,SortIndex] = sort(NpliesperLam,'descend');
end


NGuidePlies  = max(NpliesperLam);                                           % number of plies in the guide laminate (take half for Sym.)
NDropPlies   = abs(diff(NpliesperLam));                                     % number of ply drops
GuideAngles  = Individual(Nsec + [1:NGuidePlies]);
GuideAngles  = GuideAngles(~isnan(GuideAngles));



% --- organise ply drop sequences
Ndrop = length(NDropPlies);
StartIndex   = Nsec+max(cell2mat(AllowedNplies));
DropIndexes = cell(1,Ndrop);

for iDrop = 1 : Ndrop
    DropIndexesTEMP     = Individual(StartIndex + [1:NDropPlies(iDrop)]);
    DropIndexes{iDrop}  = DropIndexesTEMP(~isnan(DropIndexesTEMP));
    StartIndex          = StartIndex + NDropPlies(iDrop);
end


ConstraintVector = Constraints.Vector;
if  ConstraintVector(5)     % if DiscreteAngles
    if ~ConstraintVector(1)     % if not Damtol
        GuideAngles = -90 + GuideAngles*Constraints.DeltaAngle;
    else
        GuideAngles(2:end) = -90 + GuideAngles(2:end)*Constraints.DeltaAngle;
    end
    if ~isempty(find(abs(GuideAngles)>90,1)), error('Angle greater than 90 detected');   end
end


% --- check / enforce constraints for ply Guide 
if ConstraintVector(1) % if Damtol, first angle is +- 45
    A = [-1 1];
    GuideAngles(1) = 45*A(GuideAngles(1));
    if strcmp(LamType,'Generic') % convert to closest
        if abs(GuideAngles(end)-45) < abs(GuideAngles(end) +45) % closer to +45
            GuideAngles(end) = 45;
        else
            GuideAngles(end) = -45;
        end
    end
end

if ConstraintVector(2)                                                  % 10% rule (0,90,45)
    [GuideAngles] = Enforce_10PercentRule(GuideAngles);
end

if FEASIBLE
    [FEASIBLE] = Check_Feasibility(ConstraintVector,GuideAngles,cell2mat(DropIndexes),NGuidePlies,NDropPlies,LamType);
end




% ---
if 1    % fitness calculation

    FiberAngles  = Convert_dvAngles2FiberAngles(GuideAngles,LamType);

    LP(:,SortIndex(1))  = Convert_SS2LP(FiberAngles);            % evaluate lamination parameters for the guide
    SS{SortIndex(1)}    = FiberAngles;
    
    for iDrop = 1 : Ndrop
        index = SortIndex(iDrop+1);
        
        ply_angles = GuideAngles;
        
        DropsLoc = unique(cell2mat(DropIndexes(1:iDrop)));
        if max(DropsLoc)>NGuidePlies
            DropsLoc(DropsLoc>NGuidePlies) = [];
        end
   
        ply_angles(DropsLoc(DropsLoc<length(ply_angles))) = [];  % drop plies

        FiberAngles      = Convert_dvAngles2FiberAngles(ply_angles,LamType);
        SS{index}        = FiberAngles;
        
        LP(:,index) = Convert_SS2LP(FiberAngles);          % evaluate lamination parameters for the droped laminates
    end
end



% --- select only the LP corresponding to the NpliesperLam (doubled if any)

fitness1 = Constraints.FitnessFct(LP(:,SortedObj.sortIndex(SortIndex)),NpliesperLam(SortedObj.sortIndex(SortIndex)));


localFit = zeros(length(NpliesperLam),1);
for ilam = 1 : length(NpliesperLam)
    localFit(ilam) = norm(LP_obj(IndexLp,ilam) - LP(IndexLp,ilam));
end
fitness = sum(localFit.*ImportanceFactor);

% --- to remove
%         if 1
%             localFit = zeros(length(NpliesperLam),1);
%             for ilam = 1 : length(NpliesperLam)
%                 localFit(ilam) = rms(LP_obj(IndexLp,ilam) - LP(IndexLp,ilam));
%             end
%             fitRMS = sum(localFit.*ImportanceFactor);
%         end
% % ---


fitness = fitness ...
        + (Constraints.alpha)*sum(NpliesperLam)*Constraints.ply_t;        % include weight penalty

    if abs(fitness1-fitness)>1e-12, keyboard; end
    
if ~FEASIBLE  % add penalty if not FEASIBLE
    if isnan(fitness)
        keyboard
    else
    fitness = fitness * 4 ;
    end
end

output.SS          = SS;
output.DropIndexes = DropIndexes;
output.FEASIBLE    = FEASIBLE;

end


