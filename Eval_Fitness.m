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

function [fitness,output] = Eval_Fitness (Individual,Objectives,Constraints,Nsec,AllowedNplies,SortedTable,LamType)


FEASIBLE  = true;  
LamNumber = cell2mat(SortedTable(2:end,1)); % after 1st sorting

if Constraints.NRange == 1
    NpliesperLam = cell2mat(AllowedNplies);
else
    NpliesperLam = Individual(1:Nsec);
end

% sort ply order to have Guide First
if Constraints.ORDERED
    NpliesperLam = sort(NpliesperLam,'descend');                            % repair to ensure thickness is ordered
    SortIndex    = 1:length(NpliesperLam);
    if find(diff(NpliesperLam)>0,1)
        FEASIBLE = false;
    end
else
    [NpliesperLam,SortIndex] = sort(NpliesperLam,'descend');
end

SortedLamNumber = LamNumber(SortIndex);                                         % after 2nd sorting
NGuidePlies    = max(NpliesperLam);                                           % number of plies in the guide laminate (take half for Sym.)
NDropPlies     = abs(diff(NpliesperLam));                                     % number of ply drops
GuideAngles    = Individual(Nsec + [1:NGuidePlies]);
GuideAngles    = GuideAngles(~isnan(GuideAngles));



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
    R = [-1 1];
    GuideAngles(1) = 45*R(GuideAngles(1));
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
    
    % Guide
    FiberAngles = Convert_dvAngles2FiberAngles(GuideAngles,LamType);
    
    if strcmp(Objectives.Type,'LP')
        LP(:,1)  = Convert_SS2LP(FiberAngles);            % evaluate lamination parameters for the guide
    end
    if strcmp(Objectives.Type,'ABD')
        [A{1},B{1},D{1}] = Convert_SS2ABD (Objectives.mat(1),Objectives.mat(2),Objectives.mat(4),Objectives.mat(3),Constraints.ply_t,FiberAngles,true);
    end
    SS{1} = FiberAngles;
    
    % After Drops
    for iDrop = 1 : Ndrop
        index = iDrop + 1;
        
        ply_angles = GuideAngles;
        
        DropsLoc = unique(cell2mat(DropIndexes(1:iDrop)));
        if max(DropsLoc)>NGuidePlies
            DropsLoc(DropsLoc>NGuidePlies) = [];
        end
        
        ply_angles(DropsLoc(DropsLoc<length(ply_angles))) = [];  % drop plies
        
        FiberAngles      = Convert_dvAngles2FiberAngles(ply_angles,LamType);
        SS{index}        = FiberAngles;
        
        if strcmp(Objectives.Type,'LP')
            LP(:,index) = Convert_SS2LP(FiberAngles);          % evaluate lamination parameters for the droped laminates
        end
        if strcmp(Objectives.Type,'ABD')
            [A{index},B{index},D{index}] = Convert_SS2ABD (Objectives.mat(1),Objectives.mat(2),Objectives.mat(4),Objectives.mat(3),Constraints.ply_t,FiberAngles,true);
        end
    end
end


% revert both sorting at once
[~,RevertSort]    = sort(SortedLamNumber);          
RevertedLamNumber = SortedLamNumber(RevertSort);
SS = SS(RevertSort);

if strcmp(Objectives.Type,'LP')
    LP = LP(:,RevertSort);
    fitness = Objectives.FitnessFct(LP);
end
if strcmp(Objectives.Type,'ABD')
    A = A(RevertSort);
    B = B(RevertSort);
    D = D(RevertSort);
    fitness = Objectives.FitnessFct(A,B,D);
end
if strcmp(Objectives.Type,'SS') 
    fitness = Objectives.FitnessFct(SS);
end


if ~FEASIBLE  % add penalty if not FEASIBLE
    if isnan(fitness) || isinf(fitness) || ~isreal(fitness)
        error('Non appropriate Fitness (i.e. NaN,inf or complex) has been detected')
    end
    fitness = fitness * 4 ;
end



if strcmp(Objectives.Type,'ABD')
    output.A  = A;
    output.B  = B;
    output.D  = D;
end
if strcmp(Objectives.Type,'LP')
    output.LP = LP;
end

output.LamIndex    = RevertedLamNumber;
output.SS          = SS;
output.DropIndexes = DropIndexes;
output.FEASIBLE    = FEASIBLE;

end


