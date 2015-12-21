% =====                                                              ==== 
%                 Fitness evaluation function used by the GA              
%                                                                                                        
% [fitness,output] = Eval_Fitness (Individual,Objectives,Constraints,NpatchVar,NthetaVar,AllowedNplies,LamType)
%
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



function [fitness,output] = Eval_Fitness (Individual,Objectives,Constraints,NpatchVar,NthetaVar,AllowedNplies,LamType,Fixed)


ConstraintVector = Constraints.Vector;

FEASIBLE  = true; 
Nlam      = size(Objectives.Table,1)-1;
LamNumber = int8(zeros(Nlam,1));
for j=1:Nlam
    LamNumber(j) = Objectives.Table{j+1,1}; 
end

[SortedLamNumber,GuideAngles,ShuffleLoc,DropIndexes] = Convert_Genotype(Individual,LamNumber,Constraints,NpatchVar,NthetaVar,AllowedNplies,Fixed);
NGuidePlies = length(GuideAngles);
NDropPlies  = cellfun(@length,DropIndexes,'UniformOutput', true); % NUmber of dropped between each patch
Ndrop       = length(NDropPlies);                                       % Total number of drops between patches



%% --- Convert discrete values into Angles in degree
if ~ConstraintVector(1)     % if not Damtol
    GuideAngles = -90 + GuideAngles*Constraints.DeltaAngle;
else
    GuideAngles(2:end) = -90 + GuideAngles(2:end)*Constraints.DeltaAngle;
end
if ~isempty(find(abs(GuideAngles)>90,1)),
    error('Angle greater than 90 detected');
end



%% --- check / enforce constraints for ply Guide 
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

if ConstraintVector(2)                                                      % 10% rule (0,90,+-45)
    [GuideAngles] = Enforce_10PercentRule(GuideAngles);
end

if FEASIBLE
    [FEASIBLE] = Check_Feasibility(ConstraintVector,GuideAngles,ShuffleLoc,cell2mat(DropIndexes),NGuidePlies,NDropPlies,LamType);
end



%% fitness calculation
% pre-allocation
SS  = cell(1,Ndrop+1);

if strcmp(Objectives.Type,'LP')
    LP = zeros(12,Ndrop+1);
end

if strcmp(Objectives.Type,'ABD')
    A = cell(Ndrop+1,1);
    B = cell(Ndrop+1,1);
    D = cell(Ndrop+1,1);
end


% Guide + After Drops
for iDrop = 0 : Ndrop
    
    if iDrop == 0
        DropsLoc = [];
    else
        DropsLoc = [DropsLoc DropIndexes{iDrop}];
        DropsLoc = unique(DropsLoc);
        DropsLoc(DropsLoc>NGuidePlies) = [];                                    % remove Infeasible Drops (only for variable Nply)
    end
    
    FiberAngles = Convert_dvAngles2FiberAngles(GuideAngles,DropsLoc,ShuffleLoc,LamType);
    SS{iDrop + 1}   = FiberAngles;
    
    if strcmp(Objectives.Type,'LP')
        LP(:,iDrop + 1) = Convert_SS2LP(FiberAngles);          % evaluate lamination parameters for the droped laminates
    end
    if strcmp(Objectives.Type,'ABD')
        [A{iDrop + 1},B{iDrop + 1},D{iDrop + 1}] = Convert_SS2ABD (Objectives.mat(1),Objectives.mat(2),Objectives.mat(4),Objectives.mat(3),Constraints.ply_t,FiberAngles,true);
    end
end


% revert both sorting at once
[~,RevertSort]    = sort(SortedLamNumber);
SS = SS(RevertSort);


% --- check individual ply continuity (only if structure geometry is given)
if isfield(Constraints,'PatchConnectivity') 
    NGeoConstraints = CheckContinuity(SS,Constraints.PatchConnectivity);
    output.NGeoConstraints = NGeoConstraints;
else
    NGeoConstraints = 0;
end
% ---


if strcmp(Objectives.Type,'LP')
    LP        = LP(:,RevertSort);
    if ~Objectives.UserFct
        fitness          = Objectives.FitnessFct(LP);                       % Default Fitness Function (Do not Change)
    else
        [fitness,output] = Objectives.FitnessFct(LP);                       % User Fitness Function Calls
    end
    output.LP = LP;
end

if strcmp(Objectives.Type,'ABD')
    A = A(RevertSort);
    B = B(RevertSort);
    D = D(RevertSort);
    
    if ~Objectives.UserFct
        fitness = Objectives.FitnessFct(A,B,D);                             % Default Fitness Function (Do not Change)
    else
        [fitness,output] = Objectives.FitnessFct(A,B,D);                    % User Fitness Function Calls
    end
    
    output.A  = A;
    output.B  = B;
    output.D  = D;
end

if strcmp(Objectives.Type,'SS')
    [fitness,output] = Objectives.FitnessFct(SS);                           % User Fitness Function Calls
end

fitness = fitness * (1+NGeoConstraints);

if ~FEASIBLE  % add penalty if not FEASIBLE
    if isnan(fitness) || isinf(fitness) || ~isreal(fitness)
        error('Non appropriate Fitness (i.e. NaN,inf or complex) has been detected')
    end
    fitness = fitness * 4;
end



output.SS          = SS;
output.DropIndexes = DropIndexes;
output.FEASIBLE    = FEASIBLE;

