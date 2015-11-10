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
%                    Used as fitness fct by the GA to match LPs           %
%
%  Requires 8 inputs
%  ----------------
%  Individual : GA individual vector ( includes ply angles and drops locations)
%
%
%  Returns up to 3 outputs
%  -----------------------
%  fitness: is the RMS error between the objective LP (i.e. LP_obj) and the individual LP based on the LP given by IndexLp
%  LP: the lamination parameters of the corresponding individual     
%  Outputs: a structure containing stacking sequences, Drop indexes and FEASIBILITY boolean
%
%
% =====                                                              ==== %

function [fitness,LP,Outputs] = LPMatch_FitnessEval (Individual,Constraints,NGuidePlies,NDropPlies,IndexLp,LP_obj,ImportanceFactor,LamType)

FEASIBLE         = true;                                                    % individual feasibility
ConstraintVector = Constraints.Vector;
Ndrop            = length(NDropPlies);                      


GuideAngles = Individual(1:NGuidePlies);
if  ConstraintVector(5)     % if DiscreteAngles
    if ~ConstraintVector(1)     % if not Damtol
        GuideAngles = -90 + GuideAngles*Constraints.DeltaAngle;
    else
        GuideAngles(2:end) = -90 + GuideAngles(2:end)*Constraints.DeltaAngle;
    end
    if ~isempty(find(abs(GuideAngles)>90,1)), error('An ply angle outside the [-90 +90] range has been detected');   end
end




% --- organise ply drop sequences
StartIndex   = NGuidePlies;
DropsIndexes = cell(1,Ndrop);
for iDrop = 1 : Ndrop
    DropsIndexes{iDrop}  = Individual(StartIndex + [1:NDropPlies(iDrop)]);
    StartIndex           = StartIndex + NDropPlies(iDrop);
end




% ---
if ConstraintVector(1) % if Damtol, 1st ply is +- 45
    A = [-1 1];
    GuideAngles(1) = 45*A(GuideAngles(1));
end

if ConstraintVector(2) % 10% rule
    [GuideAngles] = Enforcing_10PercentRule(GuideAngles);
end

if FEASIBLE
    [FEASIBLE] = CheckFeasibility(ConstraintVector,GuideAngles,cell2mat(DropsIndexes),LamType);
end




% --- fitness calculation
SS = cell(Ndrop+1,1); % saved evaluated stacking sequences


[FiberAngles] = dvAngles2FiberAngles(GuideAngles,LamType);
LP(:,1)       = SS2LP(Constraints.ply_t,FiberAngles); 
SS{1}         = FiberAngles;


for iDrop = 1:Ndrop
    DropsLoc = unique(cell2mat(DropsIndexes(1:iDrop)));
 
    ply_angle = GuideAngles;
    ply_angle(DropsLoc) = [];  % drop plies
    
    [FiberAngles] = dvAngles2FiberAngles(ply_angle,LamType);
    LP(:,1+iDrop) = SS2LP(Constraints.ply_t,FiberAngles);               % Evaluate Droped Lam. LPs
    SS{1+iDrop} = FiberAngles;
end


% --- calculate fitness
NuniqueSec = Ndrop + 1;
localFit   = zeros(NuniqueSec,1);
for isec = 1 : NuniqueSec
    localFit(isec) = rms(LP_obj(IndexLp,isec) - LP(IndexLp,isec));
end
fitness = sum(localFit.*ImportanceFactor);

if ~FEASIBLE
   fitness = fitness * 4;
end

Outputs.DropsIndexes = DropsIndexes;
Outputs.FEASIBLE     = FEASIBLE;
Outputs.SS           = SS;

end