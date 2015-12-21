% =====                                                                ==== 
%                         Creates Initial Population                      
%                                                                         
% Recommended to use this function rather than the Random GA inipop.                            
%
% [IniPop] = Generate_IniPop (nvar,Npop,NpatchVar,NthetaVar,NdropVar,Constraints,AllowedNplies,LamType)
%
%
% nvar          : total # a design variables 
% Npop          : Size of the initial population
% NpatchVar     : Number of patches with variable number of plies
% NthetaVar     : Number of fibre angles used to describe the guide laminate
% NdropVar      : Vector containing the # of plies to drop per section
% Constraints   : Input Structure containing manufacturing constraints 
% AllowedNplies : Number of plies allowed for each laminate
% LamType       : Laminate type (e.g. balanced, symmetric, ...)
%
%
% IniPop        : Initial population matrix
% =====                                                                ====

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


function [IniPop] = Generate_IniPop (nvar,Npop,NpatchVar,NthetaVar,NdropVar,Constraints,Objectives,AllowedNplies,LamType,Fixed,fct_handle)


% Initialisation
IniPop = zeros(Npop,nvar);
display(' Creating IniPop ... ' )
ConstraintVector = Constraints.Vector;
DeltaAngle       = Constraints.DeltaAngle;
ipop          = 1 ;
NTried        = 0;



% Loop until IniPop is complete
while ipop < Npop + 1
    
    %%  Variable number of plies
    NpliesPerLam = zeros(1,length(AllowedNplies));
    for iply = 1:length(AllowedNplies)
        NpliesPerLam(iply) = AllowedNplies{iply}(randi([1 length(AllowedNplies{iply})]));
    end
    [NpliesPerLam,SortIndex] = sort(NpliesPerLam,'descend');
    NPliesGuide = max(NpliesPerLam);
    
    GuideAngles = zeros(1,NthetaVar); % some non-coded genes are included
    

    %% Damtol
    if ConstraintVector(1)
        A = [-1 1];
        GuideAngles(1) = 45*A(ceil(2*rand)); % 1st ply is +- 45
        if strcmp(LamType,'Generic')
            GuideAngles(end) = 45*A(ceil(2*rand)); % Last ply is +- 45
        end
    else
        GuideAngles(1)   = randi([0 length(0:DeltaAngle:180)-1],1,1)*DeltaAngle -90;
        GuideAngles(end) = randi([0 length(0:DeltaAngle:180)-1],1,1)*DeltaAngle -90;
    end
    
    
    %% Disorientation and Contiguity
    for iAngle = 2 : NthetaVar - 1 % NPliesGuide-1
        if ConstraintVector(3)
            if ConstraintVector(4)
                A = [(-45 + 40*rand) (5 + 40*rand)];
                AddedAngle = GuideAngles(iAngle-1) + A(ceil(2*rand));
            else
                AddedAngle = GuideAngles(iAngle-1) + (-45 + 90*rand);
            end
        else
            if ConstraintVector(4)
                A = [(-90 + 85*rand) (5 + 85*rand)];
                AddedAngle = GuideAngles(iAngle-1) + A(ceil(2*rand));
            else
                AddedAngle = randi([0 length(0:DeltaAngle:180)-1],1,1)*DeltaAngle-90; % no constraint - full range (use randi for uniform PDF)
            end
        end
        if AddedAngle>90,    AddedAngle = AddedAngle - 180;   end
        if AddedAngle<-90,   AddedAngle = AddedAngle + 180;   end
        
        GuideAngles(iAngle) = AddedAngle;
    end                                     % Stack plies according to constraints activated
    

    %% 10% rule Framework
    if ConstraintVector(2)
%         keyboard % might be to res
%         RuleAngles     = [45 90 -45 0]; 
%         FreeLoc        = 1:(Fixed.Ntheta+NthetaVar);
%         [~,Index]      = ismember(Fixed.thetaLoc,1:8);
%         FreeLoc(Index) = [];
%         
%         AddAngles = RuleAngles(1+rem(0:Fixed.Ntheta-1,4));
%         NewGuideAngles = zeros(1,(Fixed.Ntheta+NthetaVar));
%         NewGuideAngles(Fixed.thetaLoc) = AddAngles;
%         NewGuideAngles(FreeLoc) = GuideAngles;
%         
%         GuideAngles = NewGuideAngles;
%         
%         ShuffleLoc = zeros(1,Fixed.Ntheta+NthetaVar);
%         ShuffleLoc(Fixed.thetaLoc) = Fixed.thetaLoc +  (Fixed.Ntheta+NthetaVar)
        [GuideAngles] = Enforce_10PercentRule(GuideAngles);
    end
    
    
        %% 
        
%         keyboard
        if sum(abs(GuideAngles-round((GuideAngles)/DeltaAngle)*DeltaAngle))~= 0
            keyboard
        end
    GuideAngles = round((GuideAngles)/DeltaAngle)*DeltaAngle; % round up Angles to discrete value (to remove?)
    
    
    ShuffleLoc = randperm(NthetaVar*2,NthetaVar); % not used if not balanced
    if ~isempty(find(ShuffleLoc==0,1)), keyboard; end
    
    FEASIBLE = true;
    
    %% Drop Locations
    NdropPlies = NdropVar;
    if NdropPlies>0
        DropsIndexes = Generate_DropIndexes (GuideAngles,NdropPlies,ConstraintVector,LamType);
        if isempty(DropsIndexes)
            FEASIBLE = false;
        end
    else
        DropsIndexes = [];
    end
    
    %% Feasibility check
%     keyboard
    if FEASIBLE
        NGuideDropPlies = NPliesGuide-min(NpliesPerLam);
        [FEASIBLE] = Check_Feasibility(ConstraintVector,GuideAngles(1:NPliesGuide),ShuffleLoc(1:NPliesGuide),DropsIndexes(1:NGuideDropPlies),NPliesGuide,NGuideDropPlies,LamType);
    end
    
    
    %% Convert Angles in degree into the discrete corresponding values
    if ~ConstraintVector(1)
        GuideAngles = round((90+GuideAngles)/DeltaAngle);
    else
        if strcmp(LamType,'Generic')
            GuideAngles(2:end-1) = round((90+GuideAngles(2:end-1))/DeltaAngle);
        else
            GuideAngles(2:end) = round((90+GuideAngles(2:end))/DeltaAngle);
        end
    end
    
    
    if ConstraintVector(1) % make the first ply a discrete value with only 2 possible state
        if (GuideAngles(1) == -45), GuideAngles(1) = 1; end
        if (GuideAngles(1) == 45),  GuideAngles(1) = 2; end
        if strcmp(LamType,'Generic')
            if (GuideAngles(end) == -45), GuideAngles(end) = 1; end
            if (GuideAngles(end) == 45),  GuideAngles(end) = 2; end
        end
    end
    

    %% Final check and Add to population
    if FEASIBLE

        NpliesPerLam = NpliesPerLam(SortIndex);
        if Constraints.Balanced
            Individual = [NpliesPerLam(NpatchVar) GuideAngles ShuffleLoc DropsIndexes];
        else
            Individual = [NpliesPerLam(NpatchVar) GuideAngles DropsIndexes];
        end
        
%         keyboard
        [~,output] = fct_handle(Individual);
        if isfield(Constraints,'PatchConnectivity') 
            NGeoConstraints = CheckContinuity(output.SS,Constraints.PatchConnectivity);
        else
            NGeoConstraints = 0;
        end
%         display(NGeoConstraints)
        
        if output.FEASIBLE
            if Objectives.UserFct % Check based on user function = user fitness fct
                if output.NViolatedConst<=1
                    display(ipop)
                    IniPop(ipop,:) = Individual;
                    ipop   = ipop + 1;
                    NTried = 0;
                end
            else
                IniPop(ipop,:) = Individual;
                ipop = ipop + 1;
            end
        end
    end
    
    
    %% Stop if too hard to create the initial populations
    NTried = NTried + 1;
    if  NTried == 10000 && ipop<2
        error('Hard Constrained Problem. Difficulties Generating IniPop')
    end
    
end

display(' IniPop Created ... ' )
end