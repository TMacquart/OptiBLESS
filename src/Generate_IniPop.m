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


function [IniPop] = Generate_IniPop (nvar,Npop,NpatchVar,NthetaVar,NdropVar,Constraints,Objectives,AllowedNplies,LamType,fct_handle)


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
    
    %% Drop Location
    NdropPlies = NdropVar;
    if NdropPlies>0
        if ConstraintVector(7) % covering
            Offset = 1;
        else
            Offset = 0;
        end
        
        if strcmp(LamType,'Generic') || strcmp(LamType,'Sym')
            if strcmp(LamType,'Generic')
                Coeff = 1;
            else
                Coeff = 2;
            end
            
            DropsIndexes = randperm(NthetaVar-Offset,NdropPlies)+Offset;
            
            if ConstraintVector(6) % Internal continuity
                FEASIBLE = false;
                while ~FEASIBLE
                    FEASIBLE = true;
                    DropsLoc = sort(DropsIndexes);
                    for i = 1 : length(DropsLoc)-3
                        deltaLoc = diff(DropsLoc(i:i+3));
                        if sum(abs(deltaLoc))==3
                            FEASIBLE = false;
                            DropsIndexes = randperm(NthetaVar-Offset,NdropPlies)+Offset;
                            break
                        end
                    end
                end
            end
        end
        
        
        if strcmp(LamType,'Balanced') || strcmp(LamType,'Balanced_Sym')
            if strcmp(LamType,'Balanced')
                Coeff = 2;
            else
                Coeff = 4;
            end
            
            FEASIBLE = false;
            while ~FEASIBLE
                FEASIBLE = true;
                ShuffleLoc   = randperm(NthetaVar*2-1,NthetaVar)+1; % First ply not allowed for simplicity
                DropsIndexes = randperm(NthetaVar-Offset,NdropPlies)+Offset;
                
                if ConstraintVector(6) % check internal continuity
                    PatchNPly = sort(unique(NpliesPerLam*Coeff),'descend');
                    SymbolicGuideAngles =  [1:NthetaVar];
                    
                    % reconstuct 1st and Last
                    SymbolicSS     = num2cell(Convert_dvAngles2FiberAngles(SymbolicGuideAngles,[],ShuffleLoc,LamType)');
                    SymbolicSS_End = Convert_dvAngles2FiberAngles(SymbolicGuideAngles,DropsIndexes,ShuffleLoc,LamType)';
                    
                    for i=1:NthetaVar*Coeff
                        if ~isempty(find(SymbolicSS{1,i}==SymbolicSS_End,1))
                            SymbolicSS{2,i} = SymbolicSS{1,i};
                        end
                    end
                    
                    
                    SSDrops = sort(find(cellfun(@isempty,SymbolicSS(2,:))));
                    for i = 1 : length(SSDrops)-3
                        deltaLoc = diff(SSDrops(i:i+3));
                        if sum(abs(deltaLoc))==3
                            FEASIBLE = false;
                            break
                        end
                    end
                end
            end
        end

    else
        DropsIndexes = [];
    end
    
    
    %% Add Fibre Angles
    if 1 % reconstruct all
        PatchNPly = sort(unique(NpliesPerLam*Coeff),'descend');
        SymbolicGuideAngles =  [1:NthetaVar];
        
        index = 1;
        SymbolicSSTable = cell(length(PatchNPly),max(PatchNPly));
        for j=0:length(DropsIndexes)
            if ~isempty(find(((NthetaVar-j)*Coeff-PatchNPly)==0,1))
                
                SymbolicFiberAngles = Convert_dvAngles2FiberAngles(SymbolicGuideAngles,DropsIndexes(1:j),ShuffleLoc,LamType)';
                if index == 1
                    SymbolicSSTable(index,:) = num2cell(SymbolicFiberAngles);
                else
                    for i=1:max(PatchNPly)
                        if ~isempty(find(SymbolicSSTable{1,i}==SymbolicFiberAngles,1))
                            SymbolicSSTable{index,i} = SymbolicSSTable{1,i};
                        end
                    end
                end
                
                index = index + 1;
            end
        end
    end
            
    GuideNumber = 1:NthetaVar;
    GuideAssign = zeros(1,NthetaVar);
    GuideAngles = nan*zeros(1,NthetaVar);
    
    if ConstraintVector(1) % Damtol (and covering assumed)
        A = [-1 1];
        GuideAngles(1) = 45*A(ceil(2*rand)); % 1st ply is +- 45
        if strcmp(LamType,'Generic')
            GuideAngles(end) = 45*A(ceil(2*rand)); % Last ply is +- 45
        end
        GuideAssign(1) = 1;
    end
    
%     keyboard
    for j=length(SymbolicSSTable(:,1)):-1:1
        
        plyNumber = cell2mat(SymbolicSSTable(j,:));
        plyAngle  = nan*ones(1,length(plyNumber));
        UniqueplyNumber  = unique(abs(plyNumber));
        UnAssignedNumber = UniqueplyNumber(~ismember(UniqueplyNumber,find(GuideAssign)));
        
        for i=1:length(UniqueplyNumber) %replace by assigned angles
            if GuideAssign(UniqueplyNumber(i)) == 1
                Indexes           = find(UniqueplyNumber(i)==abs(plyNumber));
                plyAngle(Indexes) = GuideAngles(UniqueplyNumber(i)).*sign(plyNumber(Indexes));
            end
        end

        TriedSol = 0;
        FEASIBLE = false;
        clear FiberAngle
        while ~FEASIBLE && TriedSol<500
            FEASIBLE = true;
            
            FiberAngle = randi([0 length(0:DeltaAngle:180)-1],1,length(UnAssignedNumber))*DeltaAngle-90;
            
              if ConstraintVector(2) && j==length(SymbolicSSTable(:,1))
                %keyboard % need to change fibre angles accordingly
                FiberAngle = Enforce_10PercentRule( [GuideAngles(~isnan(GuideAngles)) FiberAngle]); % only enforce for the thinnest and let random for the rest
                if ConstraintVector(1)
                    FiberAngle(1) = [];
                end
              end
            
            for i=1:length(UnAssignedNumber)
                Indexes           = find(UnAssignedNumber(i)==abs(plyNumber));
                plyAngle(Indexes) = FiberAngle(i).*sign(plyNumber(Indexes));
            end
            
            DetlaAngle = ComputeDeltaAngle(plyAngle);

            if ConstraintVector(3) && ~isempty(find(DetlaAngle>45,1)) % disorientation
                FEASIBLE = false;
            end
            
            if FEASIBLE && ConstraintVector(4) % contiguity
                Contiguity=0;
                for i=1:length(DetlaAngle)
                    if DetlaAngle(i)==0
                        Contiguity = Contiguity + 1;
                    else
                        Contiguity = 0;
                    end
                    if Contiguity >= Constraints.Contiguity
                        FEASIBLE = false;
                    end
                end
            end
            
            if FEASIBLE && ConstraintVector(2) % 10% rule
                FEASIBLE = Check_10PercentRule(plyAngle);
            end
            

            
            TriedSol = TriedSol +1;
        end
        
        if TriedSol==500
            % no-feasible solution, go back one more letter
            break
        end
        
        GuideAssign(UnAssignedNumber) = 1;
        GuideAngles(UnAssignedNumber) = FiberAngle;
    end
    
    if FEASIBLE && ~isempty(find(isnan(GuideAngles),1))
        FEASIBLE = false;
    end
    
    
    
    %% OLD APPROACH
    if 0
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
                    A = round([(-45 + 40*rand) (5 + 40*rand)]/DeltaAngle)*DeltaAngle;
                    AddedAngle = GuideAngles(iAngle-1) + A(ceil(2*rand));
                else
                    AddedAngle = GuideAngles(iAngle-1) + round((-45 + 90*rand)/DeltaAngle)*DeltaAngle;
                end
            else
                if ConstraintVector(4)
                    A = round([(-90 + 85*rand) (5 + 85*rand)]/DeltaAngle)*DeltaAngle;
                    AddedAngle = GuideAngles(iAngle-1) + A(ceil(2*rand));
                else
                    AddedAngle = randi([0 length(0:DeltaAngle:180)-1],1,1)*DeltaAngle-90; % no constraint - full range (use randi for uniform PDF)
                end
            end
            if AddedAngle>90,    AddedAngle = AddedAngle - 180;   end
            if AddedAngle<-90,   AddedAngle = AddedAngle + 180;   end
            
            GuideAngles(iAngle) = AddedAngle;
        end                                     % Stack plies according to constraints activated
        
        %     keyboard
        
        %% 10% rule Framework
        if ConstraintVector(2)
            [GuideAngles] = Enforce_10PercentRule(GuideAngles);
        end
        
        
        %%
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
            [FEASIBLE] = Check_Feasibility(ConstraintVector,GuideAngles(1:NPliesGuide),ShuffleLoc(1:NPliesGuide),DropsIndexes(1:NGuideDropPlies),NPliesGuide,NGuideDropPlies,LamType,Constraints.Contiguity);
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
        
        
    end
    %% End OLD APPROACH
    
    %% Final check and Add to population
    if FEASIBLE
        keyboard
        if 0 % check purpose only
            SSTable = ComputeSSTable(unique(NpliesPerLam),GuideAngles,DropsIndexes,ShuffleLoc,LamType); 
        end
        
        GuideAnglesInt = round((90+GuideAngles)/DeltaAngle);
        if ConstraintVector(1)
            if GuideAngles(1) == -45,  GuideAnglesInt(1)=1; end
            if GuideAngles(1) ==  45,  GuideAnglesInt(1)=2; end
            if abs(GuideAngles(1)) ~= 45,  keyboard; end % should not happen
        end
        
        NpliesPerLam = NpliesPerLam(SortIndex);
        if Constraints.Balanced
            Individual = [NpliesPerLam(NpatchVar) GuideAnglesInt ShuffleLoc DropsIndexes];
        else
            Individual = [NpliesPerLam(NpatchVar) GuideAnglesInt DropsIndexes];
        end
        
        % keyboard
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
    NTried = NTried + 1
    if  NTried == 1000000 && ipop<2
        error('Hard Constrained Problem. Difficulties Generating IniPop')
    end
    
end

display(' IniPop Created ... ' )
end