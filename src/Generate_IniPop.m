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
% =====                                                                ====
%                         Creates Initial Population
%
% Recommended to use this function rather than the Random GA inipop. It
% will help (not ensure) the generation of feasible solution.
%
% [IniPop] = Generate_IniPop (NStruct,Npop,Constraints,Objectives,AllowedNplies,LamType,BCs,fct_handle)
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
% IniPop        : Initial population matrix (coded in genotype format)
% =====                                                                ====



function [IniPop] = Generate_IniPop (NStruct,NStructMin,GAoptions,Constraints,Objectives,LamType,BCs,fct_handle)

% Initialisation
Npop   = GAoptions.Npop;
IniPop = zeros(Npop,NStruct.Nvar);
display(' Creating IniPop ... ' )
ipop          = 1 ;
NTried        = 0;


if Constraints.Vector(7)
    if NStruct.NInsertVar > NStructMin.MinNplies*1.5 % Empirical limit for Constraints.NInternalCont = 3
        warning(['A Feasible initial Population will be difficult or even impossible to generate due to internal continuity constraints.' ...
            'You may which to increase the value of Constraints.NInternalCont or reduced the difference between Max and Min Nply'])
        
        str = input('Would you like to continue anyway? [Y/N]','s');
        if ~strcmp(str,'Y')
            error(' ---------   Code stopped by users during initial population generation  ---------')
        end
    end
%     if NStruct.NInsertVar > (NStructMin.MinNplies*3), % ultimate limit
%         error('Impossible to satisfy all constraints due to internal continuity')
%     end
end

% Loop until IniPop is complete
while ipop < Npop + 1
% parfor Iparallel = 1 : 100000
    
    %% First, we generate a feasible thin laminate ( stored into SSTable(1,:) ), we add ply afterwards
    %  Start by generating balanced and 10% rule angle locations as well as  mid plane angles.
    
    if 1
        InsertIndexes = [];
    else
        if NStruct.NInsertVar > 0
            InsertIndexes = randi([BCs.LB.InsertIndex(1)+2 BCs.UB.InsertIndex(1)],1,2);
        else
            InsertIndexes = [];
        end
    end
    
    if NStruct.NbalVar==0
        BalancedLoc = [];
    else
        BalancedLoc = randi([BCs.LB.BalancedLoc(1) BCs.UB.BalancedLoc(1)],1,NStruct.NbalVar);
    end
    
    if NStruct.N10percentVar ==0
        TenPercentLoc = [];
    else
        TenPercentLoc = zeros(1,NStruct.N10percentVar);
        for j=1:NStruct.N10percentVar
            if j == 1
                TenPercentLoc(j) = randi([BCs.LB.TenPercentLoc(j) BCs.UB.TenPercentLoc(j)]);
            else
                Finished = false;
                while ~Finished
                    temp = randi([BCs.LB.TenPercentLoc(j) BCs.UB.TenPercentLoc(j)]);
                    if isempty(find(TenPercentLoc==temp,1))
                        TenPercentLoc(j) = temp;
                        Finished = true;
                    end
                end
            end
        end
    end
    
    if NStruct.NMidPlane == 0
        Thetas_Mid = [];
    else
        if NStruct.NMidPlane == 1
            Thetas_Mid(1) = randi([BCs.LB.MidPlane(1) BCs.UB.MidPlane(1)]);
        end
        
        if NStruct.NMidPlane == 2 || NStruct.NMidPlane == 3
            Thetas_Mid(1) = randi([BCs.LB.MidPlane(1) BCs.UB.MidPlane(1)]);
            Thetas_Mid(2) = randi([BCs.LB.MidPlane(2) BCs.UB.MidPlane(2)]);
        end
    end
    
    
    % try various Guide Fiber Angles (no ply drops yet) until a feasible sol. is found
    
    FEASIBLE = false;
    Thetas   = zeros(1,NStruct.NthetaVar);
    TriedSol = 0;   % Keep track of the number of tried solutions (stop if too many)
    
    % start by generatng full random theta vector
    for j=1:NStruct.NthetaVar
        Thetas(j) = randi([BCs.LB.Thetas(j) BCs.UB.Thetas(j)]);
    end
    
    SSTable           = ComputeSSTable(Thetas,InsertIndexes,BalancedLoc,TenPercentLoc,Thetas_Mid,LamType,Constraints,NStruct,NStructMin,false);
    [FEASIBLE,output] = Check_Feasibility(Constraints,SSTable);
    
    % repair for disorientation until feasible or discard if impossible
    while ~FEASIBLE && TriedSol<100;
        TriedSol = TriedSol +1;
        
        SSTableSymbolic   = ComputeSSTable(Thetas,InsertIndexes,BalancedLoc,TenPercentLoc,Thetas_Mid,LamType,Constraints,NStruct,NStructMin,true);
        [FEASIBLE,output] = Check_Feasibility(Constraints,SSTable,[],SSTableSymbolic);
        
        if ~FEASIBLE
            if strcmp(output.ConstViolated,'Disorientation')
                
                if output.FailedDv<=NStructMin.NthetaVar
                    % Change thetas angle
                    Thetas(output.FailedDv) = randi([BCs.LB.Thetas(output.FailedDv) BCs.UB.Thetas(output.FailedDv)]);
                    
                else
                    
                    if output.FailedDv==(NStructMin.NthetaVar+1) && NStruct.NMidPlane>0
                        Thetas_Mid(1) = randi([BCs.LB.MidPlane(1) BCs.UB.MidPlane(1)]);
                    end
                    
                    if Constraints.Vector(4)
                        % change 10% rule angle
                        try
                            TenPercentIndex = unique(cell2mat(SSTableSymbolic(1,find(cell2mat(SSTableSymbolic(1,:))>NStructMin.NthetaVar))));
                            FailedIndex = find(TenPercentIndex==output.FailedDv);
                            
                            Finished = false;
                            while ~Finished
                                temp = randi([BCs.LB.TenPercentLoc(FailedIndex) BCs.UB.TenPercentLoc(FailedIndex)]);
                                if isempty(find(TenPercentLoc==temp,1))
                                    TenPercentLoc(FailedIndex) = temp;
                                    Finished = true;
                                end
                            end
                        catch
                            keyboard
                        end
                    end
                end
                SSTable = ComputeSSTable(Thetas,InsertIndexes,BalancedLoc,TenPercentLoc,Thetas_Mid,LamType,Constraints,NStruct,NStructMin,false);
                
                
            else
                break
            end
            
        end % end repair
    end % end while
    
    %% Add plies (a feasible thinnest laminate has been found (you can check it in SSTable) - now we add feasible plies)
    
    if FEASIBLE
        
        
%         keyboard
        if 1 % random
            if NStruct.NInsertVar>0
                InsertIndexes = [];
                InsertTried   = 0;
                InsertVar     = 1 + length(InsertIndexes);
                
                while InsertVar<=NStruct.NInsertVar && InsertTried<201%35                     
                    
                    Nloc = length(SSTable);
                    if Constraints.Vector(1), Nloc = floor(Nloc/2); end % symmetric
                    
                    InsertTrial = InsertIndexes;
                    InsertTrial = [InsertTrial randi([BCs.LB.InsertIndex(InsertVar)+2 Nloc+2])];

%                     InsertTrial = [InsertTrial randi([BCs.LB.InsertIndex(InsertVar)+2 BCs.UB.InsertIndex(InsertVar)])];
                   
                    SSTable = ComputeSSTable(Thetas,InsertTrial,BalancedLoc,TenPercentLoc,Thetas_Mid,LamType,Constraints,NStruct,NStructMin,false);
                    [FEASIBLE,output] = Check_Feasibility(Constraints,SSTable);
                      
%                     output
                    
                    if FEASIBLE
                        InsertIndexes = InsertTrial;
                        InsertVar     = InsertVar +1;
                        InsertTried   = 0;
                    else
                        InsertTried = InsertTried +1;
                    end
                    
                    if  rem(InsertTried,10) && Constraints.Vector(4) && (NStructMin.N10percentVar+InsertVar)<NStruct.N10percentVar 
%                         TenPercentLoc(NStructMin.N10percentVar+InsertVar) = randi([BCs.LB.TenPercentLoc(1) BCs.UB.TenPercentLoc(1)]);
                        TenPercentLoc(NStructMin.N10percentVar+InsertVar) = randi([BCs.LB.TenPercentLoc(1) Nloc+2]);
                    end
                    
                    if rem(InsertTried,15) && Constraints.Vector(2) && (NStructMin.NbalVar+InsertVar)<NStruct.NbalVar
%                         BalancedLoc(NStructMin.NbalVar+InsertVar) = randi([BCs.LB.BalancedLoc(1) BCs.UB.BalancedLoc(1)]);
                        BalancedLoc(NStructMin.NbalVar+InsertVar) = randi([BCs.LB.BalancedLoc(1) Nloc+2]);
                    end
                    
                    if  (rem(InsertTried,20) || strcmp(output.ConstViolated,'TenPercentRule')) && (NStructMin.NbalVar+InsertVar)<NStruct.NthetaVar
                        Thetas(NStructMin.NthetaVar+InsertVar) = randi([BCs.LB.Thetas(2) BCs.UB.Thetas(2)]);
                    end
                    
                    if FEASIBLE && length(SSTable(1,:))>=NStruct.MaxNplies
                        break
                    end
                end % end while
            end
        end
        
%         if ~exist('BestInsertVar','var')
%             BestInsertVar = InsertVar
%         elseif InsertVar>BestInsertVar
%             BestInsertVar = InsertVar
%             keyboard
%         end

        
%         keyboard
        if 0 % Enumeration
            if NStruct.NInsertVar>0
                InsertVar   = 1;
                InsertTried = 0;
                InsertIndexes = zeros(1,NStruct.NInsertVar);
                
                NonAllowed = cell(NStruct.NInsertVar,1);
                j=1;
                iEnumerate = 0;
                NEnumerate = NStruct.NInsertVar*5;
                while (j<=NStruct.NInsertVar) && (iEnumerate<NEnumerate)
                    iEnumerate = iEnumerate +1;
                    PossibleInsert = randperm(BCs.UB.InsertIndex(InsertVar)-1,BCs.UB.InsertIndex(InsertVar)-1)+1;
                    PossibleInsert(ismember(PossibleInsert,NonAllowed{j})) = [];
                    
                    for i=1:length(PossibleInsert)
                        
                        InsertIndexes(j)  = PossibleInsert(i);
                        SSTable           = ComputeSSTable(Thetas,InsertIndexes,BalancedLoc,TenPercentLoc,Thetas_Mid,LamType,Constraints,NStruct,NStructMin,false);
                        [FEASIBLE,output] = Check_Feasibility(Constraints,SSTable);
                        
                        if FEASIBLE
                            j=j+1;
                            iEnumerate = 0;
                            break % break out the for loop
                        end
                    end
                    
                    if ~FEASIBLE % no feasible solution with the current InsertIndexes, remove the last index
                        if (j-1)==0
                            break
                        end
                        InsertIndexes(j) = 0;
                        NonAllowed{j}    = [];
                        j=j-1;
                        NonAllowed{j}    = [NonAllowed{j} InsertIndexes(j)];
                        InsertIndexes(j) = 0;
                        InsertTried      = InsertTried +1;
                    end
                    
                    if FEASIBLE && length(SSTable(1,:))==NStruct.MaxNplies
                        break
                    end
                end % end while
            end
        end
        
    end
    
%         keyboard
    
    %% Finally, add number of plies for each patch until feasible (or discard) and Add to population
    
    
    if FEASIBLE
        
        % --- generate the Number of plies for each patch
        NvarPatch     = length(BCs.UB.Nply);            % number of patches with variable thickness
        if NvarPatch>0
            NPly_Choices  = BCs.UB.Nply-BCs.LB.Nply;
            % random number of plies
            NpliesPerLam_Coded = zeros(1,NvarPatch);
            for iply = 1:NvarPatch
                NpliesPerLam_Coded(iply) = randi([1 NPly_Choices(iply)],1,1);
            end
            
        else
            NpliesPerLam_Coded = [];
        end
        % ---
        
        
        % --- Construct the Individual
        Individual = [  NpliesPerLam_Coded, ...
            Thetas, ...
            BalancedLoc, ...
            TenPercentLoc, ...
            Thetas_Mid...
            ones(1,NStruct.NDV_NMidPlane-length(Thetas_Mid))...
            InsertIndexes zeros(1,NStruct.NInsertVar-length(InsertIndexes))];
        % ---
        
        
        % --- Internal validation for balanced solutions (can be commented for fastest inipop generation)
        if Constraints.Vector(2) % balanced
            SSTable = ComputeSSTable(Thetas,InsertIndexes,BalancedLoc,TenPercentLoc,Thetas_Mid,LamType,Constraints,NStruct,NStructMin,false);
            for irow = size(SSTable,1):-1:1
                temp = cell2mat(SSTable(irow,:));
                temp(abs(temp)==90) = [];
                if sum(temp)~=0
                    keyboard
                end
            end
        end
        % ---
        
        
        if GAoptions.IniPopFEASIBLE == 1
            % add to populations (satisfy all design guidelines)
            IniPop(ipop,:) = Individual';
            display(ipop)
            ipop = ipop + 1;
            NTried = 0;
            
        elseif Objectives.UserFct % this option is only available for user-based fitness function
            % check that it is also feasible w.r.t user defined fitness function
            
            TriedSol = 1;
            while TriedSol<100
                
                [~,output] = fct_handle(Individual);
                
                if output.NViolatedConst == 0,    break;   end
                
                %%  Increase thickness of all failed panels
                FailedLamIndex = find(output.BucklingFactor>0);             % replace BucklingFactor by the field you would like to correspond to your user function
                
                NpliesPerLam_Coded(FailedLamIndex) = NpliesPerLam_Coded(FailedLamIndex) + 1;
                BoundedIndex = NpliesPerLam_Coded'>NPly_Choices;
                NpliesPerLam_Coded(BoundedIndex) = NPly_Choices(BoundedIndex);
                Individual(1:NvarPatch) = NpliesPerLam_Coded;                    % update new Ply numbers
                
                TriedSol = TriedSol +1;
            end
            
            
            if output.FEASIBLE && output.NViolatedConst==0
                display(ipop)
                IniPop(ipop,:) = Individual';
                ipop   = ipop + 1;
                NTried = 0;
            end
        else
            error('You cannot set GAoptions.IniPopFEASIBLE = 2  if Objectives.UserFct is not set to true')
        end
    end
    
    
    
    % --- Check geometric continuity of ply over patches
    %         if FEASIBLE
    %             NGeoConstraints = 0;
    %             if isfield(Constraints,'PatchConnectivity')
    %                 NGeoConstraints = CheckContinuity(output.SS_Patch,Constraints.PatchConnectivity);
    %             end
    %             display(NGeoConstraints)
    %         end
    % ---
    
    %% Stop if too hard to create the initial populations
    NTried = NTried + 1;
%     if NTried ==1000
%         warning('Difficulties generating feasible initial guess! This may take some times')
%     end
    if  NTried == 100000 && ipop<2
        error('Hard Constrained Problem. Difficulties Generating IniPop')
    end
end




display(' IniPop Created ... ' )

