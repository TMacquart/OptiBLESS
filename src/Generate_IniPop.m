% =====                                                                ====
%                         Creates Initial Population
%
% Recommended to use this function rather than the Random GA inipop. It
% will help (not ensure) the generation of feasible solution
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


function [IniPop] = Generate_IniPop (NStruct,Npop,Constraints,Objectives,AllowedNplies,LamType,BCs,fct_handle)


% Initialisation
IniPop = zeros(Npop,NStruct.Nvar);
display(' Creating IniPop ... ' )
DeltaAngle    = Constraints.DeltaAngle;
ipop          = 1 ;
NTried        = 0;



% Loop until IniPop is complete
while ipop < Npop + 1
    
    %%  Variable number of plies
    NpliesPerLam = zeros(length(AllowedNplies),1);
    for iply = 1:length(AllowedNplies)
        NpliesPerLam(iply) = AllowedNplies{iply}(randi([1 length(AllowedNplies{iply})]));
    end
    [NpliesPerLam,SortIndex] = sort(NpliesPerLam,'descend');
    NPliesGuide = max(NpliesPerLam);
    
    
    %% Drop Location
    
    NStruct0   = AttributeNply(NpliesPerLam,Constraints,LamType);
    NdropPlies = NStruct0.NdropVar;
    
    if NdropPlies>0
        
        % create a symbolic SSsequence including BalancedLoc, 10% rule and
        % PlyDrops becaue they all affect internal continuity
        
        Thetas_Symbolic    = (1:NStruct.NthetaVar)';
        
        FEASIBLE = false;
        
        while  ~FEASIBLE
            FEASIBLE = true;
            
            DropsIndexes       = randperm(BCs.UB.PlyDrop(1)-BCs.LB.PlyDrop(1) , NdropPlies) + BCs.LB.PlyDrop(1);
            
            if NStruct.NbalVar==0
                BalancedLoc = [];
            else
                BalancedLoc = (randperm(BCs.UB.BalancedLoc(1)-BCs.LB.BalancedLoc(1) , NStruct.NbalVar) + BCs.LB.BalancedLoc(1))';
            end
            
            if NStruct.N10percentVar ==0
                TenPercentLoc = [];
            else
                TenPercentLoc = (randperm(BCs.UB.BalancedLoc(1)-BCs.LB.BalancedLoc(1) , NStruct.N10percentVar) + BCs.LB.BalancedLoc(1))';
            end
            
            if NStruct.NMidPlane == 0
                Thetas_Mid = [];
            else
                Thetas_Mid = (randi([BCs.LB.MidPlane(1) BCs.UB.MidPlane(1)] ,1 , NStruct.NMidPlane))';
            end
            
            SSTable_Symbolic = NewComputeSSTable(Thetas_Symbolic,DropsIndexes,BalancedLoc,TenPercentLoc,Thetas_Mid,LamType,Constraints);
            
            if Constraints.Vector(6)
                FEASIBLE = Check_InternalContinuity(SSTable_Symbolic);
            end
            
        end
        
        %% Add Fiber Angles
        
        FEASIBLE = false;
        TriedSol = 1;
        while ~FEASIBLE
            Thetas = (randperm(BCs.UB.Thetas(1)-BCs.LB.Thetas(1) , NStruct.NthetaVar) + BCs.LB.Thetas(1))';
            
            % repaired solution - change NpliesPerLam to match SS values (avoid infeasibility)
            if 1
                for j = 1:size(SSTable_Symbolic,1)
                    NplySS(j) = length(find(~cellfun(@isempty,SSTable_Symbolic(j,:))));
                end
                for j = 1:length(NpliesPerLam(:,1))
                    if isempty(find(NpliesPerLam(j,1)==NplySS,1))
                        [~,NIndex] = min(abs(NplySS-NpliesPerLam(j,1)));
                        NpliesPerLam(j,1) = NplySS(NIndex);
                    end
                end
            end
            %
            
            SSTable = NewComputeSSTable(Thetas,DropsIndexes,BalancedLoc,TenPercentLoc,Thetas_Mid,LamType,Constraints);
            
            FEASIBLE = Check_Feasibility(Constraints,NpliesPerLam,SSTable);
            
            TriedSol = TriedSol +1;
            if TriedSol==2000,   break;  end     % no-feasible solution, go back one more level
           
        end
        
        
        %% Final check and Add to population
        if FEASIBLE
            keyboard

            % keyboard
            
            SSTableLocal = SSTable
            
            Individual = [%NpliesPerLam(SortIndex);  % no Nplies if fixed variables
                          Thetas;
                          BalancedLoc;
                          TenPercentLoc;
                          %Thetas_Mid;
                          0; % add thetamid if zero
                          0; % add thetamid if zero
                          DropsIndexes ];
                          
                           
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
end

display(' IniPop Created ... ' )