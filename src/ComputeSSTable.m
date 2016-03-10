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
% =====                                                              ====== 
%           Convert the coded genotype into a Stacking sequence table, 
%                       and the number of plies corresponding
%
% [SSTable,NplySS] = ComputeSSTable(ThetasCoded,DropsIndexes,BalancedLoc,TenPercentLoc,Thetas_MidCoded,LamType,Constraints,NStruct,SymbolicTable)   
% =====                                                              ====== 

function [SSTable,NplySS] = ComputeSSTable(ThetasCoded,InsertIndexes,BalancedLoc,TenPercentLoc,Thetas_MidCoded,LamType,Constraints,NStruct,NStructMin,SymbolicTable)   

%% Split Design variables

% Design variables added to the thinnest laminate patch in order to obtain
% the full stacking sequence table
ThetasCoded_Add   = ThetasCoded(NStructMin.NthetaVar+1:end);
BalancedLoc_Add   = BalancedLoc(NStructMin.NbalVar+1:end);
TenPercentLoc_Add = TenPercentLoc(NStructMin.N10percentVar+1:end);

% Design variable describing the thinnest laminate patch
ThetasCoded   = ThetasCoded(1:NStructMin.NthetaVar);
BalancedLoc   = BalancedLoc(1:NStructMin.NbalVar);
TenPercentLoc = TenPercentLoc(1:NStructMin.N10percentVar);

% keyboard

%% Calculate fibre angles values
if SymbolicTable
    % only used for inipop generation and repair
    Thetas = 1:length(ThetasCoded);
    
else
    % convert coded Theta values into fibre angles
    Thetas = ThetasCoded;
    
    if ~Constraints.Vector(3)  
        % if not Damtol
        Thetas = -90 + Thetas*Constraints.DeltaAngle;
    else
        % first angle is +- 45
        R = [-1 1];
        Thetas(1) = 45*R(Thetas(1));
        
        if strcmp(LamType,'Generic')
            % last angle is also +- 45
            Thetas(2:end-1) = -90 + Thetas(2:end-1)*Constraints.DeltaAngle;
            Thetas(end) = 45*R(Thetas(end));
        else
            Thetas(2:end) = -90 + Thetas(2:end)*Constraints.DeltaAngle;
        end
    end
    
end


%% Insert Balanced Angles (if any)
% For Balanced Lam. we reconstruct the SS by inserting pairs ply angle (+-)

% ThetasBalanced = [Thetas;                 % Angles
%                  [1:length(Thetas)]];     % Initial Position in the thinnest lam. patch

ThetasBalanced = [Thetas; [1:length(Thetas)]; [1:length(Thetas)]];

if strcmp(LamType,'Balanced_Sym') || strcmp(LamType,'Balanced')
%     BalancedAngles = [-Thetas' BalancedLoc'];
    BalancedAngles = [-Thetas' BalancedLoc' [1:length(BalancedLoc)]'];
    BalancedAngles = sortrows(BalancedAngles,2)';
    
    for j = 1:length(BalancedAngles(1,:))
        if BalancedAngles(2,j)>size(ThetasBalanced,2)
            ThetasBalanced = [ThetasBalanced BalancedAngles(:,j)];
        else
            ThetasBalanced = [ThetasBalanced(:,1:BalancedAngles(2,j)-1) BalancedAngles(:,j)  ThetasBalanced(:,BalancedAngles(2,j):end)];
        end
    end
end



%% Insert 10% rule fibre angles (if any)
% For 10% angles need to reconstruct SS by inserting ply angle
if Constraints.Vector(4)
    
    NAngles = length(TenPercentLoc);
    TenPercentAngles = repmat([0 45 90 -45],1,ceil(NAngles/4));
    TenPercentAngles = TenPercentAngles(1:NAngles);
    TenPercentAngles = [TenPercentAngles' TenPercentLoc'];
    
    index = NStruct.NthetaVar;
    TenPercentAngles(TenPercentAngles(:,1)==0,3)   = index + [1:NAngles/4];
    TenPercentAngles(TenPercentAngles(:,1)==90,3)  = index + [NAngles/4+1:NAngles/2];
    TenPercentAngles(TenPercentAngles(:,1)==-45,3) = index + NAngles/2 + [1:NAngles/4];
    TenPercentAngles(TenPercentAngles(:,1)==+45,3) = index + NAngles/2 + [1:NAngles/4];
    
    TenPercentAngles = sortrows(TenPercentAngles,2)';


%     keyboard
    if ~SymbolicTable
        for j = 1:NAngles
            if TenPercentAngles(2,j)>size(ThetasBalanced,2)
                ThetasBalanced = [ThetasBalanced TenPercentAngles(:,j)];
            else
                ThetasBalanced = [ThetasBalanced(:,1:TenPercentAngles(2,j)-1) TenPercentAngles(:,j)  ThetasBalanced(:,TenPercentAngles(2,j):end)];
            end
        end
    else
        for j = 1:NAngles
            if TenPercentAngles(2,j)>size(ThetasBalanced,2)
                ThetasBalanced = [ThetasBalanced  ones(3,1)*TenPercentAngles(3,j) ];
            else
                ThetasBalanced = [ThetasBalanced(:,1:TenPercentAngles(2,j)-1)  ones(3,1)*TenPercentAngles(3,j)  ThetasBalanced(:,TenPercentAngles(2,j):end)];
            end
        end
    end
end

% keyboard


%% Add Mid plane angles
if ~isempty(Thetas_MidCoded)
 
    MidAngle = true;
    
    LastIndex(1,1) = max(ThetasBalanced(2,:));
    LastIndex(2,1) = max(ThetasBalanced(3,:));
    
    if strcmp(LamType,'Generic')
        Thetas_Mid = (-90 + Thetas_MidCoded*Constraints.DeltaAngle)';
%         Thetas_Mid = [Thetas_Mid; zeros(1,length(Thetas_Mid))];
        Thetas_Mid = [Thetas_Mid; zeros(1,length(Thetas_Mid)); [1:length(Thetas_Mid(1,:))] + LastIndex(2,1)];
    end
    
    if strcmp(LamType,'Balanced')
        if NStruct.NMidPlane >= NStruct.NDV_NMidPlane
            keyboard % should never happen, internal check
        end
        
        R = [0 90];
        Thetas_Mid = R(Thetas_MidCoded);
        
%         Thetas_Mid = [Thetas_Mid; zeros(1,length(Thetas_Mid))];
        Thetas_Mid = [Thetas_Mid; zeros(1,length(Thetas_Mid)); [1:length(Thetas_Mid(1,:))] + LastIndex(2,1)];
    end
    
    if strcmp(LamType,'Sym')
        if NStruct.NMidPlane >= NStruct.NDV_NMidPlane
            keyboard % should never happen, internal check
        end
        Thetas_Mid = -90 + Thetas_MidCoded(1,1)*ones(1,length(Thetas_MidCoded))*Constraints.DeltaAngle;
        
%         Thetas_Mid = [Thetas_Mid; zeros(1,length(Thetas_Mid))];
        Thetas_Mid = [Thetas_Mid; zeros(1,length(Thetas_Mid)); ones(1,length(Thetas_Mid(1,:)))*(1+LastIndex(2,1))];
    end
    
    if strcmp(LamType,'Balanced_Sym')
        R = [0 90];
        
        if NStruct.NMidPlane==1 || NStruct.NMidPlane==2
            Thetas_Mid = R(Thetas_MidCoded(1))*ones(1,length(Thetas_MidCoded));
            
            Thetas_Mid = [Thetas_Mid; zeros(1,length(Thetas_Mid)); ones(1,length(Thetas_Mid(1,:)))*(1+LastIndex(2,1))];
%             Thetas_Mid = [Thetas_Mid; zeros(1,length(Thetas_Mid))];
        end
      
        if NStruct.NMidPlane==3
            Thetas_Mid0  = R(Thetas_MidCoded(1));
            Thetas_Mid13 = R(Thetas_MidCoded(2))*ones(1,2);
            
            Thetas_Mid = [Thetas_Mid13(1) Thetas_Mid0 Thetas_Mid13(2)];
            
            
%             Thetas_Mid = [Thetas_Mid; zeros(1,3)];
            Thetas_Mid = [Thetas_Mid; zeros(1,3); [LastIndex(2,1)+1 LastIndex(2,1)+2 LastIndex(2,1)+1]];
        end

    end
    if ~SymbolicTable
        ThetasBalanced = [ThetasBalanced Thetas_Mid];
    else
        ThetasBalanced = [ThetasBalanced [Thetas_Mid(3,:);Thetas_Mid(2:3,:)]];
    end
else
    MidAngle = false;
end




%% Build the stacking sequence Table by inserting plies
% ThetasBalanced = [ThetasBalanced; [1:size(ThetasBalanced,2)]];

ThinnestLam   = num2cell(ThetasBalanced(1,:));  % Stacking Sequence of the thinnest lam. Patch
SSTable       = ThinnestLam;                    % SSTable initialisation


% Remove un-used Insert Indexes (=0) and if damtol or covering desig guidelines 
% are active also remove 1st ply insertion
if Constraints.Vector(3) || Constraints.Vector(8)
    InsertIndexes(InsertIndexes<=1) = [];
else
    InsertIndexes(InsertIndexes==0) = [];
end

% keyboard
if ~isempty(InsertIndexes)
  
    TenPercentAngles    = repmat([0 45 90 -45],1,ceil(length([TenPercentLoc TenPercentLoc_Add])/4)); % Ten % angle orders
    TenPercentAngles    = TenPercentAngles(length(TenPercentLoc)+1:end);
    
    % Number of angles missing to reach the thickest laminate patch
    Delta_NthetaVar     = NStruct.NthetaVar     - NStructMin.NthetaVar;
    Delta_NbalVar       = NStruct.NbalVar       - NStructMin.NbalVar;
    Delta_N10percentVar = NStruct.N10percentVar - NStructMin.N10percentVar;
    SumDelta            = Delta_NthetaVar+Delta_NbalVar+Delta_N10percentVar;
    
    i  = 1; % Theta's index
    ii = 1; % Ten % rule index
    NIndex = length(InsertIndexes);
    
    while (NIndex)>=i || (Delta_NthetaVar==0 &&  SumDelta>0) %SumDelta>0 && NIndex>=i %( (NIndex>=i && ~GuideLam) || GuideLam)
        NewLine = SSTable(1,:);     % New line of the stacking sequence table
      
        % ---
            if Delta_NthetaVar>0  % add theta's design variables and balanced (if any)

                j = InsertIndexes(i); % index at each the angle is inserted
                NewTheta = -90 + ThetasCoded_Add(i)*Constraints.DeltaAngle;

                [NewLine,SSTable] = InsertAngle(NewTheta,InsertIndexes(i),NewLine,SSTable,LamType,Constraints,MidAngle);
                Delta_NthetaVar = Delta_NthetaVar-1;


                if Constraints.Vector(2) % balanced
                    % Add theta's corresponding balanced angles                   
                    [NewLine,SSTable] = InsertAngle(-NewTheta,BalancedLoc_Add(i),NewLine,SSTable,LamType,Constraints,MidAngle);
                    Delta_NbalVar     = Delta_NbalVar - 1;
                end
                
                SSTable = [NewLine; SSTable];  % Concatenate the New line with the table
            end 
        % ---
        
        
        % ---
%         keyboard
            if Delta_N10percentVar>0 && length(TenPercentAngles)>=ii
                % Add New lines with 10% fibre angle (if any)
                
%                 Added10Percent = true;
                Added10Percent = false;
                if strcmp(LamType,'Balanced_Sym'),                Added10Percent = true; end
                if strcmp(LamType,'Sym')      && rem(i,2) == 0,   Added10Percent = true; end
                if strcmp(LamType,'Balanced') && rem(i,2) == 0,   Added10Percent = true; end
                if strcmp(LamType,'Generic')  && rem(i,4) == 0,   Added10Percent = true; end

                if Added10Percent
                    NewLine = SSTable(1,:);
                    
                    if ~Constraints.Vector(2) % not balanced
                        [NewLine,SSTable] = InsertAngle(TenPercentAngles(ii),TenPercentLoc_Add(ii),NewLine,SSTable,LamType,Constraints,MidAngle);
                        Delta_N10percentVar = Delta_N10percentVar -1;
                        
                    else
                        if TenPercentAngles(ii)==0 || TenPercentAngles(ii)==90
                            [NewLine,SSTable] = InsertAngle(TenPercentAngles(ii),TenPercentLoc_Add(ii),NewLine,SSTable,LamType,Constraints,MidAngle);
                            Delta_N10percentVar = Delta_N10percentVar -1;
                        end
                        
                        if TenPercentAngles(ii)==-45
                            [NewLine,SSTable] = InsertAngle(TenPercentAngles(ii),TenPercentLoc_Add(ii),NewLine,SSTable,LamType,Constraints,MidAngle);
                            
%                             Delta_N10percentVar = Delta_N10percentVar -1;
                            [NewLine,SSTable] = InsertAngle(TenPercentAngles(ii-2),TenPercentLoc_Add(ii-2),NewLine,SSTable,LamType,Constraints,MidAngle);
                            Delta_N10percentVar = Delta_N10percentVar -2;
                        end
                        
%                         if TenPercentAngles(ii)==45
%                             [NewLine,SSTable] = InsertAngle(TenPercentAngles(ii),TenPercentLoc_Add(ii),NewLine,SSTable,LamType,Constraints,MidAngle);
%                             
%                             Delta_N10percentVar = Delta_N10percentVar -1;
% %                             [NewLine,SSTable] = InsertAngle(TenPercentAngles(ii-2),TenPercentLoc_Add(ii-2),NewLine,SSTable,LamType,Constraints,MidAngle);
% %                             Delta_N10percentVar = Delta_N10percentVar -2;
%                         end
                        
                    end
                    SSTable = [NewLine; SSTable];
                    ii = ii+1;
                end       
                
            end
        % ---
        
        SumDelta = Delta_NthetaVar+Delta_NbalVar+Delta_N10percentVar;
        i  = i+1;
    end
end


%% Apply symmetry
if Constraints.Vector(1) % symmetric
    if isempty(Thetas_MidCoded)
        SSTable = [SSTable fliplr(SSTable)]; % remove midplane from symmetry
    else
        SSTable = [SSTable fliplr(SSTable(:,1:end-NStruct.NMidPlane  ))]; % remove midplane from symmetry
    end
end


%% Remove redundant lines in SSTable
% if size(SSTable,1)>1
    NplySS = nan*ones(size(SSTable,1),1);
    for j = 1:size(SSTable,1)
        NplySS(j) = length(find(~cellfun(@isempty,SSTable(j,:))));
    end
    [NplySS,Index] = unique(NplySS,'stable');
    SSTable = SSTable(Index,:);
% end


%% Remove empty Columns (ply not used at all)

Nrows = length(NplySS);
Temp  = sum(cellfun(@isempty,SSTable));

EmptyCol_Index = find(Nrows==Temp);
if ~isempty(EmptyCol_Index)
    SSTable(:,EmptyCol_Index) = [];
end

end



% Local Insertion angle function
function [NewLine,SSTable] = InsertAngle(NewAngle,InsertIndex,NewLine,SSTable,LamType,Constraints,MidAngle)

    if nargin == 6
        MidAngle = false;
    end

    Nrows = size(SSTable,1);
    II    = InsertIndex; % Short name for sake of clarity

    if length(SSTable)>II
        NewLine = [NewLine(1:II-1)       NewAngle           NewLine(II:end)];
        SSTable = [SSTable(:,1:II-1)     cell(Nrows,1)      SSTable(:,II:end)];
    else
        
        if MidAngle || strcmp(LamType,'Generic') && (Constraints.Vector(3) || Constraints.Vector(8)) 
            % special case where last angle from the thinnest laminate must remain the last angle of the SSTable
            NewLine = [NewLine(1:end-1)     NewAngle        NewLine(end)];
            SSTable = [SSTable(:,1:end-1)   cell(Nrows,1)   SSTable(:,end)];
            
        else
            NewLine = [NewLine NewAngle ];
            SSTable = [SSTable cell(Nrows,1)];
        end
        
    end
end