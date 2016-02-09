function SSTable = NewComputeSSTable(ThetasCoded,DropsIndexes,BalancedLoc,TenPercentLoc,Thetas_MidCoded,LamType,Constraints)   


% convert Theta's into angles
Thetas = ThetasCoded;
if ~Constraints.Vector(1)  % if not Damtol
    Thetas = -90 + Thetas*Constraints.DeltaAngle;
else
    Thetas(2:end) = -90 + Thetas(2:end)*Constraints.DeltaAngle;
    R = [-1 1];
    Thetas(1) = 45*R(Thetas(1));
%     if strcmp(LamType,'Generic') % convert to closest
%         if abs(Thetas(end)-45) < abs(Thetas(end) +45) % closer to +45
%             Thetas(end) = 45;
%         else
%             Thetas(end) = -45;
%         end
%     end
end



% For Balanced Lam. need to reconstruct SS by inserting ply angle (+-)pairs
ThetasBalanced = [Thetas'; [1:length(Thetas)]; [1:length(Thetas)]];

if strcmp(LamType,'Balanced_Sym') || strcmp(LamType,'Balanced')
    BalancedAngles = [-Thetas BalancedLoc [1:length(BalancedLoc)]'];
    BalancedAngles = sortrows(BalancedAngles,2)';
    
    for j = 1:length(BalancedAngles(1,:))
        if BalancedAngles(2,j)>size(ThetasBalanced,2)
            ThetasBalanced = [ThetasBalanced BalancedAngles(:,j)];
        else
            ThetasBalanced = [ThetasBalanced(:,1:BalancedAngles(2,j)-1) BalancedAngles(:,j)  ThetasBalanced(:,BalancedAngles(2,j):end)];
        end
    end
end



% For 10% angles need to reconstruct SS by inserting ply angle
if Constraints.Vector(2)
    NAngles = length(TenPercentLoc);
    TenPercentAngles = repmat([-45 45 0 90],1,ceil(NAngles/4));
    TenPercentAngles = TenPercentAngles(1:NAngles);
    TenPercentAngles = [TenPercentAngles' TenPercentLoc];
    index = max(ThetasBalanced(3,:)) + 1;
    for j = 1:NAngles
        TenPercentAngles(j,3) = index;
        
        if Constraints.Balanced
            if j~=1 && rem(j,5)~=0, 
                index = index +1; 
            end
        else
            index = index +1; 
        end
    end
    TenPercentAngles = sortrows(TenPercentAngles,2)';
    
    for j = 1:NAngles
        if TenPercentAngles(2,j)>size(ThetasBalanced,2)
            ThetasBalanced = [ThetasBalanced TenPercentAngles(:,j)];
        else
            ThetasBalanced = [ThetasBalanced(:,1:TenPercentAngles(2,j)-1) TenPercentAngles(:,j)  ThetasBalanced(:,TenPercentAngles(2,j):end)];
        end
    end
end



% Add Mid plane angles
if ~isempty(Thetas_MidCoded)
    
    LastIndex(1,1) = max(ThetasBalanced(2,:));
    LastIndex(2,1) = max(ThetasBalanced(3,:));
    
    if strcmp(LamType,'Generic')
        Thetas_Mid = (-90 + Thetas_MidCoded*Constraints.DeltaAngle)';
        Thetas_Mid = [Thetas_Mid; zeros(1,length(Thetas_Mid)); [1:length(Thetas_Mid(1,:))] + LastIndex(2,1)];
    end
    
    if strcmp(LamType,'Balanced')
        R = [0 90];
        Thetas_Mid = R(Thetas_MidCoded);
        Thetas_Mid = [Thetas_Mid; zeros(1,length(Thetas_Mid)); [1:length(Thetas_Mid(1,:))] + LastIndex(2,1)];
    end
    
    if strcmp(LamType,'Sym')
        Thetas_Mid = -90 + Thetas_MidCoded(1,1)*ones(1,length(Thetas_MidCoded))*Constraints.DeltaAngle;
        Thetas_Mid = [Thetas_Mid; zeros(1,length(Thetas_Mid)); ones(1,length(Thetas_Mid(1,:)))*(1+LastIndex(2,1))];
    end
    
    if strcmp(LamType,'Balanced_Sym')
        R = [0 90];
        Thetas_Mid = R(Thetas_MidCoded(1))*ones(1,length(Thetas_MidCoded));
        Thetas_Mid = [Thetas_Mid; zeros(1,length(Thetas_Mid)); ones(1,length(Thetas_Mid(1,:)))*(1+LastIndex(2,1))];
    end
   
    ThetasBalanced = [ThetasBalanced Thetas_Mid];
else
    
end
    

% Constructing the stacking sequence
ThetasBalanced = [ThetasBalanced; [1:length(ThetasBalanced)]];

Guide   = num2cell(ThetasBalanced(1,:));
SSTable = Guide;

DropsIndexes = unique(DropsIndexes,'stable');

for iDrop = 1:length(DropsIndexes)
    if ~isempty(ThetasBalanced(3,:)==DropsIndexes(iDrop))
        SSTable = [SSTable ; Guide];
        for j = 1:iDrop
            DropLoc = find(DropsIndexes(j) == ThetasBalanced(3,:));
            SSTable(end,DropLoc) = cell(1,length(DropLoc));
            
        end
    end
end

if Constraints.Sym
    if isempty(Thetas_MidCoded)
         SSTable = [SSTable fliplr(SSTable)]; % remove midplane from symmetry
    else
        SSTable = [SSTable fliplr(SSTable(:,1:end-length(Thetas_MidCoded(1,:)) ))]; % remove midplane from symmetry
    end
end


end