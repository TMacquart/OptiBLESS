function SSTable = NEWComputeSSTable(UniqueNPly,Thetas,DropsIndexes,BalancedLoc,TenPercentLoc,Thetas_Mid,LamType)   

keyboard

% if strcmp(LamType,'Generic'), Coeff=1; end
% if strcmp(LamType,'Balanced') || strcmp(LamType,'Sym'), Coeff=2; end
% if strcmp(LamType,'Balanced_Sym'), Coeff=4; end
% PatchNPly = sort(UniqueNPly*Coeff,'descend');

PatchNPly     = sort(UniqueNPly,'descend');
SymbolicThetas =  [1:length(Thetas)];

index = 1;
SymbolicSSTable = cell(length(PatchNPly),max(PatchNPly));
for j=0:length(DropsIndexes)
    if ~isempty(find(((length(Thetas)-j)*Coeff-PatchNPly)==0,1))
        
        SymbolicFiberAngles = Convert_dvAngles2FiberAngles(SymbolicThetas,DropsIndexes(1:j),BalancedLoc,LamType)';
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

SSTable = cell(size(SymbolicSSTable));
for irow = 1:size(SymbolicSSTable,1)
    for icol = 1:size(SymbolicSSTable,2)
        SSTable{irow,icol} = Thetas(abs(SymbolicSSTable{irow,icol}))*sign(SymbolicSSTable{irow,icol});
    end
end


end