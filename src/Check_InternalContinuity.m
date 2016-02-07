function FEASIBLE = Check_InternalContinuity(SSTable)

SSDrops = sort(find(cellfun(@isempty,SSTable(end,:)))); % check last line of SS Table
for i = 1 : length(SSDrops)-3
    deltaLoc = diff(SSDrops(i:i+3));
    if sum(abs(deltaLoc))==3
        FEASIBLE = false;
        break
    end
end

end