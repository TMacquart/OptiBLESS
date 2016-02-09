function FEASIBLE = Check_InternalContinuity(SSTable)
FEASIBLE = true;

SSDrops  = sort(find(cellfun(@isempty,SSTable(end,:)))); % check last line of SS Table
if SSDrops>3
    for i = 1 : length(SSDrops)-3
        deltaLoc = diff(SSDrops(i:i+3));
        if sum(abs(deltaLoc))==3
            FEASIBLE = false;
            break
        end
    end
end

end