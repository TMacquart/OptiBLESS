function NGeoConstraints = CheckContinuity(SS_Patch,PatchConnectivity)

NPatch = size(SS_Patch,1);
NplySS = nan*ones(size(SS_Patch,1),1);
for j = 1:size(SS_Patch,1)
    NplySS(j) = length(find(~cellfun(@isempty,SS_Patch(j,:))));
end
[NGuide,GuideIndex] = max(NplySS);


% check the geometrical continuity of each ply
NGeoConstraints = 0;
for iply = 1:NGuide
    
    Ply        = SS_Patch(:,iply);
    PlyPatches = find(~cellfun(@isempty, Ply)); % Patches which the ply is covering
    
    if length(PlyPatches) ~= NPatch                                         % possible discontinuity
        
        CoveredPatch  = PlyPatches;                                         % we are progressively removing reachable patches 
        ContinuousPly = GuideIndex;                                         % patch numbers that are connected and in the ply
        CoveredPatch(CoveredPatch == ContinuousPly) = [];                          
        ConnectedPatch = find(PatchConnectivity(ContinuousPly,:));          % patches next to the guide              
        
       
        while ~isempty(CoveredPatch)                                        % as long as some patches are not accounted for
            OldConnectedPatch = ConnectedPatch;
            
            % identify connected patch that are aslo covered by the ply
            for j=1:length(ConnectedPatch)
                PatchIndexes = find(CoveredPatch==ConnectedPatch(j)); % find Covered patch reachable from the guide extending ply        
                ContinuousPly = [ContinuousPly CoveredPatch(PatchIndexes)]; %#ok<AGROW>
                CoveredPatch(PatchIndexes)=[];
            end
            
            % compute new connected patches
            NewConnection = []; 
            for j=1:length(ContinuousPly)
                NewConnection = [NewConnection find(PatchConnectivity(ContinuousPly(j),:))];
            end
            ConnectedPatch = unique(NewConnection);
            
            % remove patches already accounted for
            for j=1:length(ContinuousPly)
                ConnectedPatch(ConnectedPatch == ContinuousPly(j)) = []; 
            end
            
            % break from while loop if no change 
            if length(OldConnectedPatch) == length(ConnectedPatch)
                if sum(sort(OldConnectedPatch) - sort(ConnectedPatch)) == 0
                    NGeoConstraints = NGeoConstraints +1;
                    break
                end
            end
            
        end
    end
end

end

