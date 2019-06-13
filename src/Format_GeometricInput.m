function [PatchConnectivity] = Format_GeometricInput(Patch)


%% Number Edges
NPatch      = length(Patch);         % Total number of patches
PatchEdgeId = zeros(NPatch,4);          % save the 4 Edge Identifier for all patches
Edges = zeros(4*NPatch,6);              % Pre-allocation

PatchEdgeId(1,:) = [1 2 3 4];
[Edges(1:4,:)]   = ReturnEdges(Patch{1}); % 1st patch = 4 new edges
EdgeId = 5;                                     % Edge Id Increment


for iPatch = 2:NPatch  % check if edges does not already exist before accounting it
    
    [PatchEdges]     = ReturnEdges(Patch{iPatch});
    RevertPatchEdges = [PatchEdges(:,4:6) PatchEdges(:,1:3)];
    
    for i = 1:4
        
        [Member1,Index1] = ismember(PatchEdges(i,:),Edges,'rows');
        [Member2,Index2] = ismember(RevertPatchEdges(i,:),Edges,'rows');
        
        if ~Member1 && ~Member2 % new edge
            Edges(EdgeId,:) = PatchEdges(i,:);
            PatchEdgeId(iPatch,i) =  EdgeId;
            EdgeId = EdgeId+1;
            
        else
            if Member1
                PatchEdgeId(iPatch,i) = Index1;
            else
                PatchEdgeId(iPatch,i) = Index2;
            end
        end
    end
    
    
end
EdgeId = EdgeId-1;
Edges  = Edges(1:EdgeId,:);


if 0 % Edge Visual Check
    figure
    hold all
    for iEdge= 1:EdgeId
        % check geometry by plotting edges
        plot3([Edges(iEdge,[1 4])],[Edges(iEdge,[2 5])],[Edges(iEdge,[3 6])])
        text(mean([Edges(iEdge,[1 4])])+0.0,mean([Edges(iEdge,[2 5])])+0.0,mean([Edges(iEdge,[3 6])]),num2str(iEdge))
    end
end



%% Normalised Edge Vector
EdgeDir = Edges(:,1:3)-Edges(:,4:6);
RowNorm = arrayfun(@(idx) norm(EdgeDir(idx,:)), 1:size(EdgeDir,1))';
for i=1:size(EdgeDir,1)
    EdgeDir(i,:) = EdgeDir(i,:)/RowNorm(i,:);
end

UniqueEdgeDir      = unique(EdgeDir,'rows');
UniqueEdgeDir(:,4) = 1:size(UniqueEdgeDir,1);

for j = 1:size(UniqueEdgeDir,1)
    [Bool,Index] = ismember(-UniqueEdgeDir(j,1:3),UniqueEdgeDir(:,1:3),'rows');
    if  Bool
        UniqueEdgeDir(Index,4) = UniqueEdgeDir(j,4);
    end
end


[~,EdgeDir(:,4)] = ismember(EdgeDir,UniqueEdgeDir(:,1:3),'rows');
EdgeDir(:,4)     = UniqueEdgeDir(EdgeDir(:,4),4);

if 1 % Edge dir visual check
    colorslist = [{'blue'},{'red'},{'green'},{'black'}];
    figure
    hold all
    for iEdge= 1:EdgeId
        plot3([Edges(iEdge,[1 4])],[Edges(iEdge,[2 5])],[Edges(iEdge,[3 6])],colorslist{EdgeDir(iEdge,4)})
        text(mean([Edges(iEdge,[1 4])])+0.0,mean([Edges(iEdge,[2 5])])+0.0,mean([Edges(iEdge,[3 6])]),num2str(iEdge))
    end
end




%% Check the non-trivial Edge connectivity (subpart of other vertices)
EdgeConnectivity= [];
for j= 1:max(UniqueEdgeDir(:,4))
    EdgeIndex = find(EdgeDir(:,4)==j); % all edge with same vector direction
    
    for i = 1:length(EdgeIndex)
        p1  = Edges(EdgeIndex(i),1:3);
        p2  = Edges(EdgeIndex(i),4:6);
        Dir = find (p2-p1,1);
        
        for ii = i+1:length(EdgeIndex)
            p3 = Edges(EdgeIndex(ii),1:3);
            p4 = Edges(EdgeIndex(ii),4:6);
            
            p3Hat = interp1([p1(Dir) p2(Dir)],[p1; p2],p3(Dir));
            p4Hat = interp1([p1(Dir) p2(Dir)],[p1; p2],p4(Dir));

            if sum(abs(p3-p3Hat))<1e-10 && isnan(p4Hat(1)) % p4 outside the edge
                
                if abs(norm(p4)-norm(p1)) < abs(norm(p4)-norm(p2))
                    overlap = norm(p3-p1);
                else
                    overlap = norm(p3-p2);
                end

                if overlap >= 0.25*norm(p2-p1) &&  overlap >= 0.25*norm(p4-p3)
                    EdgeConnectivity =[EdgeConnectivity; EdgeIndex(i) EdgeIndex(ii)];
                end
            end
            if sum(abs(p4-p4Hat))<1e-10 && isnan(p3Hat(1))
                
                if abs(norm(p3)-norm(p1)) < abs(norm(p3)-norm(p2))
                    overlap = norm(p4-p1);
                else
                    overlap = norm(p4-p2);
                end
                
                if overlap >= 0.25*norm(p2-p1) &&  overlap >= 0.25*norm(p4-p3)
                    EdgeConnectivity =[EdgeConnectivity; EdgeIndex(i) EdgeIndex(ii)];
                end
            end 
            
            if sum(abs(p3-p3Hat))<1e-10 && sum(abs(p4-p4Hat))<1e-10 % Edge subpart
                EdgeConnectivity =[EdgeConnectivity; EdgeIndex(i) EdgeIndex(ii)];
            end
        end
    end
end




%% Compute connection matrix between patches (same edge == connected or EdgeConnectivity exist)
PatchConnectivity = zeros(NPatch);
for iPatch = 1:NPatch
    for iEdge = 1:4
        % find patch with same edges
        [irow,~]=find(PatchEdgeId==PatchEdgeId(iPatch,iEdge));
        for j=1:length(irow)
            if irow(j)~=iPatch
                PatchConnectivity(iPatch,irow(j))=1;
            end
        end
        
        % find connected edge
        [ii,jj]    = find(PatchEdgeId(iPatch,iEdge)==EdgeConnectivity);
        
        for j=1:length(ii)
            if jj(j) ==1
                [irow,~] = find(PatchEdgeId==EdgeConnectivity(ii(j),2));
            else
                [irow,~] = find(PatchEdgeId==EdgeConnectivity(ii(j),1));
            end
            PatchConnectivity(iPatch,irow)=1;
        end
    end
end
end

function [PatchEdges] = ReturnEdges(Patch)
    PatchEdges = zeros(4,6);
    for iEdge = 1:4
        if iEdge == 4
            PatchEdges(iEdge,:) =  [Patch.X(iEdge) Patch.Y(iEdge) Patch.Z(iEdge) Patch.X(1) Patch.Y(1) Patch.Z(1)];
        else
            PatchEdges(iEdge,:) =  [Patch.X(iEdge) Patch.Y(iEdge) Patch.Z(iEdge) Patch.X(iEdge+1) Patch.Y(iEdge+1) Patch.Z(iEdge+1)];
        end
    end
    
    PatchEdges = round(PatchEdges*1e8)*1e-8;
end