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


% =====                                                              ==== %
%               Calculates N Feasible ply drops from a laminate           %
%
%  Returns a possible droped laminate and the drop plies indexes          %
%  Terence Macquart (26/04/2015)                                          %
% =====                                                              ==== %

function [Indexes,Laminate] = Generating_DropIndexes (ply_angle,Ndrop,ConstraintVector)

% -- ply angles defined from bottom (left) to top (right)
NoSolution = false;
Laminate = num2cell(ply_angle');

for iDrop = 1 : Ndrop
    
    NonEmptyCell = ~cellfun(@isempty,Laminate(:,iDrop));
    FeasibleDrop = find(NonEmptyCell);
    IndexArray   = [FeasibleDrop(1:end-2) FeasibleDrop(3:end) FeasibleDrop(2:end-1) zeros(length(FeasibleDrop)-2,1)];  % 1col - 2col = Delta, 3col = index
    
    % ---
    if ConstraintVector(3) == 1          % apply disorientation constrainst to PlyDrop
        DeltaDrops   = abs( cell2mat(Laminate(IndexArray(:,1),iDrop)) - cell2mat(Laminate(IndexArray(:,2),iDrop)) );
        IndexArray   = IndexArray(DeltaDrops<=45,:);
    end
    % ---
    
    % ---
    if ConstraintVector(4) == 1          % apply contiguity constraint
        DeltaDrops   = abs( cell2mat(Laminate(IndexArray(:,1),iDrop)) - cell2mat(Laminate(IndexArray(:,2),iDrop)) );
        IndexArray   = IndexArray(DeltaDrops>=5,:);
    end
    % ---
    
    % ---
    if ConstraintVector(6) && iDrop>1    % apply internal continuity constraints
        for i = 1 : length(DropIndex)
            for j = 1 : length(IndexArray(:,1))
                if DropIndex(i)>IndexArray(j,1) && DropIndex(i)<IndexArray(j,2)
                    IndexArray(j,4) = IndexArray(j,4) +1;
                end
            end
        end
        IndexArray = IndexArray(IndexArray(:,4)<3,:);
    end
    % ---
    
    % ---
    if ConstraintVector(2)               % apply 10% rule constraint
        lam_vect   = cell2mat(Laminate(:,iDrop));
        N10percent = round(0.1*length(lam_vect));
        N0Plies    = length(find(lam_vect==0));
        N90Plies   = length(find(abs(lam_vect)==90));
        N45Plies   = length(find(abs(lam_vect)==45));
        
        if N0Plies <= N10percent % cannot remove a 0 ply
            temp = [abs(cell2mat(Laminate(IndexArray(:,3)))) IndexArray(:,3)];
            IndexArray(temp(:,1)==0,:) = [];
        end
        
        if N45Plies <= N10percent
            temp = [abs(cell2mat(Laminate(IndexArray(:,3)))) IndexArray(:,3)];
            IndexArray(temp(:,1)==45,:) = [];
        end
        
        if N90Plies <= N10percent
            temp = [abs(cell2mat(Laminate(IndexArray(:,3)))) IndexArray(:,3)];
            IndexArray(temp(:,1)==90,:) = [];
        end
    end
    % ---
    
    if isempty(IndexArray),
        display('Imposible to remove any more ply (probably due to constraints)');
        NoSolution = true;
        break
    end
    
    FeasibleDrop        = IndexArray(:,3);
    DropIndex(iDrop)    = FeasibleDrop( ceil(length(FeasibleDrop)*rand) );
    Laminate(:,iDrop+1) = Laminate(:,iDrop);
    Laminate(DropIndex(iDrop),iDrop+1) = {[]};
    
    if length(unique(DropIndex)) ~=length(DropIndex), keyboard; end
end


if NoSolution
    Laminate = [];
    Indexes  = [];
else
    Laminate = Laminate(:,end);
    Indexes  = DropIndex;
end

end










