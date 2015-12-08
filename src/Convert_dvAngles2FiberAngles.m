% =====                                                              ==== %
%               Converts the ply angle design variable in the full
%                           length staking sequence
% =====                                                              ==== %

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


function [FiberAngles] = Convert_dvAngles2FiberAngles(GuideAngles,DropsLoc,ShuffleLoc,LamType)

%% Input Check
if ~isnumeric(GuideAngles) || ~isreal(GuideAngles) || ~isvector(GuideAngles)
    error('Input angles (GuideAngles) must be numeric and real')
end
if size(GuideAngles,2) == 1
    GuideAngles = GuideAngles';
end


StringList = {'Balanced_Sym' 'Sym' 'Balanced' 'Generic'};
Index      = find(strncmp(LamType,StringList,12),1);

if isempty(Index), 
    error('Non-recognised LamType used as input');
end

%% Concatenation
ply_angle = GuideAngles;
if strcmp(LamType,'Generic'),  
    ply_angle(DropsLoc) = [];
    FiberAngles = ply_angle';    
end


if strcmp(LamType,'Sym'), 
    ply_angle(DropsLoc) = []; 
    FiberAngles = [ply_angle, fliplr(ply_angle)]'; 
end


% For Balanced Lam. need to reconstruct SS by inserting ply angle (+-)pairs
if strcmp(LamType,'Balanced_Sym') || strcmp(LamType,'Balanced')
    FiberAngles    = [ply_angle; [1:length(ply_angle)]; [1:length(ply_angle)]];
    BalancedAngles = [-ply_angle' ShuffleLoc' [1:length(ShuffleLoc)]'];
    BalancedAngles = sortrows(BalancedAngles,2)';
    
    for j = 1:length(BalancedAngles)
        if BalancedAngles(2,j)>size(FiberAngles,2)
            FiberAngles = [FiberAngles BalancedAngles(:,j)];
        else
            FiberAngles = [FiberAngles(:,1:BalancedAngles(2,j)-1) BalancedAngles(:,j)  FiberAngles(:,BalancedAngles(2,j):end)];
        end
    end
    
    for iDrop = 1:length(DropsLoc)
        NewDropsLoc = find(FiberAngles(3,:)==DropsLoc(iDrop));
        FiberAngles(:,NewDropsLoc) = [];                                    %#ok<FNDSB> % left for clarity
    end
    
    FiberAngles = FiberAngles(1,:);
    
    if strcmp(LamType,'Balanced_Sym')
        FiberAngles = [FiberAngles, fliplr([FiberAngles])]';
    end
end


end