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
%               Converts the ply angle design variable in the full
%                           length staking sequence
% =====                                                              ==== %

function [FiberAngles] = Convert_dvAngles2FiberAngles(ply_angle,ShuffleLoc,LamType)


%% Input Check
if ~isnumeric(ply_angle) || ~isreal(ply_angle) || ~isvector(ply_angle)
    error('Input angles (ply_angle) must be numeric and real')
end
if size(ply_angle,2) == 1
    ply_angle = ply_angle';
end

StringList = {'Balanced_Sym' 'Sym' 'Balanced' 'Generic'};
Index      = find(strncmp(LamType,StringList,12),1);

if isempty(Index), 
    error('Non-recognised LamType used as input');
end

%% Concatenation
if strcmp(LamType,'Generic'),        FiberAngles = ply_angle';                          end
if strcmp(LamType,'Sym'),            FiberAngles = [ply_angle, fliplr(ply_angle)]';     end

try
if strcmp(LamType,'Balanced_Sym') || strcmp(LamType,'Balanced')
    FiberAngles    = ply_angle;
    BalancedAngles = [-ply_angle' ShuffleLoc'];
    BalancedAngles = sortrows(BalancedAngles,2);
    
    for j = 1:length(BalancedAngles)
        if BalancedAngles(j,2)>length(FiberAngles)
            FiberAngles = [FiberAngles BalancedAngles(j,1)];
        else
            FiberAngles = [FiberAngles(1:BalancedAngles(j,2)-1) BalancedAngles(j,1)  FiberAngles(BalancedAngles(j,2):end)];
        end
    end
    if strcmp(LamType,'Balanced_Sym')
        FiberAngles = [FiberAngles, fliplr([FiberAngles])]';
    end
end
catch
    keyboard
end

end