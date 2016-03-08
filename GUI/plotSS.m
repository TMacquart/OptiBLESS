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
% ===                                                                   === 
%                     Stacking Sequence Plotting function                       
%
% plotSS(Output,PatchXYZ,varargin)
%
% - Output: structure returned by RetrieveSS
% - Dim:    Either set to 2 or 3, corresponding to 2D or 3D plot
% - PatchXYZ: (optional) 3D location of laminate patches
% ===                                                                   === 

function [] = plotSS(Output,Scale,PatchXYZ,varargin)

Fig = figure;
FigPosition    = [0.0,0.0,0.8,0.8]*Scale;
FigPosition(1) = 0.98-(FigPosition(end-1)); %  |0.98|   are used to fixed the window on the top right corner (you can change these values acccording to your preferences)
FigPosition(2) = 0.95-(FigPosition(end));   %  |0.95|

set(Fig,'Units','normalized','Position',FigPosition)

Axis1 = axes('Position',[0.05,0.075,0.5,0.85],'XTick',[0:1:length(Output.SS_Patch)]);
xlabel('X')
ylabel('Y')
zlabel('Ply #')
hold all

Map = colormap(pmkmp(19,'CubicYF'));
caxis([-90 90])
colormap(Map); 


Axis2 = axes('Position',[0.6,0.075,0.3,0.85],'XTick',[0:1:size(Output.SS_Patch,1)],'xlim', [0 size(Output.SS_Patch,1)+1]);
xlabel('Patch #')
zlabel('Ply #')
hold all

Map = colormap(pmkmp(19,'CubicYF'));
caxis([-90 90])
colormap(Map); 
h = colorbar;
set(h, 'ticks', [-90:10:90])
set(h, 'ylim', [-90 90])


SS = Output.SS_Patch;
NplySS = nan*ones(size(SS,1),1);
for j = 1:size(SS,1)
    NplySS(j) = length(find(~cellfun(@isempty,SS(j,:))));
end
[~,GuideIndex] = max(NplySS);
GuideIndex = GuideIndex(1);
Opacity = ones(size(SS,1),1);


if nargin~=3 % no location data have been given
    for i=1:size(SS,1)
        PatchXYZ{i}.X = i+[-1 0 0 -1];
        PatchXYZ{i}.Y = [0 0 1 1];
        PatchXYZ{i}.Z = [0 0 0 0];
    end
end

axes(Axis1)
for j=1:size(SS,1)
    X = PatchXYZ{j}.X;
    Y = PatchXYZ{j}.Y;
    Z = PatchXYZ{j}.Z;
     
    for i = 1:length(SS(j,:))
        if ~isempty(SS{j,i})
            F{j,i} = fill3(X,Y,Z+i,SS{j,i},'facealpha',Opacity(1) );
            if j == GuideIndex
                T{i} = text(mean(X),mean(Y),mean(Z+i)+0.1,num2str(SS{j,i}));
            end
        end
    end    
end
view(3)

for i=1:size(SS,1)
    PatchXYZ{i}.X = i+[-1 0 0 -1];
    PatchXYZ{i}.Y = [0 0 1 1];
    PatchXYZ{i}.Z = [0 0 0 0];
end
    
axes(Axis2)
for j=1:size(SS,1)
    X = PatchXYZ{j}.X;
    Y = PatchXYZ{j}.Y;
    Z = PatchXYZ{j}.Z;

    X(3) = X(2);     X(4) = X(1);
    Y(3) = Y(2);     Y(4) = Y(1);
    Z(3) = Z(2)+0.5; Z(4) = Z(1)+0.5;

    for i = 1:size(SS,2)
        if ~isempty(SS{j,i})
            F2{j,i} = fill3(X,Y,Z+i,SS{j,i},'facealpha',Opacity(1) );
        end
    end  
end

for  i = 1:length(SS(GuideIndex,:))
    X = PatchXYZ{end}.X+1;
    Y = PatchXYZ{end}.Y;
    Z = PatchXYZ{end}.Z+0.5;
    T2{i} = text(mean(X),mean(Y),mean(Z+i),num2str(SS{GuideIndex,i}));
end
view(0,0)


    
% keyboard
% mytxt = uicontrol('Parent',Fig,'Style','text','String',['Total # of plies:' num2str(length(SS{j}))],'units','normalized','Position',[0.75,0.9,0.25,0.05]);

b = uicontrol('Parent',Fig,'Style','slider','units','normalized','Position',[0.95,0.125,0.025,0.75],...
              'value',length(SS(j,:)), 'min',0, 'max',length(SS(j,:)),'Callback', @(es,ed) SetOpacity(es,ed,F,F2));
          
ShowAngles = uicontrol('Parent',Fig,'Style','checkbox','String','Display Angles','Value',1,'units','normalized','Position',[0.85,0.95,0.25,0.05],'Callback', @(es,ed) SetAngle(es,ed,T));
          