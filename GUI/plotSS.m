% ===                                                                   === 
%                     Stacking Sequence Plotting function                       
%
% plotSS(Output,Dim,PatchXYZ,varargin)
%
% - Output: structure returned by RetrieveSS
% - Dim:    Either set to 2 or 3, corresponding to 2D or 3D plot
% - PatchXYZ: (optional) 3D location of laminate patches
% ===                                                                   === 

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


function [] = plotSS(Output,PatchXYZ,varargin)


[~,GuideIndex] = max(cellfun(@length,Output.SS));
GuideIndex     = GuideIndex(1);

Fig = figure;
set(Fig,'Units','normalized','Position',[0.1,0.1,0.8,0.8])

Axis1 = axes('Position',[0.05,0.075,0.5,0.85],'XTick',[0:1:length(Output.SS)]);
xlabel('X')
ylabel('Y')
zlabel('Ply #')
hold all

Map = colormap(pmkmp(19,'CubicYF'));
caxis([-90 90])
colormap(Map); 


Axis2 = axes('Position',[0.6,0.075,0.3,0.85],'XTick',[0:1:length(Output.SS)],'xlim', [0 length(Output.SS)+1]);
xlabel('Patch #')
zlabel('Ply #')
hold all

Map = colormap(pmkmp(19,'CubicYF'));
caxis([-90 90])
colormap(Map); 
h = colorbar;
set(h, 'ticks', [-90:10:90])
set(h, 'ylim', [-90 90])


% Reconstruct Stacking Sequences
NUniqueLam  = length(Output.SS);
NPliesLam   = cellfun(@length,Output.SS);
[NGuide,GuideIndex] = max(NPliesLam);
GuideAngles = num2cell(Output.SS{GuideIndex});


for i = 1:NUniqueLam
    if i == GuideIndex || NPliesLam(i) == NGuide 
        SS{i} = GuideAngles;    
    else 
        index = 1;
        for j=1:NGuide
            if index<=NPliesLam(i) && GuideAngles{j} == Output.SS{i}(index)
                SS{i}{j} = Output.SS{i}(index);
                index = index +1;
            else
                SS{i}{j} = [];
            end
        end
    end
end

if strcmp(Output.LamType,'Sym') || strcmp(Output.LamType,'Balanced_Sym')
    HalfNum = length(SS{1})/2;
    for i = 1:NUniqueLam
        SS{i} = SS{i}(1:HalfNum);
    end
end


if nargin~=2 % no location data have been given
    for i=1:NUniqueLam
        PatchXYZ{i}.X = i+[-1 0 0 -1];
        PatchXYZ{i}.Y = [0 0 1 1];
        PatchXYZ{i}.Z = [0 0 0 0];
    end
end


Opacity = ones(length(SS{1}),1);



axes(Axis1)
for j=1:length(SS)
    X = PatchXYZ{j}.X;
    Y = PatchXYZ{j}.Y;
    Z = PatchXYZ{j}.Z;
     
    for i = 1:length(SS{j})
        if ~isempty(SS{j}{i})
            F{j,i} = fill3(X,Y,Z+i,SS{j}{i},'facealpha',Opacity(1) );
            if j == GuideIndex
                T{i} = text(mean(X),mean(Y),mean(Z+i)+0.1,num2str(SS{j}{i}));
            end
        end
    end    
end
view(3)


for i=1:NUniqueLam
    PatchXYZ{i}.X = i+[-1 0 0 -1];
    PatchXYZ{i}.Y = [0 0 1 1];
    PatchXYZ{i}.Z = [0 0 0 0];
end
    
axes(Axis2)
for j=1:length(SS)
    X = PatchXYZ{j}.X;
    Y = PatchXYZ{j}.Y;
    Z = PatchXYZ{j}.Z;

    X(3) = X(2);     X(4) = X(1);
    Y(3) = Y(2);     Y(4) = Y(1);
    Z(3) = Z(2)+0.5; Z(4) = Z(1)+0.5;

    for i = 1:length(SS{j})
        if ~isempty(SS{j}{i})
            F2{j,i} = fill3(X,Y,Z+i,SS{j}{i},'facealpha',Opacity(1) );
%             if j == GuideIndex
%                 T2{j,i} = text(mean(X),mean(Y),mean(Z+i),num2str(SS{j}{i}));
%             end
        end
    end  
end

for  i = 1:length(SS{GuideIndex})
    X = PatchXYZ{end}.X+1;
    Y = PatchXYZ{end}.Y;
    Z = PatchXYZ{end}.Z+0.5;
    T2{i} = text(mean(X),mean(Y),mean(Z+i),num2str(SS{GuideIndex}{i}));
end
view(0,0)


    
% keyboard
% mytxt = uicontrol('Parent',Fig,'Style','text','String',['Total # of plies:' num2str(length(SS{j}))],'units','normalized','Position',[0.75,0.9,0.25,0.05]);

b = uicontrol('Parent',Fig,'Style','slider','units','normalized','Position',[0.95,0.125,0.025,0.75],...
              'value',length(SS{j}), 'min',0, 'max',length(SS{j}),'Callback', @(es,ed) SetOpacity(es,ed,F,F2));
          
ShowAngles = uicontrol('Parent',Fig,'Style','checkbox','String','Display Angles','Value',1,'units','normalized','Position',[0.85,0.95,0.25,0.05],'Callback', @(es,ed) SetAngle(es,ed,T));
          