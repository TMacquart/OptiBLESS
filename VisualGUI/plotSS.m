function [] = plotSS(Output,Dim,PatchXYZ,varargin)

Fig = figure;

axes('Position',[0.125,0.125,0.75,0.75],'XTick',[0:1:length(Output.SS)])
xlabel('Patch #')
zlabel('Ply #')
hold all

Map = colormap(pmkmp(19,'CubicYF'));
Map(:,3) = -90:180/(size(Map,1)-1):90;

% Reconstruct Stacking Sequences
NUniqueLam  = length(Output.SS);

for i = 1:NUniqueLam
    if i == 1
        SS{1} = num2cell(Output.SS{1});    
    else 
        index = 1;
        for j=1:length(SS{1})
            if SS{1}{j} == Output.SS{i}(index)
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


if nargin~=3 % no location data have been given
    for i=1:NUniqueLam
        PatchXYZ{i}.X = i+[-1 0 0 -1];
        PatchXYZ{i}.Y = [0 0 1 1];
        PatchXYZ{i}.Z = [0 0 0 0];
    end
end


Opacity = ones(length(SS{1}),1);

if Dim==2 % 2D
    for j=1:length(SS)
        X = PatchXYZ{j}.X;
        Y = PatchXYZ{j}.Y;
        Z = PatchXYZ{j}.Z;
        
        
        X(3) = X(2);     X(4) = X(1);
        Y(3) = Y(2);     Y(4) = Y(1);
        Z(3) = Z(2)+0.5; Z(4) = Z(1)+0.5;
        
        for i = 1:length(SS{j})
            if ~isempty(SS{j}{i})
                F{j,i} = fill3(X,Y,Z+i,SS{j}{i},'facealpha',Opacity(1) );
                T{j,i} = text(mean(X),mean(Y),mean(Z+i),num2str(SS{j}{i}));
            end
            
        end
    end
    view(0,0)
end



if Dim == 3 % 3D
    for j=1:length(SS)
        X = PatchXYZ{j}.X;
        Y = PatchXYZ{j}.Y;
        Z = PatchXYZ{j}.Z;
        
        for i = 1:length(SS{j})
            if ~isempty(SS{j}{i})
                F{j,i} = fill3(X,Y,Z+(i-1),SS{j}{i},'facealpha',Opacity(1) );
                T{j,i} = text(mean(X),mean(Y),mean(Z+(i-1))+0.1,num2str(SS{j}{i}));
            end
            
        end
    end
    view(3)
end

h=colorbar;
set(h, 'ylim', [-90 90])
set(h, 'ticks', [-90:10:90])
b = uicontrol('Parent',Fig,'Style','slider','units','normalized','Position',[0.9,0.125,0.025,0.75],...
              'value',length(SS{j}), 'min',0, 'max',length(SS{j}),'Callback', @(es,ed) setOpacity(es,ed,Fig,F));
          