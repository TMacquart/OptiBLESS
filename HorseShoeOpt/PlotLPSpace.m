% =====                                                              ==== %
%                      Plots lamination parameters space                  %
%                                                                         %
%  Terence Macquart Copyright (01/04/2015 )                               %
% LP1 and LP2 are the lamination parameters space you want to plot
% =====                                                              ==== %


function [] = PlotLPSpace(LP1,LP2,FigNumber)

theta  = -90:1:90;
Ntheta = length(theta);
V1_A   = zeros(Ntheta,1);
V2_A   = zeros(Ntheta,1);
V3_A   = zeros(Ntheta,1);
V4_A   = zeros(Ntheta,1);

for itheta = 1:Ntheta
    V1_A(itheta) = cosd(2*theta(itheta));
    V2_A(itheta) = sind(2*theta(itheta));
    V3_A(itheta) = cosd(4*theta(itheta));
    V4_A(itheta) = sind(4*theta(itheta));
end

if strcmp(LP1,'V1A'),
    LP1 = V1_A;
end
if strcmp(LP1,'V2A'),
    LP1 = V2_A;
end
if strcmp(LP1,'V3A'),
    LP1 = V3_A;
end
if strcmp(LP1,'V4A'),
    LP1 = V4_A;
end

if strcmp(LP2,'V1A'),
    LP2 = V1_A;
end
if strcmp(LP2,'V2A'),
    LP2 = V2_A;
end
if strcmp(LP2,'V3A'),
    LP2 = V3_A;
end
if strcmp(LP2,'V4A'),
    LP2 = V4_A;
end

figure(FigNumber); hold on;
plot(LP1,LP2,'blue --+')

end