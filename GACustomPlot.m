function state = GACustomPlot(options,state,flag)
%GAPLOTBESTF Plots the best score and the mean score.
%   STATE = GAPLOTBESTF(OPTIONS,STATE,FLAG) plots the best score as well
%   as the mean of the scores.
%
%   Example:
%    Create an options structure that will use GAPLOTBESTF
%    as the plot function
%     options = gaoptimset('PlotFcns',@gaplotbestf);

%   Copyright 2003-2007 The MathWorks, Inc.

global Pop

Pop{state.Generation+1}=state.Population;

if size(state.Score,2) > 1
    title('Best Fitness Plot: not available','interp','none');
    return;
end

switch flag
    case 'init'
        hold on;
        set(gca,'xlim',[0,options.Generations]);
        xlabel('Generation','interp','none');
        ylabel('Fitness value','interp','none');
        plotBest = plot(state.Generation,min(state.Score),'.k');
        set(plotBest,'Tag','gaplotbestf');
        plotMean = plot(state.Generation,meanf(state.Score),'.b');
        set(plotMean,'Tag','gaplotmean');
        title(['Best: ',' Mean: '],'interp','none')
    case 'iter'
        best = min(state.Score);
        m    = meanf(state.Score);
        plotBest = findobj(get(gca,'Children'),'Tag','gaplotbestf');
        plotMean = findobj(get(gca,'Children'),'Tag','gaplotmean');
        newX = [get(plotBest,'Xdata') state.Generation];
        newY = [get(plotBest,'Ydata') best];
        set(plotBest,'Xdata',newX, 'Ydata',newY);
        newY = [get(plotMean,'Ydata') m];
        set(plotMean,'Xdata',newX, 'Ydata',newY);
        set(get(gca,'Title'),'String',sprintf('Best: %g Mean: %g',best,m));
    case 'done'
        LegnD = legend('Best fitness','Mean fitness');
        set(LegnD,'FontSize',8);
        hold off;
end

%------------------------------------------------
function m = meanf(x)
nans = isnan(x);
x(nans) = 0;
n = sum(~nans);
n(n==0) = NaN; % prevent divideByZero warnings
% Sum up non-NaNs, and divide by the number of non-NaNs.
m = sum(x) ./ n;

