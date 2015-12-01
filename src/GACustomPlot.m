function state = GACustomPlot(options,state,flag,PlotInterval,SaveInterval)

persistent Scores fileID                        
    
if state.Generation == 0
    Scores          = zeros(options.Generations,4);
    
    if ~isempty(SaveInterval)
        delete ('Results.txt');
        fileID = fopen('Results.txt','w');
        fprintf(fileID,'%4s \t %6s \t %6s \t %6s \n','Gen#','MinFit','AvgFit','MaxFit');
    end
    
    if ~isempty(PlotInterval)
        figure(1)
        hold on
    end
end

Scores(1+state.Generation,:) = [state.Generation min(state.Score) mean(state.Score) max(state.Score)];
if ~isempty(PlotInterval) && rem(state.Generation,PlotInterval) == 0
    figure(1)
    plot(Scores(1+state.Generation,1),Scores(1+state.Generation,2),'black x')
    plot(Scores(1+state.Generation,1),Scores(1+state.Generation,3),'blue +')
    plot(Scores(1+state.Generation,1),Scores(1+state.Generation,4),'red o')
    if state.Generation == 0
        legend('Best Fit','Avg Fit','Max Fit')
    end
end


if ~isempty(SaveInterval) && rem(state.Generation,SaveInterval) == 0
    fprintf(fileID,'%i \t %1.6f \t %1.6f \t %1.6f \n',Scores(1+state.Generation,:));
end

end