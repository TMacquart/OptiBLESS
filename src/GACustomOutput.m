function [state,options,optchanged] = GACustomOutput(options,state,flag,SaveInterval)

optchanged = false;
persistent Scores fileID ButtonHandle                    
    
if state.Generation == 0
    Scores = zeros(options.Generations,4);
    
    ButtonHandle = findall(0,'String','Stop'); % find the handle of the push button to stop GA
    set(ButtonHandle,'Callback',@SetValue)
    
    if ~isempty(SaveInterval)
        fclose('all');
        fileID = fopen('Results.txt','wt');
        fprintf(fileID,'%4s \t %6s \t %6s \t %6s \n','Gen#','MinFit','AvgFit','MaxFit');
    end
end

Scores(1+state.Generation,:) = [state.Generation min(state.Score) mean(state.Score) max(state.Score)];
if ~isempty(SaveInterval) && rem(state.Generation,SaveInterval) == 0
    fprintf(fileID,'%3.0f \t %1.6f \t %1.6f \t %1.6f \n',Scores(1+state.Generation,:));
end


if strcmp(get(ButtonHandle,'String'),'HasBeenPushed')
    state.StopFlag = 'true';
end

end

function SetValue(hObject,callbackdata)
   set(hObject,'String','HasBeenPushed')
end