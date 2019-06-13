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
% =====                                                              ==== 
%                   Custom GA output Function.
%
% If GAoptions.SaveInterval is not empty this function is called during GA 
% and population scores are saved into an external text file (see main
% folder Results.txt). 
%
%
% [state,options,optchanged] = GACustomOutput(options,state,flag,SaveInterval)
%  
% =====                                                              ==== 

function [state,options,optchanged] = GACustomOutput(options,state,flag,SaveInterval,ObjType)

optchanged = false;
persistent fileID ButtonHandle     
% fileID contains the indentier to the Result.txt file
    

if state.Generation == 0

    ButtonHandle = findall(0,'String','Stop'); % find the handle of the push button to stop GA
    set(ButtonHandle,'Callback',@SetValue)     % Add stop callback to the button
    
    if ~isempty(SaveInterval)
        fclose('all');
        fileID = fopen('Results.txt','wt');
        fprintf(fileID,'%4s \t %6s \t %6s \t %6s \n','Gen#','MinFit','AvgFit','MaxFit');
    end
end

Scores = [state.Generation min(state.Score) mean(state.Score) max(state.Score)];
if ~isempty(SaveInterval) && rem(state.Generation,SaveInterval) == 0
    if strcmp(ObjType,'lp')
        fprintf(fileID,'%3.0f \t %1.6f \t %1.6f \t %1.6f \n',Scores);
        
    elseif strcmp(ObjType,'ABD')
        fprintf(fileID,'%3.0f \t %1.5E \t %1.5E \t %1.5E \n',Scores);
        
    else % User Function 
        fprintf(fileID,'%3.0f \t %1.5E \t %1.5E \t %1.5E \n',Scores); % Change here if needed for user function output precision in txt file.
        
    end
end


if strcmp(get(ButtonHandle,'String'),'HasBeenPushed')
    state.StopFlag = 'true';
end

end

function SetValue(hObject,callbackdata)
   set(hObject,'String','HasBeenPushed')
end