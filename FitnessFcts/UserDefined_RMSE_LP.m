% =====                                                                ====
%     Example of user defined function using RMSE_LP as fitness function 
%     Average of the Root Mean Square Error between lamination paramters
% 
% In addition to Fitness, user defined fitness function must also return an
% output structure which must at least contain 1 field corresponding to the
% number of violated constraints (i.e. output.NViolatedConst). This is used
% during initial population generation to check the feasibility of initial
% individuals.
% 
% [Fitness,output]  = UserDefined_RMSE_LP(LP,Objectives)
% =====                                                                ====

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

function [Fitness,output] = UserDefined_RMSE_LP(LP,Objectives)

LP2Match    = reshape(cell2mat(Objectives.Table(2:end,3)),12,size(Objectives.Table,1)-1);
ScalingCoef = reshape(cell2mat(Objectives.Table(2:end,4)),12,size(Objectives.Table,1)-1);

Nlam = size(LP2Match,2);
localFit = zeros(Nlam,1);
for ilam = 1 : Nlam
    Error = (LP2Match(:,ilam) - LP(:,ilam)).*ScalingCoef(:,ilam);
    localFit(ilam) = rms(Error);
end
Fitness = sum(localFit)/Nlam;

output.NViolatedConst = 0; % used only for Initial Population
end