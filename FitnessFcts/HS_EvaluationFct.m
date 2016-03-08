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
%                        HorseShoe Fitness function                        
% ===                                                                   === 

function [Fitness,output] = HS_EvaluationFct(SSTable,Parameters)


for i=1:size(SSTable,1)
    Nplies(i)  = length(cell2mat(SSTable(i,:)));
end
% Nplies = [34 28 22 18 16 22 18 26 38 36 30 28 22 18 26 38 18 22] % Adams et al.
% Nplies = [34 28 22 19 16 22 19 26 38 35 30 28 22 19 26 32 19 24] % Irisarri et al.
tplies  = Nplies*Parameters.ply_t; % in inch
                   
Vol = zeros(18,1);
for i=1:18
    Vol(i,1) = tplies(i)*Parameters.Dim{Parameters.PanelDim(i)}(1)*Parameters.Dim{Parameters.PanelDim(i)}(2);
end
Fitness = sum(Vol)*Parameters.rho;


% ---
% Buckling constraints for all 18 panels
OptConstraint.FoS = 1;
Constraints = zeros(18,1);  

NORMALISED = false;
rowIndex   = 0;
for i = 1 : 18

    ply_angle = cell2mat(SSTable(i,:));
    [~,~,D]   = Convert_SS2ABD (Parameters.E1,Parameters.E2,Parameters.v12,Parameters.G12,Parameters.ply_t,ply_angle,NORMALISED);

    a = Parameters.Dim{Parameters.PanelDim(i)}(1); % longitudinal panel length
    b = Parameters.Dim{Parameters.PanelDim(i)}(2); % transversal  panel length
    
    for m = 1 : Parameters.mMax
        for n = 1 : Parameters.nMax
            
            rowIndex = rowIndex + 1;
            den      = (Parameters.Nx(i)*(m/a)^2 + Parameters.Ny(i)*(n/b)^2);
            Constraints(rowIndex) = OptConstraint.FoS*1.0 - pi^2* ((D(1,1)*((m/a)^4) + 2*(D(1,2)+2*D(3,3))*((m/a)^2)*((n/b)^2) + D(2,2)*(n/b)^4)) / den;  % should be 1 by default instead of 1.1
            
        end
    end
end

% ---
output.BucklingFactor = Constraints;
output.NViolatedConst = length(find(Constraints>0));

Fitness = Fitness + output.NViolatedConst * 25;

end