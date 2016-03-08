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
% ====                                                                 ====
%                       Custom Function (TM 21/05/2015)
%
% Sbar_ComplianceMatrix Returns the compliance matrices of a composite Mat.
% 
% SBar = Sbar_ComplianceMatrix (E1,E2,v12,G12,theta)
% Scalar: real values [E1,E2,v12,G12,theta] (theta must be in degree)
%
% Example: Sbar_ComplianceMatrix (82e9,4e9,0.25,2.8e9,30)
%
% ====                                                                 ====

function SBar = Sbar_ComplianceMatrix (E1,E2,v12,G12,theta)

v21    = E2*v12/E1;
thetha = theta*pi/180;

S = [1/E1     -v21/E2      0;            % strain = S * stress
     -v12/E1   1/E2        0;
     0          0       1/G12];

% Rotation matrices
c_2 = cos(thetha)^2;
s_2 = sin(thetha)^2;
sc  = cos(thetha)*sin(thetha);

A = [c_2   s_2      2*sc;
     s_2   c_2     -2*sc;
     -sc   sc     c_2-s_2];
 
 
R    = [1 0 0;0 1 0; 0 0 2];
Rinv = [1 0 0;0 1 0; 0 0 0.5];

% rotated compliance matrix
SBar = R*(A^-1)*(Rinv)*S*A;

end
