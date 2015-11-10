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


% % ===================================================================== %
%                       Custom Function (TM 21/05/2015)
% 
% SS2LP Returns the lamination parameters of a given stacking sequence.
% Ply thickness is assumed constant through the laminate.
% 
% SS2LP(ply_t,ply_angle)
% [LP,V1A,V2A,V3A,V4A,V1B,V2B,V3B,V4B,V1D,V2D,V3D,V4D] = SS2LP(ply_t,ply_angle)
% Scalar: real values [ply_t]
% Vector: real values [ply_angle] (must be in degree)
% 
% ply_angles are defined from bottom (1st element) to top (last element)
%
% Example: LP = SS2LP(0.000127,[42 -40  19 -38 -38 18 59 55 -47 -6 -47 32 37 -47 39 -24]' )
% % ===================================================================== %

function [LP,V1A,V2A,V3A,V4A,V1B,V2B,V3B,V4B,V1D,V2D,V3D,V4D] = SS2LP(ply_t,ply_angle)

Nplies    = length(ply_angle);
h         = Nplies*ply_t;
lamparvec = zeros(12,1);
z1        = zeros(Nplies,1);
z2        = zeros(Nplies,1);

for i = 1 : Nplies
    z1(i) = -h/2  + ply_t * (i-1);
    z2(i) = z1(i) + ply_t ;
end

lamparvec(1)  = sum(cosd(2*ply_angle));
lamparvec(2)  = sum(sind(2*ply_angle));
lamparvec(3)  = sum(cosd(4*ply_angle));
lamparvec(4)  = sum(sind(4*ply_angle));
lamparvec(5)  = sum(1/2*cosd(2*ply_angle).*(z2.^2-z1.^2));
lamparvec(6)  = sum(1/2*sind(2*ply_angle).*(z2.^2-z1.^2));
lamparvec(7)  = sum(1/2*cosd(4*ply_angle).*(z2.^2-z1.^2));
lamparvec(8)  = sum(1/2*sind(4*ply_angle).*(z2.^2-z1.^2));
lamparvec(9)  = sum(1/3*cosd(2*ply_angle).*(z2.^3-z1.^3));
lamparvec(10) = sum(1/3*sind(2*ply_angle).*(z2.^3-z1.^3));
lamparvec(11) = sum(1/3*cosd(4*ply_angle).*(z2.^3-z1.^3));
lamparvec(12) = sum(1/3*sind(4*ply_angle).*(z2.^3-z1.^3));

% - 12 lamination parameters
V1A = 1/Nplies  * lamparvec(1); 
V2A = 1/Nplies  * lamparvec(2);
V3A = 1/Nplies  * lamparvec(3); 
V4A = 1/Nplies  * lamparvec(4); 
V1B = 4/h^2     * lamparvec(5);
V2B = 4/h^2     * lamparvec(6);
V3B = 4/h^2     * lamparvec(7);
V4B = 4/h^2     * lamparvec(8);
V1D = 12/h^3    * lamparvec(9);
V2D = 12/h^3    * lamparvec(10);
V3D = 12/h^3    * lamparvec(11);
V4D = 12/h^3    * lamparvec(12);

LP  = [V1A, V2A, V3A, V4A, V1B, V2B, V3B, V4B, V1D, V2D, V3D, V4D]';

end