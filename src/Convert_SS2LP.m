% =====                                                              ====
% SS2LP returns the lamination parameters of a given stacking sequence.
% Ply thickness is assumed constant through the laminate.
%
% [LP,V1A,V2A,V3A,V4A,V1B,V2B,V3B,V4B,V1D,V2D,V3D,V4D] = SS2LP(ply_angle,unit,ply_t,varargin)
% SS2LP Necessary Input
% 'ply_angle' is a vector of real values (degree by default)
% The ply angles are defined from bottom (1st element) to top (last element)
%
% SS2LP Optional Inputs
% 'unit' is a string specifying the angles units either 'deg' or 'rad'
% 'ply_t' is a real scalar value corresponding the ply thickness
%
% SS2LP Example:
% LP = SS2LP([42 -40  19 -38 -38 18 59 55 -47 -6 -47 32 37 -47 39 -24]','deg', 0.000127)
% =====                                                              ====


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


function [LP,V1A,V2A,V3A,V4A,V1B,V2B,V3B,V4B,V1D,V2D,V3D,V4D] = Convert_SS2LP(ply_angle,unit,ply_t,varargin)
%% Input Check


if ~isnumeric(ply_angle) || ~isreal(ply_angle) || ~isvector(ply_angle)
    error('Input angles (ply_angle) must be numeric and real vector')
end

if nargin > 1
    if strcmp(unit,'rad')
        ply_angle = rad2deg(ply_angle);
    else
        if ~strcmp(unit,'deg')
            error('Unallowed unit type used as input to SS2LP. Only accepts -- rad -- or  -- deg --')
        end
    end
end

if nargin == 3
    if ~isnumeric(ply_t) || ~isreal(ply_t)
        error('Input ply thickness (ply_t) must be numeric and real')
    end
end



%% input format
if size(ply_angle,1) == 1
    ply_angle = ply_angle';
end
ply_angle = ply_angle*pi/180;
Nplies    = length(ply_angle);

%% classic sum from -1/2 to 1/2
if nargin == 1 || nargin == 2
    z1 = (-0.5:1/Nplies:0.5-1/Nplies)';
    z2 = (-0.5+1/Nplies:1/Nplies:0.5)';
    
    Cos2Angle = cos(2*ply_angle);
    Sin2Angle = sin(2*ply_angle);
    Cos4Angle = cos(4*ply_angle);
    Sin4Angle = sin(4*ply_angle);
    
    DeltaZ  =  z2-z1;
    DeltaZ2 = (z2.^2-z1.^2);
    DeltaZ3 = (z2.^3-z1.^3);
    
    V1A  = sum(Cos2Angle.*DeltaZ);      % = sum(cos(2*ply_angle).*(z2-z1))
    V2A  = sum(Sin2Angle.*DeltaZ);      % = sum(sin(2*ply_angle).*(z2-z1));
    V3A  = sum(Cos4Angle.*DeltaZ);      % = sum(cos(4*ply_angle).*(z2-z1));
    V4A  = sum(Sin4Angle.*DeltaZ);      % = sum(sin(4*ply_angle).*(z2-z1));
    
    V1B  = 2*sum(Cos2Angle.*DeltaZ2);   % = 4*sum(1/2*cos(2*ply_angle).*(z2.^2-z1.^2));
    V2B  = 2*sum(Sin2Angle.*DeltaZ2);   % = 4*sum(1/2*sin(2*ply_angle).*(z2.^2-z1.^2));
    V3B  = 2*sum(Cos4Angle.*DeltaZ2);   % = 4*sum(1/2*cos(4*ply_angle).*(z2.^2-z1.^2));
    V4B  = 2*sum(Sin4Angle.*DeltaZ2);   % = 4*sum(1/2*sin(4*ply_angle).*(z2.^2-z1.^2));
    
    V1D  = 4*sum(Cos2Angle.*DeltaZ3);   % = 12*sum(1/3*cos(2*ply_angle).*(z2.^3-z1.^3));
    V2D  = 4*sum(Sin2Angle.*DeltaZ3);   % = 12*sum(1/3*sin(2*ply_angle).*(z2.^3-z1.^3));
    V3D  = 4*sum(Cos4Angle.*DeltaZ3);   % = 12*sum(1/3*cos(4*ply_angle).*(z2.^3-z1.^3));
    V4D  = 4*sum(Sin4Angle.*DeltaZ3);   % = 12*sum(1/3*sin(4*ply_angle).*(z2.^3-z1.^3));
end


%% sum from -N/2 to N/2 (not Used)
if 0 && nargin == 2
    Zi = - Nplies/2 + [0:Nplies];
    z1 = Zi(1:end-1)';
    z2 = Zi(2:end)';
    
    V1A = 1/Nplies * sum(cos(2*ply_angle));
    V2A = 1/Nplies * sum(sin(2*ply_angle));
    V3A = 1/Nplies * sum(cos(4*ply_angle));
    V4A = 1/Nplies * sum(sin(4*ply_angle));
    
    V1B = 2/(Nplies^2) * sum(cos(2*ply_angle).*(z2.^2-z1.^2));
    V2B = 2/(Nplies^2) * sum(sin(2*ply_angle).*(z2.^2-z1.^2));
    V3B = 2/(Nplies^2) * sum(cos(4*ply_angle).*(z2.^2-z1.^2));
    V4B = 2/(Nplies^2) * sum(sin(4*ply_angle).*(z2.^2-z1.^2));
    
    V1D = 4/(Nplies^3) * sum(cos(2*ply_angle).*(z2.^3-z1.^3));
    V2D = 4/(Nplies^3) * sum(sin(2*ply_angle).*(z2.^3-z1.^3));
    V3D = 4/(Nplies^3) * sum(cos(4*ply_angle).*(z2.^3-z1.^3));
    V4D = 4/(Nplies^3) * sum(sin(4*ply_angle).*(z2.^3-z1.^3));
end


%% sum form -h/2 to h/2
if nargin == 3
    z1 = zeros(Nplies,1);
    z2 = zeros(Nplies,1);
    h  = Nplies*ply_t;
    
    for i = 1 : Nplies
        z1(i) = -h/2  + ply_t * (i-1);
        z2(i) = z1(i) + ply_t ;
    end
    
    V1A = 1/Nplies  * sum(cos(2*ply_angle));
    V2A = 1/Nplies  * sum(sin(2*ply_angle));
    V3A = 1/Nplies  * sum(cos(4*ply_angle));
    V4A = 1/Nplies  * sum(sin(4*ply_angle));
    V1B = 4/h^2     * sum(1/2*cos(2*ply_angle).*(z2.^2-z1.^2));
    V2B = 4/h^2     * sum(1/2*sin(2*ply_angle).*(z2.^2-z1.^2));
    V3B = 4/h^2     * sum(1/2*cos(4*ply_angle).*(z2.^2-z1.^2));
    V4B = 4/h^2     * sum(1/2*sin(4*ply_angle).*(z2.^2-z1.^2));
    V1D = 12/h^3    * sum(1/3*cos(2*ply_angle).*(z2.^3-z1.^3));
    V2D = 12/h^3    * sum(1/3*sin(2*ply_angle).*(z2.^3-z1.^3));
    V3D = 12/h^3    * sum(1/3*cos(4*ply_angle).*(z2.^3-z1.^3));
    V4D = 12/h^3    * sum(1/3*sin(4*ply_angle).*(z2.^3-z1.^3));
end


%%
LP  = [V1A, V2A, V3A, V4A, V1B, V2B, V3B, V4B, V1D, V2D, V3D, V4D]';

