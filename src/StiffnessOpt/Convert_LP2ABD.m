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
%                       Custom Function (TM 21/05/2015)
%
% LP2ABD Returns the stiffness matrices of the given lamination parameters.
% Ply thickness must be constant through the laminate.
% [A] and [D] matrices have been checked, need to verify [B] !!!
%
% [A,B,D] = LP2ABD (E1,E2,v12,G12,ply_t,LP,normalised)
% Scalar: real    values [E1,E2,v12,G12,ply_t,LP]
% Scalar: logical values [normalised] % if true, returns normalised stiffness
%
% LPs are defined as [V1A,V2A,V3A,V4A,V1B,V2B,V3B,V4B,V1D,V2D,V3D,V4D]
%
% Example:
% 1 - [A,B,D] = LP2ABD (181e9,10.3e9,0.28,7.17e9,0.000127*20,[0.24422  -0.00064598   -0.56232  0.0073377   0.010772 -0.0075078 -0.0017752  -0.012958   0.31104  0.00055537  -0.64423 0.0032003],true)
% 2 - [A,B,D] = LP2ABD (181e9,10.3e9,0.28,7.17e9,0.000127*16,SS2LP(0.000127,[42 -40  19 -38 -38 18 59 55 -47 -6 -47 32 37 -47 39 -24]' ),true)
%
% ===                                                                   ===

function [A,B,D,dAdLP,dBdLP,dDdLP] = Convert_LP2ABD (E1,E2,v12,G12,h,LP,NORMALISED)

V1A = LP(1);
V2A = LP(2);
V3A = LP(3);
V4A = LP(4);
V1B = LP(5);
V2B = LP(6);
V3B = LP(7);
V4B = LP(8);
V1D = LP(9);
V2D = LP(10);
V3D = LP(11);
V4D = LP(12);

% --
v21 = v12*E2/E1;
Q11 = E1/(1-v12*v21);
Q22 = E2/(1-v12*v21);
Q12 = v12*E2/(1-v12*v21);
Q66 = G12;

% Material invariants
U1 = 1/8*(3*Q11+3*Q22+2*Q12+4*Q66);
U2 = 1/2*(Q11-Q22);
U3 = 1/8*(Q11+Q22-2*Q12-4*Q66);
U4 = 1/8*(Q11+Q22+6*Q12-4*Q66);
U5 = 1/8*(Q11+Q22-2*Q12+4*Q66);

% Stiffness matrix
T0 = [U1, U4, 0 ;
    U4, U1, 0 ;
    0,  0,  U5];

T1 = [U2, 0,  0;
    0, -U2, 0;
    0,  0,  0];

T2 = [0,    0,    U2/2;
    0,    0,    U2/2;
    U2/2, U2/2,  0 ];

T3 = [U3, -U3, 0;
    -U3,  U3, 0;
    0,   0,  -U3];

T4 = [0,   0,  U3;
    0,   0, -U3;
    U3, -U3, 0 ];


dAdLP  = cell(3,3);
dBdLP  = cell(3,3);
dDdLP  = cell(3,3);


if NORMALISED
    A =  (T0 + T1*V1A + T2*V2A + T3*V3A + T4*V4A);
    B =  (     T1*V1B + T2*V2B + T3*V3B + T4*V4B);
    D =  (T0 + T1*V1D + T2*V2D + T3*V3D + T4*V4D);
    
    for i = 1:3
        for j = 1:3
            dAdLP{i,j} =  [T1(i,j), T2(i,j), T3(i,j), T4(i,j), 0, 0, 0, 0, 0, 0, 0, 0];
            dBdLP{i,j} =  [0, 0, 0, 0, T1(i,j), T2(i,j), T3(i,j), T4(i,j), 0, 0, 0, 0];
            dDdLP{i,j} =  [0, 0, 0, 0, 0, 0, 0, 0, T1(i,j), T2(i,j), T3(i,j), T4(i,j)];
        end
    end
else
    A =  h      * (T0 + T1*V1A + T2*V2A + T3*V3A + T4*V4A);
    B =  h^2/4  * (     T1*V1B + T2*V2B + T3*V3B + T4*V4B);
    D =  h^3/12 * (T0 + T1*V1D + T2*V2D + T3*V3D + T4*V4D);
    
    for i = 1:3
        for j = 1:3
            dAdLP{i,j} =  h      *[T1(i,j), T2(i,j), T3(i,j), T4(i,j), 0, 0, 0, 0, 0, 0, 0, 0];
            dBdLP{i,j} =  h^2/4  *[0, 0, 0, 0, T1(i,j), T2(i,j), T3(i,j), T4(i,j), 0, 0, 0, 0];
            dDdLP{i,j} =  h^3/12 *[0, 0, 0, 0, 0, 0, 0, 0, T1(i,j), T2(i,j), T3(i,j), T4(i,j)];
        end
    end
end



end