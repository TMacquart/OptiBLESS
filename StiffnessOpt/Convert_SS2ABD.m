% % ===================================================================== %
%                       Custom Function (TM 21/05/2015)
%
% SS2ABD Returns the stiffness matrices of a given stacking sequence.
% Ply thickness must be constant through the laminate.
% [A] and [D] matrices have been checked, need to verify [B] !!!
%
% [A,B,D] = SS2ABD (E1,E2,v12,G12,ply_t,ply_angle,NORMALISED)
% Scalar: real values [E1,E2,v12,G12,ply_t]
% Scalar: logical     [NORMALISED] % set to true for NORMALISED stiffness
% Vector: real values [ply_angle] (must be in degree)
%
% ply_angles are defined from bottom (1st element) to top (last element)
%
% Example: [A,B,D] = SS2ABD(181e9,10.3e9,0.28,7.17e9,0.000127,[42 -40  19 -38 -38 18 59 55 -47 -6 -47 32 37 -47 39 -24]',true )
% % ===================================================================== %


function [A,B,D] = Convert_SS2ABD (E1,E2,v12,G12,ply_t,ply_angle,NORMALISED)

ply_t      = ply_t*ones(length(ply_angle),1);
h          = sum(ply_t);
Nplies     = length(ply_angle);
lam_center = h/2;               % laminate center

SBar  = cell(Nplies,1);
DBar  = cell(Nplies,1);
A     = zeros(3,3);
B     = zeros(3,3);
D     = zeros(3,3);
z     = zeros(Nplies+1,1);

for i = 1:Nplies+1
    if      i == 1,         z(i) = -lam_center;
    elseif  i == Nplies+1,  z(i) = lam_center;
    else                    z(i) = -lam_center + sum(ply_t(1:i-1));   end
end


for i = 1:Nplies
    SBar{i} = Sbar_ComplianceMatrix(E1,E2,v12,G12,ply_angle(i)); % Sbar_matrix (E1,E2,v12,G12,theta)
    DBar{i} = SBar{i}^-1;
   
    A = A + DBar{i} * ply_t(i);
    B = B + 1/2 * DBar{i} * (z(i+1)^2 - z(i)^2);
    D = D + 1/3 * DBar{i} * (z(i+1)^3 - z(i)^3);
end;


if NORMALISED
    A = A/h;  
    B = 4/h^2 *B;
    D = 12/h^3*D;
end
end