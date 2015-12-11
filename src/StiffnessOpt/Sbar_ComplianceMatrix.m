% % ===================================================================== %
%                       Custom Function (TM 21/05/2015)
%
% Sbar_ComplianceMatrix Returns the compliance matrices of a composite Mat.
% 
% SBar = Sbar_ComplianceMatrix (E1,E2,v12,G12,theta)
% Scalar: real values [E1,E2,v12,G12,theta] (theta must be in degree)
%
% Example: Sbar_ComplianceMatrix (82e9,4e9,0.25,2.8e9,30)
% % ===================================================================== %

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
 
 
R = [1 0 0;0 1 0; 0 0 2];

% rotated compliance matrix
SBar = R*(A^-1)*(R^-1)*S*A;

end
