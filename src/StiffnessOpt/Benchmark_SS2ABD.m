% Benchmark SS2ABD

clear all; clc; format short g; format compact; close all force;

E1  = 181e9;
E2  = 10.3e9;
G12 = 7.170e9;
v12 = 0.28;
ply_t = 0.000127;
NORMALISED = true;
% ply_angle = [0 13.68 -12.10 -16.64 7.13 -15.50 17.46 5.82 5.97 -12.63];
% ply_angle = [42 -40 19 -38 -38 18 59 55 -47 -6 -47 32 37 -47 39 -24]; % (2)
ply_angle = fliplr([0 -6 -84 -5 42 4 -1 5 -72 -22 65 -84 5 -14 5 4]); % (22)
[A,B,D] = SS2ABD (E1,E2,v12,G12,ply_t,ply_angle,NORMALISED);

% E1 = 230e9;
% E2 = 6.6e9;
% G12 = 4.8e9;
% v12 = 0.25;
% ply_t = 0.5;
% ply_angle = [90 30 50]
% NORMALISED = false;
% [A,B,D] = SS2ABD (E1,E2,v12,G12,ply_t,ply_angle,NORMALISED);

% E1 = 82e9;
% E2 = 4e9;
% G12 = 2.8e9;
% v12 = 0.25;
% ply_t = 0.5;
% ply_angle = [30]
% NORMALISED = false;
% [A,B,D] = SS2ABD (E1,E2,v12,G12,ply_t,ply_angle,NORMALISED);

A = A*1e-9;
B = B*1e-9;
D = D*1e-9;

display(A)
display(B)
display(D)
