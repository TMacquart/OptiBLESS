% Benchmark ABD2LP

clear all; clc; format short g; format compact; close all;

Lp2Match = [
% LP2Match1 LP2Match2  LP2Match3
    0.1821	 1	 0.3000   % V1A
   -0.3643	 0	-0.1732   % V2A
    0.0667	 1	 0.1000   % V3A
   -0.1155	 0	-0.1732   % V4A
    0.0000	 0	 0.1000   % V1B
    0.0000	 0	 -0.500   % V2B
    0.0000	 0	 0.2000   % V3B
    0.0000	 0	 0.0000   % V4B
    0.1699	 1	-0.0741   % V1D
   -0.2584	 0	 0.1131   % V2D
    0.2261	 1	 0.4120   % V3D
   -0.3177	 0	-0.3811]; % V4D


E1  = 13.0e9;
E2  = 72.0e9;
G12 = 26.9e9;
v12 = 0.33;
h   = 0.000127 * 10;

[A2Match,B2Match,D2Match] = LP2ABD (E1,E2,v12,G12,h,Lp2Match(:,2),true);
[LP,Abar,Bbar,Dbar]       = ABD2LP (E1,E2,v12,G12,h,A2Match,B2Match,D2Match,true);
