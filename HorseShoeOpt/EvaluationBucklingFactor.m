function [SumBucklingFactors,Gradient] = EvaluationBucklingFactor(dv)

global Parameters OptConstraint

BucklingFactors = zeros(200,1);                 % inequality constraints PreAllocation (over-sized)
Jacobian    = zeros(200,length(dv));        % Jacobian of the inequality constraints PreAllocation (over-sized)
rowIndex    = 0;

% ---
if 1    % Buckling constraints for all 18 panels
    
    t = OptConstraint.FixedThick * Parameters.ply_tMax;
    
    for i = 1 : 18
        
        LP      = [zeros(1,8) dv(1 + (i-1)*2) 0 dv(2 + (i-1)*2) 0 ]';
        [~,~,D] = LP2ABD (Parameters.E1,Parameters.E2,Parameters.v12,Parameters.G12,t(i),LP,false);
        
        a = Parameters.Dim{Parameters.PanelDim(i)}(1);
        b = Parameters.Dim{Parameters.PanelDim(i)}(2);
        
        for m = 1 : Parameters.mMax
            for n = 1 : Parameters.nMax
                
                den = ((m/a)^2*Parameters.Nx(i) + Parameters.Ny(i)*(n/b)^2);
                rowIndex = rowIndex + 1;
                
                BucklingFactors(rowIndex) =   OptConstraint.FoS*1.0 - pi^2* ((D(1,1)*(m/a)^4 + 2*(D(1,2)+2*D(3,3))*(m/a)^2*(n/b)^2 + D(2,2)*(n/b)^4)) / den;  % should be 1 by default instead of 1.1
                
                if OptConstraint.Analytical % --- Gradient
                    %  D =  (dv(i)*Coeff)^3/12 * (T0 + T1*V1D + T3*V3D);
                    
                    dlbd_dD11 = - pi^2 /den *(m/a)^4;
                    dlbd_dD12 = - pi^2 /den *2*(m/a)^2*(n/b)^2;
                    dlbd_dD22 = - pi^2 /den *(n/b)^4;
                    dlbd_dD66 = - pi^2 /den *4*(m/a)^2*(n/b)^2;
                    
                    %                     dD11_dt   = Parameters.ply_tMax^3 * dv(i)^2/4 *(Parameters.T0(1,1) + Parameters.T1(1,1)*dv(19 + (i-1)*2) + Parameters.T3(1,1)*dv(20 + (i-1)*2));
                    dD11_dV1D = t(i)^3/12*(Parameters.T1(1,1));
                    dD11_dV3D = t(i)^3/12*(Parameters.T3(1,1));
                    
                    %                     dD12_dt   = Parameters.ply_tMax^3 * dv(i)^2/4 *(Parameters.T0(1,2) + Parameters.T1(1,2)*dv(19 + (i-1)*2) + Parameters.T3(1,2)*dv(20 + (i-1)*2));
                    dD12_dV1D = t(i)^3/12*(Parameters.T1(1,2));
                    dD12_dV3D = t(i)^3/12*(Parameters.T3(1,2));
                    
                    %                     dD22_dt   = Parameters.ply_tMax^3 * dv(i)^2/4 *(Parameters.T0(2,2) + Parameters.T1(2,2)*dv(19 + (i-1)*2) + Parameters.T3(2,2)*dv(20 + (i-1)*2));
                    dD22_dV1D = t(i)^3/12*(Parameters.T1(2,2));
                    dD22_dV3D = t(i)^3/12*(Parameters.T3(2,2));
                    
                    %                     dD66_dt   = Parameters.ply_tMax^3 * dv(i)^2/4 *(Parameters.T0(3,3) + Parameters.T1(3,3)*dv(19 + (i-1)*2) + Parameters.T3(3,3)*dv(20 + (i-1)*2));
                    dD66_dV1D = t(i)^3/12*(Parameters.T1(3,3));
                    dD66_dV3D = t(i)^3/12*(Parameters.T3(3,3));
                    
                    %                     Jacobian(rowIndex,i)              = dlbd_dD11 * dD11_dt    + dlbd_dD12*dD12_dt   + dlbd_dD22*dD22_dt   + dlbd_dD66*dD66_dt;
                    Jacobian(rowIndex,(1 + (i-1)*2)) = dlbd_dD11 * dD11_dV1D  + dlbd_dD12*dD12_dV1D + dlbd_dD22*dD22_dV1D + dlbd_dD66*dD66_dV1D;
                    Jacobian(rowIndex,(2 + (i-1)*2)) = dlbd_dD11 * dD11_dV3D  + dlbd_dD12*dD12_dV3D + dlbd_dD22*dD22_dV3D + dlbd_dD66*dD66_dV3D;
                end
            end
        end
    end
end
% ---

BucklingFactors(rowIndex+1:end) = [];
Jacobian(rowIndex+1:end,:)      = [];

SumBucklingFactors = sum(BucklingFactors);
Jacobian           = sum (Jacobian,1);
Gradient           = Jacobian';
end