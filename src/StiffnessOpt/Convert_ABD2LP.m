function [LPOpt,AOpt,BOpt,DOpt] = Convert_ABD2LP (E1,E2,v12,G12,h,A2Match,B2Match,D2Match,NORMALISED)


optionsOpt  = optimset('Algorithm','interior-point','GradConstr','on','GradObj','off','DerivativeCheck','on','DiffMaxChange',1e-5,...
    'DiffMinChange',1e-6,'TolFun',1e-8,'TolX',1e-6,'TolCon',1e-6,'PlotFcns',{},... % @optimplotfval @optimplotconstrviolation
    'FinDiffType','central','MaxFunEvals',1e6,'MaxIter',1e5);

EvaluationFct = @(LP) WrapperLP2ABD(E1,E2,v12,G12,h,LP,NORMALISED,A2Match,B2Match,D2Match);

LPini = zeros(12,1);
LPmin = -1*ones(12,1);
LPmax = ones(12,1);


display('Rapid matching in progress, please wait ...')
LPOpt = fmincon(EvaluationFct,LPini,[],[],[],[],LPmin,LPmax,[],optionsOpt);

if isempty(find(LPConstraints(LPOpt)>0,1))
    display('Rapid matching succeeded. The following results have been achieved:')
else
    display('Rapid matching failed. Constraints matching in progress...')
    LPOpt = fmincon(EvaluationFct,LPOpt,[],[],[],[],LPmin,LPmax,@LPConstraints,optionsOpt);
end

[AOpt,BOpt,DOpt] = Convert_LP2ABD (E1,E2,v12,G12,h,LPOpt,NORMALISED);
display('A Matrix Matching Error')
display([A2Match-AOpt])
display('B Matrix Matching Error')
display([B2Match-BOpt])
display('D Matrix Matching Error')
display([D2Match-DOpt])




    function [Obj,dObj] = WrapperLP2ABD(E1,E2,v12,G12,h,LP,NORMALISED,A2Match,B2Match,D2Match)
        
        [A,B,D,dAdLP,dBdLP,dDdLP] = Convert_LP2ABD (E1,E2,v12,G12,h,LP,NORMALISED);

        Obj = (A-A2Match).^2 +(B-B2Match).^2 +(D-D2Match).^2;
        Obj = sum(Obj(:));

        dObjdA = +2*(A-A2Match);
        dObjdA = dObjdA(:);
        
        dObjdB = +2*(B-B2Match);
        dObjdB = dObjdB(:);
        
        dObjdD = +2*(D-D2Match);
        dObjdD = dObjdD(:);
        
        dObj = 0;
        for jj = 1 : 9
            dObj  = dObj + dObjdA(jj)*dAdLP{jj} + dObjdB(jj)*dBdLP{jj} + dObjdD(jj)*dDdLP{jj};
        end
    end

    function [Constraints,Ceq,gradc,gradceq]  = LPConstraints(LP)
        
        Constraints = [];
        gradc   = [];
        Ceq     = [];
        gradceq = [];
        
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
        
        Constraints = [Constraints ; 1-2*V1A^2*(1-V3A)-2*V2A^2*(1+V3A)-V3A^2-V4A^2+4*V1A*V2A*V4A];
        Constraints = [Constraints ; 1- V1A^2 - V2A^2];
        Constraints = [Constraints ; 1-2*V1B^2*(1-V3B)-2*V2B^2*(1+V3B)-V3B^2-V4B^2+4*V1B*V2B*V4B];
        Constraints = [Constraints ; 1- V1B^2 - V2B^2];
        Constraints = [Constraints ; 1-2*V1D^2*(1-V3D)-2*V2D^2*(1+V3D)-V3D^2-V4D^2+4*V1D*V2D*V4D];
        Constraints = [Constraints ; 1- V1D^2 - V2D^2];
        
        Constraints = - Constraints;
        
        gradc = [gradc; [-4*V1A*(1-V3A) + 4*V2A*V4A, -4*V2A*(1+V3A)+4*V1A*V4A, 2*V1A^2-2*V2A^2-2*V3A, -2*V4A+4*V1A*V2A, zeros(1,8)]];
        gradc = [gradc; [-2*V1A,-2*V2A,0,0,zeros(1,8)]];
        gradc = [gradc; [zeros(1,4),-4*V1B*(1-V3B)+4*V2B*V4B,-4*V2B*(1+V3B)+4*V1B*V4B,2*V1B^2-2*V2B^2-2*V3B,-2*V4B+4*V1B*V2B,zeros(1,4)]];
        gradc = [gradc; [zeros(1,4),-2*V1B,-2*V2B,0,0,zeros(1,4)]];
        gradc = [gradc; [zeros(1,8),-4*V1D*(1-V3D)+4*V2D*V4D,-4*V2D*(1+V3D)+4*V1D*V4D,2*V1D^2-2*V2D^2-2*V3D,-2*V4D+4*V1D*V2D]];
        gradc = [gradc; [zeros(1,8),-2*V1D,-2*V2D,0,0]];
        
        gradc = -gradc';
    end
end