function [Fitness,output] = EvaluationFct(SS,Parameters)

Nplies  = cellfun(@length,SS);
tplies  = Nplies*Parameters.ply_t;
Fitness = sum(tplies);                    

% ---

% Buckling constraints for all 18 panels
OptConstraint.FoS = 1;
Constraints = zeros(18,1);  

NORMALISED = false;
rowIndex   = 0;
for i = 1 : 18

    ply_angle = SS{i};
    [~,~,D] = Convert_SS2ABD (Parameters.E1,Parameters.E2,Parameters.v12,Parameters.G12,Parameters.ply_t,ply_angle,NORMALISED);

    a = Parameters.Dim{Parameters.PanelDim(i)}(1);
    b = Parameters.Dim{Parameters.PanelDim(i)}(2);
    
    for m = 1 : Parameters.mMax
        for n = 1 : Parameters.nMax
            
            rowIndex = rowIndex + 1;
            den      = ((m/a)^2*Parameters.Nx(i) + Parameters.Ny(i)*(n/b)^2);
            Constraints(rowIndex) = OptConstraint.FoS*1.0 - pi^2* ((D(1,1)*(m/a)^4 + 2*(D(1,2)+2*D(3,3))*(m/a)^2*(n/b)^2 + D(2,2)*(n/b)^4)) / den;  % should be 1 by default instead of 1.1
            
        end
    end
end

% ---
output.BucklingFactor = Constraints;
output.NViolatedConst = length(find(Constraints>0));

Fitness = Fitness + output.NViolatedConst *2;

end