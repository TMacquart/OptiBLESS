function Fitness = SumRMSABD(A,B,D,Objectives)

Nlam = size(A,2);
localFit = zeros(Nlam,1);
for ilam = 1 : Nlam
    
    AScaling = Objectives.Table{ilam+1,6};
    BSacling = Objectives.Table{ilam+1,7};
    DScaling = Objectives.Table{ilam+1,8};
    
    localFit(ilam) = rms( AScaling(:).*(A{ilam}(:) - Objectives.Table{ilam+1,3}(:)) ) ...
                   + rms( BSacling(:).*(B{ilam}(:) - Objectives.Table{ilam+1,4}(:)) ) ...
                   + rms( DScaling(:).*(D{ilam}(:) - Objectives.Table{ilam+1,5}(:)) );
end
Fitness = sum(localFit)/Nlam; %

end