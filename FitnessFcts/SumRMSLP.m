function Fitness = SumRMSLP(LP,Objectives)

LP2Match    = reshape(cell2mat(Objectives.Table(2:end,3)),12,size(Objectives.Table,1)-1);
ScalingCoef = reshape(cell2mat(Objectives.Table(2:end,4)),12,size(Objectives.Table,1)-1);

localFit = zeros(size(LP2Match,2),1);
for ilam = 1 : size(LP2Match,2)
    localFit(ilam) = rms( (LP2Match(:,ilam) - LP(:,ilam)).*ScalingCoef(:,ilam) );
end
Fitness = sum(localFit);

end