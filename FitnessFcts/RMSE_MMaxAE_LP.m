function Fitness = RMSE_MMaxAE_LP(LP,Objectives)

LP2Match    = reshape(cell2mat(Objectives.Table(2:end,3)),12,size(Objectives.Table,1)-1);
ScalingCoef = reshape(cell2mat(Objectives.Table(2:end,4)),12,size(Objectives.Table,1)-1);

Nlam = size(LP2Match,2);
localFit = zeros(Nlam,1);
for ilam = 1 : Nlam
    Error = (LP2Match(:,ilam) - LP(:,ilam)).*ScalingCoef(:,ilam);
    localFit(ilam) = rms(Error) + max(abs(Error));
end
Fitness = mean(localFit);

end