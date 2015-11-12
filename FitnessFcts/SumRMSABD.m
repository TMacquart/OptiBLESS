function Fitness = SumRMSABD(A,B,D,Objectives)

% keyboard 

% IndexLP     = Objectives.IndexLP;
% LP2Match    = reshape(cell2mat(Objectives.Table(2:end,3)),12,size(Objectives.Table,1)-1);
% ScalingCoef = cell2mat(Objectives.Table(2:end,4));

localFit = zeros(size(A,2),1);
for ilam = 1 : size(A,2)
    localFit(ilam) = rms(A{ilam}(:) - Objectives.A{ilam}(:)) + rms(B{ilam}(:) - Objectives.B{ilam}(:)) + rms(D{ilam}(:) - Objectives.D{ilam}(:));
end
Fitness = sum(localFit); %sum(localFit.*ScalingCoef);

end