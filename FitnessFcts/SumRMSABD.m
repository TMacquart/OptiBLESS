function Fitness = SumRMSABD(A,B,D,Objectives)

% keyboard 

% IndexLP     = Objectives.IndexLP;
% LP2Match    = reshape(cell2mat(Objectives.Table(2:end,3)),12,size(Objectives.Table,1)-1);
% ScalingCoef = cell2mat(Objectives.Table(2:end,4));

localFit = zeros(size(A,2),1);
for ilam = 1 : size(A,2)
    localFit(ilam) = rms(A{ilam}(Objectives.IndexAStiff) - Objectives.A{ilam}(Objectives.IndexAStiff)) ...
                   + rms(B{ilam}(Objectives.IndexBStiff) - Objectives.B{ilam}(Objectives.IndexBStiff)) ...
                   + rms(D{ilam}(Objectives.IndexDStiff) - Objectives.D{ilam}(Objectives.IndexDStiff));
end
Fitness = sum(localFit); %sum(localFit.*ScalingCoef);

end