function Fitness = SumNormLP(LP,Nplies,LP2Match,NPliesIni,IndexLp,ImportanceFactor,Constraints)

localFit = zeros(length(NPliesIni),1);
for ilam = 1 : length(Nplies)
    localFit(ilam) = norm(LP2Match(IndexLp,ilam) - LP(IndexLp,ilam));
end
Fitness = sum(localFit.*ImportanceFactor');


Fitness = Fitness ...
        + (Constraints.alpha)*sum(Nplies)*Constraints.ply_t;        % include weight penalty

end