function [NthetaVar,NbalVar,N10percentVar] = AttributeDesignVariable(MaxNplies,Constraints)

if Constraints.Vector(2)  % 10% rule explicit constrait handling genotype
    N10percentVar  = ceil(floor(0.4*MaxNplies)/2)*2;
else
    N10percentVar = 0;
end

NthetaVar = MaxNplies - N10percentVar;

if Constraints.Balanced
    if rem(NthetaVar,2)~=0, keyboard; end
    NthetaVar = NthetaVar/2;
    NbalVar   = NthetaVar;                                           % Max number of balanced angles
else
    NbalVar = 0;
end

if Constraints.Sym
    N10percentVar = N10percentVar/2;
    NthetaVar     = NthetaVar/2;
    NbalVar       = NbalVar/2;
end

RemovedN = N10percentVar - floor(N10percentVar/4)*4;

if rem(RemovedN,2)~=0
    keyboard
end
N10percentVar = N10percentVar-RemovedN;
NthetaVar     = NthetaVar + RemovedN/2;
NbalVar       = NbalVar + RemovedN/2;

NtotalPly = (NthetaVar+NbalVar+N10percentVar);

if Constraints.Sym && (MaxNplies/2)~= NtotalPly
    keyboard
end
if ~Constraints.Sym && (MaxNplies)~= NtotalPly
    keyboard
end
end