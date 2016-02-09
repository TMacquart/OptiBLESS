function NStruct = AttributeNply(Nplies,Constraints,LamType)

MaxNplies = max(Nplies(:)); 

if strcmp(LamType,'Generic')
    NbalVar   = 0;
    NMidPlane = 0;
    if Constraints.Vector(2)==0 % no 10% rule
        N10percentVar = 0; 
        NthetaVar = MaxNplies;
    else
        N10percentVar = round((0.4*MaxNplies)/4)*4;
        if N10percentVar > MaxNplies
            keyboard % should not happen
        end
        NthetaVar = MaxNplies-N10percentVar;
    end
    NtotalPly = NthetaVar + NbalVar + N10percentVar;
end



if strcmp(LamType,'Sym')
    NbalVar   = 0;
    NMidPlane = 0;
    if Constraints.Vector(2)==0 % no 10% rule
        N10percentVar = 0; 
        if rem(MaxNplies,2) == 0
            NthetaVar = MaxNplies/2;
        else
            NthetaVar = floor(MaxNplies/2);
            NMidPlane = 1;
        end
    else
        N10percentVar = round((0.4*MaxNplies/2)/4)*4;
        if rem(MaxNplies,2) == 0
            NthetaVar = floor(MaxNplies/2-N10percentVar);
        else
            NthetaVar = floor(MaxNplies/2-N10percentVar);
            NMidPlane = 1;
        end
    end
    NtotalPly = 2*(NthetaVar + NbalVar + N10percentVar) + NMidPlane;
end



if strcmp(LamType,'Balanced')
    NMidPlane = 0;
    if Constraints.Vector(2)==0 % no 10% rule
        N10percentVar = 0; 
        if rem(MaxNplies,2) == 0
            NthetaVar = MaxNplies/2;
            NbalVar   = MaxNplies/2;
        else
            NthetaVar = floor(MaxNplies/2);
            NbalVar   = floor(MaxNplies/2);
            NMidPlane = 1;
        end
    else
        N10percentVar = round((0.4*MaxNplies)/4)*4;
        if rem(MaxNplies-N10percentVar,2) == 0
            NthetaVar = (MaxNplies-N10percentVar)/2;
            NbalVar   = (MaxNplies-N10percentVar)/2;
        else
            NthetaVar = floor((MaxNplies-N10percentVar)/2);
            NbalVar   = floor((MaxNplies-N10percentVar)/2);
            NMidPlane = 1;
        end
    end
    NtotalPly = NthetaVar + NbalVar + N10percentVar + NMidPlane;
end



if strcmp(LamType,'Balanced_Sym') % enforce even number of plies
    NMidPlane = 0;
    if Constraints.Vector(2)==0 % no 10% rule
        N10percentVar = 0; 
        if rem(MaxNplies,4) == 0
            NthetaVar = MaxNplies/4;
            NbalVar   = MaxNplies/4;
        else
            NthetaVar = floor(MaxNplies/4);
            NbalVar   = floor(MaxNplies/4);
            NMidPlane = MaxNplies-NthetaVar*2-NbalVar*2;
        end
    else

        N10percentVar = round((0.4*MaxNplies)/4)*4/2;
        NMissingPly = MaxNplies-N10percentVar;
        NthetaVar = 0;
        NbalVar   = 0;
        
        while NMissingPly>2
            NthetaVar = NthetaVar + 1;
            NbalVar   = NbalVar + 1;
            NtotalPly = 2*(NthetaVar + NbalVar + N10percentVar);
            NMissingPly = MaxNplies-NtotalPly;
        end  
        NMidPlane = NMissingPly;
    end
    NtotalPly = 2*(NthetaVar + NbalVar + N10percentVar)+ NMidPlane;
end

NDvs = NthetaVar + NbalVar + N10percentVar + NMidPlane;

if NtotalPly~=MaxNplies
    keyboard % should not happen
end




%% NdropVar

NdropPlies = NtotalPly - min(Nplies(:));

if strcmp(LamType,'Generic')
    NdropVar = NdropPlies;
end

if strcmp(LamType,'Sym') || strcmp(LamType,'Balanced')
    if rem(NdropPlies,2) == 0
        NdropVar = NdropPlies/2;
    else
        NdropVar = floor(NdropPlies/2) + 1;
    end
end

if strcmp(LamType,'Balanced_Sym')
    if rem(NdropPlies,4) == 0
        NdropVar = NdropPlies/4;
    else
        NdropVar = floor(NdropPlies/4) + 3;
    end
end

NdropVar = NdropVar + ceil(NMidPlane/2);

NStruct.NDvs            = NDvs;
NStruct.NthetaVar       = NthetaVar;
NStruct.NbalVar         = NbalVar;
NStruct.N10percentVar   = N10percentVar;
NStruct.NMidPlane       = NMidPlane;
NStruct.NtotalPly       = NtotalPly;
NStruct.NdropVar        = NdropVar;
