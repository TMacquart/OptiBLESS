function DetlaAngle = ComputeDeltaAngle(FiberAngles)

DetlaAngle = zeros(numel(FiberAngles)-1,1);

for iply = 1:numel(FiberAngles)-1
    if FiberAngles(iply)>=-45 && FiberAngles(iply)<=45
        if abs(FiberAngles(iply)==45) && abs(FiberAngles(iply+1))==90
            DetlaAngle(iply) = 45;
        else
            DetlaAngle(iply) = abs(FiberAngles(iply)-FiberAngles(iply+1));
        end
    elseif FiberAngles(iply)>45
        if FiberAngles(iply+1)>-45
            DetlaAngle(iply) = abs(FiberAngles(iply)-FiberAngles(iply+1));
        else
            DetlaAngle(iply) = abs(FiberAngles(iply)-(180+FiberAngles(iply+1)));
        end
    elseif FiberAngles(iply)<-45
        if FiberAngles(iply+1)<45
            DetlaAngle(iply) = abs(FiberAngles(iply)-FiberAngles(iply+1));
        else
            DetlaAngle(iply) = abs(FiberAngles(iply)-(-180+FiberAngles(iply+1)));
        end
    end
end