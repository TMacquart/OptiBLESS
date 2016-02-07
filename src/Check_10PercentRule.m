function FEASIBLE = Check_10PercentRule(GuideAngles)

FEASIBLE = true;

TenpercentDV = round(length(GuideAngles)*0.1);
N0Plies      = length(find(GuideAngles==0));
N90Plies     = length(find(abs(GuideAngles)==90));
N45Plies     = length(find(GuideAngles==45));
NM45Plies    = length(find(GuideAngles==-45));

if N0Plies<TenpercentDV || N90Plies<TenpercentDV || N45Plies<TenpercentDV || NM45Plies<TenpercentDV
    FEASIBLE = false;
end


end