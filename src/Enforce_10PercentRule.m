% =====                                                              ==== 
% Enforcing_10PercentRule converts the ply angles close to [+-45 90 0] such as 
%   at least 10% of each +45/-45/0/90 plies are present in the laminate.
%   Note that this may not always be possible.
%
% [GuideAngles] = Enforcing_10PercentRule(GuideAngles)
% =====                                                              ==== 



% ----------------------------------------------------------------------- %
% Copyright (c) <2015>, <Terence Macquart>
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% The views and conclusions contained in the software and documentation are those
% of the authors and should not be interpreted as representing official policies,
% either expressed or implied, of the FreeBSD Project.
% ----------------------------------------------------------------------- %




function [GuideAngles] = Enforce_10PercentRule(GuideAngles,unit)

TenpercentDV = round(length(GuideAngles)*0.1);
N0Plies      = length(find(GuideAngles==0));
N90Plies     = length(find(abs(GuideAngles)==90));
N45Plies     = length(find(GuideAngles==45));
NM45Plies    = length(find(GuideAngles==-45));
indexPly     = 1:length(GuideAngles);


if N0Plies < TenpercentDV
    PlyTable  = sortrows([abs(GuideAngles)' indexPly'],1);
    PlyIndex0 = PlyTable(1:TenpercentDV,2);
    GuideAngles(PlyIndex0) = 0;
end


if N45Plies < TenpercentDV
    PlyTable   = sortrows([abs(GuideAngles-45)' indexPly'],1);
    PlyIndex45 = PlyTable(1:TenpercentDV,2);
    GuideAngles(PlyIndex45) = 45;
end


if NM45Plies < TenpercentDV
    PlyTable   = sortrows([abs(GuideAngles+45)' indexPly'],1);
    PlyIndex45 = PlyTable(1:TenpercentDV,2);
    GuideAngles(PlyIndex45) = -45;
end


if N90Plies < TenpercentDV
    PlyTable   = sortrows([abs(GuideAngles)' indexPly'],1);
    PlyIndex90 = PlyTable(end-TenpercentDV+1:end,2);
    A = [-1 1];
    GuideAngles(PlyIndex90) = 90*A(ceil(rand*2));
end

end