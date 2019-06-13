function rmsValue=MYrms(x)
%% This function use to calculating rms with out using matlab rms functions
N=length(x);
rmsValue=sqrt(1/N.*(sum(x(:).^2)));
end 