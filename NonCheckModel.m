function [Err, model]=NonCheckModel(model);

% This function checks the dimensions of the system parameters to see if
% they are consistent with each other
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

Err=1;

erA=CheckA(model);
eraff=Checkaff(model);
erD=CheckDomain(model);
[erfx, model]=Checkfx(model);
erBx=CheckBx(model);
erxcl=Checkxcl(model);
Err=any([erA eraff erD erfx erBx erxcl]); %if every thing is Ok Err==0, otherwise it remains 1
