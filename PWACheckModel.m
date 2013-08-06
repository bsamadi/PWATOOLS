function [Err, model]=PWACheckModel(model);

% This function checks the dimensions of the system parameters to see if
% they are consistent with each other
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

Err=1;

model=DI_identifier(model);

if ~strcmp(model.pwatype, 'null')
    erA=CheckA_PWA(model);
    eraff=Checkaff_PWA(model);
    erD=CheckDomain_PWA(model);
    erBx=CheckB_PWA(model);
    erxcl=Checkxcl_PWA(model);
    erK=CheckK_PWA(model);
    Err=any([erA eraff erD erBx erxcl]); %if every thing is Ok Err==0, otherwise it remains 1
end