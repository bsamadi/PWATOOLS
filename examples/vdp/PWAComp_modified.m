function [pwainc, pwasys, nlsys]=PWAComp(model)
% This function gets the nonlinear model and calculates the PWA approximation
%
% input: nonlinear model: model
% output: PWA model     : pwasys
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

[nA, mA]=size(model.A);
nlsys.A      = model.A;
nlsys.B      = model.Bx;
nlsys.a      = model.aff;
nlsys.NonlinearFunction = model.fx;
nlsys.NonlinearStateEquations=model.NonlinearStateEquations;
nlsys.NonlinearDomain=model.NonlinearDomain;
nlsys.Domain = model.Domain;
nlsys.Method = model.mtd;
nlsys.UGR    = model.NR;
nlsys.TNR    = model.NR;
if isfield(model, 'Rstar') && ~isempty(model.Rstar)
    [nr, mr]=size(model.Rstar);
    nlsys.Rstar  = reshape(model.Rstar, max(nr, mr), min(nr, mr));
end

%% computing PWA Defferential Inclusions

nlsys.Resolution =5;
nlsys.ObjFun = 'L2';
nlsys.x0=rand(nA, 1);
nlsys.AbarLin = [model.A zeros(nA,1);zeros(1,nA+1)];   % Linear approximation of the nonlinear system
nlsys.xcl = model.xcl;

% PWA approximation
pwasys = pwa_approximate(nlsys);
pwasys.type='lower-envelope';

% PWA Envelope
nlsys.X = pwasys.X;
nlsys.Method = 'UpperBound';
pwasysU = pwa_approximate(nlsys);

nlsys.Method = 'LowerBound';
pwasysL = pwa_approximate(nlsys);

pwainc = pwasysU;

pwainc.Abar = [pwasysL.Abar  pwasysU.Abar];
pwainc.Bbar = [pwasysL.Bbar  pwasysU.Bbar];

pwainc.X = pwasysL.X;
pwainc.W = pwasysL.W;

pwainc.Z = {pwasysL.Z, pwasysU.Z};
pwainc.Y = {pwasysL.Y, pwasysU.Y};
pwainc.type='pwadi';

%% nlsys
clc
if isfield(nlsys, 'Rstar')
   nlsys=rmfield(nlsys, {'Rstar'});
end

nlsys=rmfield(nlsys, {'X', 'Resolution' , 'Method', 'UGR' , 'ObjFun' , 'x0',...
     'NonlinearDomain', 'Domain' , 'TNR',...
    'AbarLin'});
assignin('base', 'nlsys', nlsys);
%% writting the data of 'pwasys' in vector space R^n

pwa{1}=pwasys;
pwa{2}=pwainc;

for h=1:2
    pwasys=pwa{h};
    if (isfield(pwasys, 'L'))
        pwasys=rmfield(pwasys, 'L')
    end
    if (isfield(pwasys, 'l'))
        pwasys=rmfield(pwasys, 'l')
    end
    if (isfield(pwasys, 'Index'))
        pwasys=rmfield(pwasys, 'Index')
    end
    if (isfield(pwasys, 'X'))
        pwasys=rmfield(pwasys, 'X')
    end
    if (isfield(pwasys, 'Y'))
        pwasys=rmfield(pwasys, 'Y')
    end
    if (isfield(pwasys, 'T'))
        pwasys=rmfield(pwasys, 'T')
    end
    if (isfield(pwasys, 'W'))
        pwasys=rmfield(pwasys, 'W')
    end
    if (isfield(pwasys, 'Z'))
        pwasys=rmfield(pwasys, 'Z')    
    end
    if (isfield(pwasys, 'Err'))
        pwasys=rmfield(pwasys, 'Err')
    end
    pwatype=pwasys.type;
    if strcmp(pwatype, 'lower-envelope')
        col_index=[1];
    elseif strcmp(pwatype, 'pwadi')
        col_index=[1 2];
    elseif strcmp(pwatype, 'null')
        col_index=[];
    end
    
    Abar = pwasys.Abar;
    Bbar = pwasys.Bbar;
    Ebar = pwasys.Ebar;
    Fbar = pwasys.Fbar;
    xclbar = [pwasys.xcl;1];
    pwasys.istar = [];
    NR=size(Abar, 1);
    pwasys.NR=NR;
    
    for i=1:NR,
        pwasys.E{i}=pwasys.Ebar{i}(1:end-1, 1:end-1);
        pwasys.e{i}=pwasys.Ebar{i}(1:end-1, end);
        for j=col_index
            pwasys.A{i,j} = Abar{i,j}(1:end-1,1:end-1);
            pwasys.a{i,j} = Abar{i,j}(1:end-1,end);
            pwasys.B{i,j} = Bbar{i,j}(1:end-1,:);
        end
        xcl_is_inside_Ri = all(Ebar{i}*xclbar>=0-1e-7);
        if xcl_is_inside_Ri,
            pwasys.istar = union(pwasys.istar,i);                         % Center region(s)
        end
    end
    
    for i=1:NR
        for j=1:NR
            if ~isempty(Fbar{i,j})
                pwasys.F{i,j} = Fbar{i,j}(1:end-1,1:end-1);
                pwasys.f{i,j} = Fbar{i,j}(1:end-1,end);
            end
        end
    end
    pwa{h}=pwasys;
end
pwasys=pwa{1};
pwainc=pwa{2};
end


