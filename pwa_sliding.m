function pwactrl = pwa_lin2pwa(pwasys, Param)

% pwasys
% xbarDot = Abar{i}*xbar + Bbar{i}*u
%
% pwasys.Abar = {Abar1, ..., Abarn}
% pwasys.Bbar = {Bbar1, ..., Bbarn}
%
% Polytopic region i: {x| Ebar{i}*x > 0 }
%
%
% Boundary information table:
% Boundary between Ri and Rj : {x| xbar=Fbar{i,j}*sbar }
%

xcl = pwasys.xcl;
xclbar = [xcl;1];

NR = length(pwasys.Abar);       % Number of regions
n = size(pwasys.Abar{1},1)-1;   % Number of state variables
m = size(pwasys.Bbar{1},2);     % Number of inputs

Abar = pwasys.Abar;
Bbar = pwasys.Bbar;
Ebar = pwasys.Ebar;
Fbar = pwasys.Fbar;

istar = [];
for i=1:NR,
    A{i} = Abar{i}(1:end-1,1:end-1);
    a{i} = Abar{i}(1:end-1,end);
    B{i} = Bbar{i}(1:end-1,:);
    Echeck{i} = [zeros(1,n) 1; Ebar{i}];
    xcl_is_inside_Ri = all(Ebar{i}*xclbar>=0);
    if xcl_is_inside_Ri,
        istar = union(istar,i);
    end
end



Lin2PWA = 1;
if exist('Poles'),
    for i = istar,
        K{i} = -acker(A{i},B{i},Poles);
    end
elseif exist('KLin'),
    for i = istar,
        K{i} = KLin;
    end
elseif exist('QLin'),
    for i = istar,
        K{i} = -lqr(A{i},B{i},QLin,RLin);
    end
else
    Lin2PWA = 0;
end

if isfield(Param,'Sbar'),
    Sbar = Param.Sbar;
else,
    error('Sliding mode is not defined.');
end

if isfield(Param,'Sbar'),
    Sbar = Param.Sbar;
else,
    error('Sliding mode is not defined.');
end

if isfield(Param,'alpha'),
    alpha=Param.alpha;
else
    alpha = 5;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Sliding Mode Control %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:NR,
    Kbar{i} = -(Sbar*Bbar{i})\(Sbar*Abar{i}+alpha*Sbar);
end

pwactrl.xcl = xcl;
pwactrl.istar = istar;

for i=1:NR,
    pwactrl.Kbar{i} = double(Kbar{i});
end

pwactrl.Index = pwasys.Index;
pwactrl.X = pwasys.X;
pwactrl.T = pwasys.T;
pwactrl.xcl = pwasys.xcl;
pwactrl.Ebar = pwasys.Ebar;