function L2Gain = pwa_l2gain(pwasys, pwactrl)

% pwasys
% xbarDot = Abar{i}*xbar + Bbar{i}*u
%
% pwasys.Abar = {Abar1, ..., Abarn}
% pwasys.Bbar = {Bbar1, ..., Bbarn}
% pwasys.Cbar = Cbar
% pwasys.D = D
%
% Polytopic region i: {x| Ebar{i}*x > 0 }
%
%
% Boundary information table:
% Boundary between Ri and Rj : {x| xbar=Fbar{i,j}*sbar }
%

NR = length(pwasys.Abar);       % Number of regions
n = size(pwasys.Abar{1},1)-1;   % Number of state variables
m = size(pwasys.Bbar{1},2);     % Number of control inputs
q = size(pwasys.Bwbar{1},2);    % Number of disturbance inputs

Lyapunov = 'Global';

xcl = pwasys.xcl;
xclbar = [xcl;1];

if exist('pwactrl'),
    for i = 1:NR,
        Abar{i} = pwasys.Abar{i} + pwasys.Bbar{i}*pwactrl.Kbar{i};
    end
else
    Abar = pwasys.Abar;
end

Bwbar = pwasys.Bwbar;
Cbar = pwasys.Cbar;
Dw = pwasys.Dw;
Ebar = pwasys.Ebar;
Fbar = pwasys.Fbar;

istar = [];
for i=1:NR,
    A{i} = Abar{i}(1:end-1,1:end-1);
    Bw{i} = Bwbar{i}(1:end-1,:);
    C{i} = Cbar{i}(:,1:end-1);
    Echeck{i} = [zeros(1,n) 1; Ebar{i}];
    xcl_is_inside_Ri = all(Ebar{i}*xclbar>=0);
    if xcl_is_inside_Ri,
        istar = union(istar,i);
    end
end

istar = istar(end);

if Cbar{istar}*xclbar~=0,
    error('Cbar*xcbar is not zero!');
end

I = eye(n);
One = zeros(n+1);
One(end,end)=1;

yalmip('clear');
constraints=set([]);
p = size(Ebar{1},1);

gamma = sdpvar(1);

% Lyapunov function
if strcmp(Lyapunov,'Global'),
    Pg = sdpvar(n,n);
    for i=1:NR,
        P{i} = Pg;
        r{i} = 0;
    end
else,
    for i=1:NR,
        if i==istar,
            P{i} = sdpvar(n,n);
            r{i} = 0;
        else
            P{i} = sdpvar(n,n);
            r{i} = sdpvar(1,1);
        end
    end
end

for i=1:NR,
    % Constrainted Pbar
    Pbar{i} = [P{i} -P{i}*xcl; -xcl'*P{i} xcl'*P{i}*xcl+r{i}];
end

constraints=constraints+set(P{istar} > 0,['P' num2str(istar) '>0']);
L2Mat{istar} = [P{istar}*A{istar}+A{istar}'*P{istar}+C{istar}'*C{istar} P{istar}*Bw{istar}+C{istar}'*Dw{istar}; (P{istar}*Bw{istar}+C{istar}'*Dw{istar})' Dw{istar}'*Dw{istar}-gamma*eye(q)];
constraints=constraints+set( L2Mat{istar} < 0,['L2Mat' num2str(istar) '<0']);

for i=setdiff(1:NR,istar),
    %%%%%%%%%%%%%%%%%%%%%%%% Inequality Constraints %%%%%%%%%%%%%%%%%%%
    L{i} = sdpvar(p+1,p+1);
    constraints=constraints+set(0 < L{i}(:) < 1e5,['L{' num2str(i) '}(:)>0']);
    L2Mat{i} = [Pbar{i}*Abar{i}+Abar{i}'*Pbar{i}+Cbar{i}'*Cbar{i}+Echeck{i}'*L{i}*Echeck{i} Pbar{i}*Bwbar{i}+Cbar{i}'*Dw{i}; (Pbar{i}*Bwbar{i}+Cbar{i}'*Dw{i})' Dw{i}'*Dw{i}-gamma*eye(q)];
    constraints=constraints+set( L2Mat{i} < 0,['L2Mat' num2str(i) '<0']);
end

solution=solvesdp(constraints,gamma,sdpsettings('usex0',1,'shift',1e-7));
checkset(constraints);
gamma = double(gamma);
L2Gain = sqrt(gamma);