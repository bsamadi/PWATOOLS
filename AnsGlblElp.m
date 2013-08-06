function [pwasys,pwactrl]=AnsGlblElp(pwasys, option)
% This function analyzes the stability of the system by approximating the
% regions with ellipsoid. The regions can be slab or polytopic.
% Refrences:
% 1) L. Rodrigues and S. Boyd. Piecewise-affine state feedback for
%    piecewise-affine slab systems using convex optimization. Systems and
%    Control Letters, 54:835-853, 200

% 2) A. Hassibi and S. Boyd. Quadratic stabilization and control of
%    piecewise-linear systems. Proceedings of the American Control Conference,
%    6:3659-3664, 1998
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%
xcl = pwasys.xcl;
alpha=option.alpha;

[NR NS] = size(pwasys.Abar);    % Number of Systems, Number of Regions
n = size(pwasys.Abar{1},1)-1;   % Number of state variables
m = size(pwasys.Bbar{1},2);     % Number of inputs

istar = [];
A = pwasys.A;
a = pwasys.a;
B = pwasys.B;
K=pwasys.K;
k=pwasys.k;

pwatype=pwasys.type;
if strcmp(pwatype, 'lower-envelope')
    col_index=[1];
elseif strcmp(pwatype, 'upper-envelope')
    col_index=2;
elseif strcmp(pwatype, 'pwadi')
    col_index=[1 2];
elseif strcmp(pwatype, 'null')
    col_index=[];
end


%%  Regions equations with ellipsiodal approximations
% R=norm(E*z +e) < 1

if ~isfield(pwasys, 'EpA') | ~isfield(pwasys, 'Epb')
    pwasys=MinElp(pwasys);
end

for i=1:NR
    E{i}=pwasys.EpA{i};
    e{i}=pwasys.Epb{i};
end

%% extracting system data in seperated form from pwasys cell
for i=1:NR,
    xcl_is_inside_Ri = all(pwasys.E{i}*xcl+pwasys.e{i}>=0-1e-7);
    if xcl_is_inside_Ri,
        istar = union(istar,i);                         % Center region(s)
    end
end

%% shifting the equations with respect to the equilibrium point
% shifting a{i,j}
for i=1:NR,
    for j=col_index
        a{i,j} = a{i,j}+ A{i,j}*xcl;
    end
end
% shifting the ellipsidal regions equations. only e changes, E remains the
% same.
for i=1:NR
    e{i}=e{i}+ E{i}*xcl;
end


%% main part: closed loop

yalmip('clear');
constraints=set([]);

%% Definig variables in YALMIP
region_excluded_istar=setdiff([1:NR], istar);

P = sdpvar(n,n);

for i=region_excluded_istar
    miu{i}=sdpvar(1);
end

%% Central Region equations: DV

for i=istar
    for j=col_index
        Abarij=A{i,j}+B{i,j}*K{i};
        P_11 = P*Abarij+Abarij'*P+alpha*P;
        DV{i,j}=[P_11];
    end
end

%%  Regions exculded istar:  DV
for i=region_excluded_istar
    for j=col_index
        Abarij=A{i,j}+B{i,j}*K{i};
        abarij=a{i,j}+B{i,j}*k{i};
        P_11 = P*Abarij+Abarij'*P+alpha*P+miu{i}*E{i}'*E{i};
        P_12 = P*abarij+miu{i}*E{i}'*e{i};
        P_21 = P_12';
        P_22 = -miu{i}*(1-e{i}'*e{i});
        DV{i,j}=[P_11 P_12; P_21 P_22];
    end
end

%% constraints: Q, miu, DV for all regions

constraints=constraints+set(P>0);

for i=1:NR
    for j=col_index
        constraints=constraints+set(DV{i,j}<0);
    end
end

for i=region_excluded_istar
    constraints=constraints+set(miu{i}<0);
end

%% solution


u=solvesdp(constraints);
[u1, u2]=checkset(constraints);
pwactrl.problem=u.problem;
pwactrl.u1=u1;
pwactrl.u2=u2;
pwactrl.constraints=constraints;

pwactrl.P=double(P);
for i=region_excluded_istar
    pwactrl.miu{i}=double(miu{i});
end

end

