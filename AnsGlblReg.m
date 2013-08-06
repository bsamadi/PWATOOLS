function [pwasys,pwactrl]=AnsGlblReg(pwasys, option)
% This function analyses the stability of the system by approximating the
% regions with polytopic expression E'*Z*E for a matrix Z of proper size
% whose elements are non-negative.
%
% References
% 1) B. Samadi and L. Rodrigues. Extension of local linear controllers to
%    global piecewise affine controllers for uncertain non-linear systems.
%    International Journal of Systems Science, 39(9):867-879, 2008.
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
E = pwasys.E;
e = pwasys.e;
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



%% extracting system data in seperated form from pwasys cell
for i=1:NR,
    xcl_is_inside_Ri = all(E{i}*xcl+e{i}>=0-1e-7);
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



%% Definig variables in YALMIP
region_excluded_istar=setdiff([1:NR], istar);
yalmip('clear');
constraints=set([]);

P = sdpvar(n,n);

dE=size(E{1}, 1);
for i=region_excluded_istar
    Z{i}=sdpvar(dE);
end

%% Central Region equations:  DV

for i=istar
    for j=col_index
        Abarij=A{i,j}+B{i,j}*K{i};
        P_11 = P*Abarij+Abarij'*P+alpha*P;
        DV{i,j}=[P_11];
    end
end

%%  Regions excluded istar:  DV
for i=region_excluded_istar
    for j=col_index
        Abarij=A{i,j}+B{i,j}*K{i};
        abarij=a{i,j}+B{i,j}*k{i};
        P_11 = P*Abarij+Abarij'*P+alpha*P+E{i}'*Z{i}*E{i};
        P_12 = P*abarij+E{i}'*Z{i}*e{i};
        P_21 = P_12';
        P_22 = e{i}'*Z{i}*e{i};
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
    constraints=constraints+set(Z{i}(:)>0);
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
    pwactrl.Z{i}=double(Z{i});
end
