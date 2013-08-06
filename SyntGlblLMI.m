function [pwasys,pwactrl] = SyntGlblLMI(pwasys, option)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

pwactrl=[];

[NR NS] = size(pwasys.Abar);    % Number of Systems, Number of Regions
n = size(pwasys.Abar{1},1)-1;   % Number of state variables
m = size(pwasys.Bbar{1},2);     % Number of inputs

istar = [];
A = pwasys.A;
a = pwasys.a;
B = pwasys.B;
F = pwasys.F;
f = pwasys.f;
xcl=pwasys.xcl;
alpha=option.alpha;

pwatype=pwasys.type;
if strcmp(pwatype, 'lower-envelope')
    col_index=[1];
elseif strcmp(pwatype, 'pwadi')
    col_index=[1 2];
elseif strcmp(pwatype, 'null')
    col_index=[];
end

L=length(col_index);

%%  Regions equations with ellipsiodal approximations
% R=norm(E*z +e) < 1

if ~isfield(pwasys, 'EpA') | ~isfield(pwasys, 'Epb')
    pwasys=MinElp(pwasys);
    assignin('base', 'pwasys', pwasys);
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
    for j=1:NS,
        a{i,j} = a{i,j}+ A{i,j}*xcl;
    end
end
% shifting the intersection equations R_i \cap R_j= F_{ij}*s+f_{ij}
% F{i,j} remains the same. Only f{i,j} changes
for i=1:NR
    for j=1:NR
        if ~isempty(f{i,j})
            f{i,j} = f{i,j}-xcl;
        end
    end
end

% shifting the ellipsoidal regions equations. only e changes, E remains the
% same.
for i=1:NR
    e{i}=e{i}+ E{i}*xcl;
end

%% calculating a robust LQR gain for the R_istar region

K{istar(1)}=robust_LQR_piecewise_quad(istar(1), pwasys, option);

%% YALMIP initialzation
yalmip('clear');
constraints=set([]);

%% Definig variables in YALMIP
% Y_i= K_i * Q
Q = sdpvar(n,n);
etta=sdpvar(1);
region_excluded_istar_1=setdiff([1:NR], istar(1));
region_excluded_istar=setdiff([1:NR], istar);

for i=region_excluded_istar_1
    Y{i} = sdpvar(m,n);
end
Y{istar(1)} = K{istar(1)}*Q;

for i=region_excluded_istar
    miu{i}=sdpvar(1);
    Z{i}= sdpvar(m,1);
end

for i=istar
    miu{i}=-1;
    Z{i}=zeros(m,1);
end

%% Central Region equations: Y, Z and DV
% inequality for negative definitness vdot < alpha * v
for i=istar
    for j=col_index
        P_11 = A{i,j}*Q+Q*A{i,j}'+B{i,j}*Y{i}+Y{i}'*B{i,j}';
        DV{i, j}=[P_11];
    end
end

%%  Regions exculded istar:  DV

Ed=size(pwasys.E{1}, 1); % shows if the regions are slab or 

if  Ed==2  % slab regions
    for i=region_excluded_istar
        for j=col_index
            P_11 = A{i,j}*Q+Q*A{i,j}'+B{i,j}*Y{i}+Y{i}'*B{i,j}'+alpha*Q+...
                miu{i}*a{i,j}*a{i,j}'+a{i,j}*Z{i}'*B{i,j}'+...
                B{i,j}*Z{i}*a{i,j}';
            P_12 = (miu{i}*a{i,j}+B{i,j}*Z{i})*e{i}'+Q*E{i}';
            P_21 = P_12';
            P_22 = -miu{i}*(1-e{i}*e{i}');
            DV{i,j}=[P_11 P_12; P_21 P_22];
        end
    end
    
elseif Ed>2
    % polytopic regions
    % Copy Right @ Mohsen Zamani Fekri, 2010
    for i=region_excluded_istar
        for j=col_index
            P_11 = A{i,j}*Q+Q*A{i,j}'+B{i,j}*Y{i}+Y{i}'*B{i,j}'+alpha*Q;
            P_12 = (miu{i}*a{i,j}+B{i,j}*Z{i})+Q*E{i}'*e{i};
            P_21 = P_12';
            P_22 = -miu{i}*(1-e{i}'*e{i});
            DV{i, j}=[P_11 P_12; P_21 P_22];
            
        end
    end
end
%% constraints: Q, miu, DV for all regions

constraints=constraints+set(Q>0,['Q' '>0']);

for i=1:NR
    for j=col_index
        constraints=constraints+set(DV{i,j}<0,['DV' num2str(i) '-' num2str(j) '<0']);
    end
end

for i=region_excluded_istar
    constraints=constraints+set(miu{i}<0, ['miu' num2str(i) '<0']);
end

%% solution
u=solvesdp(constraints);
[u1 u2]=checkset(constraints);
pwactrl.problem=u.problem;
pwactrl.u1=u1;
pwactrl.u2=u2;
pwactrl.constraints=constraints;

Q_double=double(Q);
pwactrl.Q=Q_double;
pwactrl.P=inv(Q_double);

for i=1:region_excluded_istar
    miu_double=double(miu{i});
    pwactrl.miu{i}=miu_double;
end

% Z{i} is zero for i=istar
for i=1:NR
    miu_double=double(miu{i});
    Y_double=double(Y{i});
    Z_double=double(Z{i});
    pwactrl.Z{i}=Z_double;
    pwactrl.Y{i}=Y_double;
    pwactrl.Kbar{i} = [Y_double*inv(Q_double) Z_double/miu_double];
end
pwactrl.Ebar = pwasys.Ebar;
pwactrl.xcl = xcl;

