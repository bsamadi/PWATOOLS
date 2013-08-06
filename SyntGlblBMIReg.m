function [pwasys,pwactrl] = SyntGlblBMIReg(pwasys, option)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%


[NR NS] = size(pwasys.Abar);    % Number of Systems, Number of Regions
n = size(pwasys.Abar{1},1)-1;   % Number of state variables
m = size(pwasys.Bbar{1},2);     % Number of inputs

alpha=option.alpha;

istar = [];
A = pwasys.A;
a = pwasys.a;
B = pwasys.B;
E=pwasys.E;
e=pwasys.e;
F = pwasys.F;
f = pwasys.f;
xcl=pwasys.xcl;

pwatype=pwasys.type;
if strcmp(pwatype, 'lower-envelope')
    col_index=[1];
elseif strcmp(pwatype, 'pwadi')
    col_index=[1 2];
elseif strcmp(pwatype, 'null')
    col_index=[];
end
L=length(col_index);

%% extracting system data in seperated form from pwasys cell
for i=1:NR,
    xcl_is_inside_Ri = all(pwasys.E{i}*xcl+pwasys.e{i}>=0-1e-7);
    if xcl_is_inside_Ri,
        istar = union(istar,i);                         % Center region(s)
    end
end

region_excuded_istar=setdiff([1:NR], istar);


%% shifting the equations with respect to the equilibrium point
% shifting a{i,j}
for i=1:NR,
    for j=col_index,
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

% shifting the ellipsidal regions equations. only e changes, E remains the
% same.
for i=1:NR
    e{i}=e{i}+ E{i}*xcl;
    sys_E{i}=pwasys.E{i};
    sys_e{i}=pwasys.E{i}*xcl+pwasys.e{i};
end


%% building a shifted "pwasys" model.

for i=1:NR
    shifted_model.Ebar{i}=pwasys.Ebar{i};
    shifted_model.Ebar{i}(1:end-1, 1:end)=[sys_E{i}  sys_e{i}];
    for j=col_index
        shifted_model.Abar{i,j}=pwasys.Abar{i,j};
        shifted_model.Abar{i,j}(1:n, n+1)=a{i,j};
        
    end
end
for i=1:NR
    for j=1:NR
        shifted_model.Fbar{i,j}=[];
        if ~isempty(f{i,j})
            shifted_model.Fbar{i,j}=pwasys.Fbar{i,j};
            shifted_model.Fbar{i,j}(1:n, n)=f{i,j};
        end
    end
end

shifted_model.Bbar=pwasys.Bbar;

shifted_model.type=pwatype;

%% YALMIP initialzation
yalmip('clear');
constraints=set([]);

%% DEFINIG VARIABLES IN YALMIP DEFINING K

clc
disp('Defining gain matrix K_i and k_i')
% defining varibales [K_i k_i] for a piecewise-quadratic lyapunov
% function
NR=size(pwasys.Abar,1);
Num_S=size(pwasys.Abar{1,1},1)-1; % number of states
Num_I=size(pwasys.Bbar{1,1},2);   % number of inputs


xy_list=[];
for x=1:NR
    for y=x+1:NR
        if ~isempty(pwasys.Fbar{x,y})
            xy_list=[xy_list; x y];
        end
    end
end
Num_rltn=size(xy_list,1); % number of unique F_ij relation
rgn_lngth=(Num_S+1)*Num_I;

istar=[];
for i_I=1:NR,
    xcl_is_inside_Ri = all(pwasys.Ebar{i_I}*[xcl;1]>=0-1e-7);
    if xcl_is_inside_Ri,
        istar = union(istar,i_I);                         % Center region(s)
    end
end

% converting MXN=0 equations to a unique AY=0 equation
%%-------------------------------------------------------------------------
%%all directons (more conservative)
if option.NormalDirectionOnly
    % to force the normal direction and bypass the whole direction approach
    row_check=0;  
else
    [reshaping, input_info]=CtrlContAllDir(istar, shifted_model);
    affne_indx=NR*rgn_lngth+1: NR*rgn_lngth+NR;
    indp_affne_indx=intersect(affne_indx, input_info.d);
    row_check=1; % to make it a DEFINED variable
    for i_I=indp_affne_indx
        [r_I,c_I]=find(input_info.d==i_I);
        row_check=input_info.w(c_I, :);
        row_check(i_I)=[];
    end
    manual_istar_indx=(NR-1)*rgn_lngth+1:(NR-1)*rgn_lngth+Num_S*Num_I; %indices of varibles
    % we intend to set manually as the equilibrium region
    indp_manual_indx=intersect(manual_istar_indx, input_info.d);
end

if norm(row_check)==0 || ~isempty(indp_manual_indx)
    %normal direction (less conservative)
    [reshaping, input_info]=CtrlContNormalDir(istar, shifted_model);
    affne_indx=NR*rgn_lngth+1: NR*rgn_lngth+NR;
    indp_affne_indx=intersect(affne_indx, input_info.d);
    for i_I=indp_affne_indx
        [r_I,c_I]=find(input_info.d==i_I);
        row_check=input_info.w(c_I, :);
        row_check(i_I)=[];
        if norm(row_check)==0 %second time we check
            disp('error! equilibrium constraints cannot be satisfied')
            return
        end
    end
    manual_istar_indx=(NR-1)*rgn_lngth+1:(NR-1)*rgn_lngth+Num_S*Num_I; %indices of varibles
    % we intend to set manually as the equilibrium region
    indp_manual_indx=intersect(manual_istar_indx, input_info.d);
    if ~isempty(indp_manual_indx) %second time we check
        disp('linear relation exists between istar region gains')
        disp('warning! manual gain selecting for istar region is not allowed')
        return
    end
end


clear K_mtrx K_vec k_var k_temp

k_var=sdpvar(1, NR*rgn_lngth);

K_robust=robust_LQR_piecewise_quad(istar(1),pwasys,option);
k_var(manual_istar_indx)=reshape(K_robust, 1, Num_S*Num_I);

k_temp={};

tot_index_I=1:size(input_info.S, 2); % indices of total variables
indp_index_I=input_info.d;
depn_index_I=setdiff(tot_index_I, indp_index_I); % dependent variables


for m_I=indp_index_I
    k_var(m_I)=0;
end

for m_I=depn_index_I
    % dependent vector corresponds to an independent
    % variable. it is true even if the vector is completely zero.
    if  ismember(m_I, affne_indx)
        k_var(m_I)=1;
    elseif (~ismember(m_I, manual_istar_indx) & ~ismember(m_I, affne_indx)) % manual istar gains should not be set again
        %elseif (~ismember(m_I, affne_indx))
        k_var(m_I)=sdpvar(1);
        k_temp{size(k_temp,2)+1}=k_var(m_I);
    end
    if  norm(input_info.S(:,m_I))~=0 % if the vector is not completey zero,
        % then the vsribles corresponding to
        % its bases will depend on the varibale
        % corresponding to lyap_info.S(:,m_I)
        
        [base_index_I, non_impt, base_coef]=find(input_info.w(:,m_I));
        for k_I=1:length(base_index_I)
            b_indx_indp=base_index_I(k_I); % index_P of the independent vector
            b_indx_I=input_info.d(b_indx_indp);
            val_indx_I=base_coef(k_I); % coefficient by which the independent
            k_var(b_indx_I)=k_var(b_indx_I)-val_indx_I*(k_var(m_I)); % li
        end
    end
end



for i_I=1:NR
    m_I=reshaping(i_I);
    K_vec{m_I}=k_var((i_I-1)*rgn_lngth+1:i_I*rgn_lngth);
end


for i_I=1:NR
    K_mtrx{i_I}=sdpvar(Num_I, Num_S+1);
    K_mtrx{i_I}=reshape(K_vec{i_I}, Num_I, Num_S+1);
end
%% end of defining variables


Q=sdpvar(n);

for i=1:NR
    K{i}=K_mtrx{i}(:,1:Num_S);
    k{i}=K_mtrx{i}(:,Num_S+1);
end

dE=size(E{1}, 1);
for i=region_excuded_istar
    Z{i}=sdpvar(dE);
end


%% Central Region equations: Y, Z and DV

%inequality for negative definitness vdot < alpha * v
for i=istar
    for j=col_index
        Abarij=A{i,j}+B{i,j}*K{i};
        P_11 = Q*Abarij+Abarij'*Q+alpha*Q;
        DV{i,j}=[P_11];
    end
end

%%  Regions exculded istar:  DV


    
    %  Regions exculded istar:  DV
    for i=region_excuded_istar
        for j=col_index
            Abarij=A{i,j}+B{i,j}*K{i};
            abarij=a{i,j}+B{i,j}*k{i};
            P_11 = Q*Abarij+Abarij'*Q+alpha*Q+E{i}'*Z{i}*E{i};
            P_12 = Q*abarij+E{i}'*Z{i}*e{i};
            P_21 = P_12';
            P_22 = e{i}'*Z{i}*e{i};
            DV{i,j}=[P_11 P_12; P_21 P_22];
        end
    end
    

constraints=set([]);
constraints=constraints+set(Q>0,['Q' '>0']);

for i=1:NR
    for j=col_index
        constraints=constraints+set(DV{i,j}<0, ['DV' num2str(i) '-' num2str(j) '<0']);
    end
end

for i=region_excuded_istar
    constraints=constraints+set(Z{i}(:)>0,['Z' num2str(i) '>0']);
end

%
%% solution
u=solvesdp(constraints);
[u1 u2]=checkset(constraints);
pwactrl.problem=u.problem;
pwactrl.u1=u1;
pwactrl.u2=u2;
pwactrl.constraints=constraints;

pwactrl.xcl = xcl;
% pwactrl.istar = istar;

Q=double(Q);
pwactrl.Q = Q;
for i=1:NR
    K{i}=double(K{i});
    k{i}=double(k{i});
end

for i=region_excuded_istar
    Z{i}=double(Z{i});
    pwactrl.Z{i} = Z{i};
end

for i=1:NR,
    pwactrl.Kbar{i} = [K{i} k{i}];
end

pwactrl.Ebar = pwasys.Ebar;
