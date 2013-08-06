function [reshaping, M]=CtrlContAllDir(istar, pwasys)
% expression AXB with [m,n]=size(A); [n,k]=size(X) and [k,p]=size(B) is
% given. X is variable and A and B are known. we rewrite AXB as S*X_vec
% where X_vec= reshape(X, nk,1) and [mp, nk]=size(S).
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

pwatype=pwasys.type;
if strcmp(pwatype, 'lower-envelope')
    col_index=[1];
elseif strcmp(pwatype, 'pwadi')
    col_index=[1 2];
elseif strcmp(pwatype, 'null')
    col_index=[];
end
L=length(col_index);  


Fbar=pwasys.Fbar;
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
rgn_lngth=[(Num_S+1)*Num_I+1];
M={};

S=zeros(L*Num_rltn*Num_S^2, NR*rgn_lngth);

rltn_index=0;
for x=1:NR
    for y=x+1:NR
        if ~isempty(Fbar{x,y})
            rltn_index=rltn_index+1;
            S_temp_x=K_row(x, y,istar, pwasys);
            S_temp_y=K_row(y, x,istar, pwasys);
            row_vec=(rltn_index-1)*(L*Num_S^2)+1:rltn_index*(L*Num_S^2);
            col_vec_x=(x-1)*rgn_lngth+1:x*rgn_lngth;
            col_vec_y=(y-1)*rgn_lngth+1:y*rgn_lngth;
            S(row_vec, col_vec_x)=S_temp_x;
            S(row_vec, col_vec_y)=-S_temp_y;
        end
    end
end


S_eql=[];
for x=1:NR
    if ismember(x, istar)
        temp=zeros(L*Num_S, NR*rgn_lngth);
        Bkij=EqlCnst(x,pwasys);
        temp(:, (x-1)*rgn_lngth+1:x*rgn_lngth)=Bkij;
        S_eql=[S_eql;temp];
    end
end

S=[S; S_eql];

affne_indx=[rgn_lngth:rgn_lngth:NR*rgn_lngth];
S_affne=S(:,affne_indx);
S(:,affne_indx)=[];

% change of the region length
rgn_lngth=rgn_lngth-1;

tot_istar_indx=[];
for i=istar
    istar_indx=[(i-1)*rgn_lngth+1: i*rgn_lngth];
    tot_istar_indx=[istar_indx tot_istar_indx];
end
S_istar=S(:,tot_istar_indx);
S(:,tot_istar_indx)=[];
S=[S S_istar S_affne];
reshaping=([setdiff(1:NR, istar) fliplr(istar)]); %order of regions in modified S
 S=round(S*10^8)/10^8;

[w,d]=rref(S);
M.S=S;
M.w=w;
M.d=d;
end


%% first local function
function S_local=local_multiplier_builder(A, X, B)
[m,n]=size(A);
[n,k]=size(X);
[k,p]=size(B);
alpha=[];

% synthesizing cooeffcients of matrix S
for i=1:m
    for j=1:p
        for r=1:n
            for q=1:k
                alpha(i,j,r,q)=A(i,r)*B(q,j);
            end
        end
    end
end


S_local=[];
%synthesizing matrix S from cooeffcients alpha
for j=1:p
    for i=1:m
        v=reshape(alpha(i,j, :,:), 1, n*k);
        S_local=[S_local;v];
    end
end

end

%% second local function
function  S_temp=K_row(x,y,istar, pwasys)

pwatype=pwasys.type;
if strcmp(pwatype, 'lower-envelope')
    col_index=[1];
elseif strcmp(pwatype, 'pwadi')
    col_index=[1 2];
elseif strcmp(pwatype, 'null')
    col_index=[];
end
L=length(col_index);  

NR=size(pwasys.Abar,1);
Num_S=size(pwasys.Abar{1,1},1)-1; % number of states
Num_I=size(pwasys.Bbar{1,1},2);   % number of inputs
Fij=pwasys.Fbar{x,y}(1:Num_S, 1:Num_S-1);
fij=pwasys.Fbar{x,y}(1:Num_S, Num_S);
cij=null(Fij');

S_temp=[];
for envelope_index=col_index
    % synthesizing equivlent of BKF_{ij}
    A=pwasys.Bbar{x,envelope_index}(1:Num_S,:);
    X=rand(Num_I,Num_S);
    B=Fij;
    BKF_ij=local_multiplier_builder(A,X, B);
    
    % synthesizing equivlent of BKf_{ij}
    A=pwasys.Bbar{x,envelope_index}(1:Num_S,:);
    X=rand(Num_I,Num_S);
    B=fij;
    BKf_ij=local_multiplier_builder(A,X, B);
    
    %synthesizing Bk_ij
    A=pwasys.Bbar{x,envelope_index}(1:Num_S,:);
    X=rand(Num_I,1);
    B=eye(1);
    Bk_ij=local_multiplier_builder(A,X, B);
     
    %synthesizing S
    temp=[BKF_ij       zeros(size(BKF_ij,1), Num_I) zeros(size(BKF_ij,1),1);
          BKf_ij       Bk_ij                        zeros(size(BKf_ij,1), 1)];
    
    S_temp=[S_temp; temp];
end

end


%% Third local function: equilibrium constraints

function  Bkij=EqlCnst(x,pwasys)
pwatype=pwasys.type;
if strcmp(pwatype, 'lower-envelope')
    col_index=[1];
elseif strcmp(pwatype, 'pwadi')
    col_index=[1 2];
elseif strcmp(pwatype, 'null')
    col_index=[];
end
L=length(col_index);  

NR=size(pwasys.Abar,1);
Num_S=size(pwasys.Abar{1,1},1)-1; % number of states
Num_I=size(pwasys.Bbar{1,1},2);   % number of inputs

%synthesizing Bk_ij
Bkij=[];
for envelope_index=col_index
    A=pwasys.Bbar{x,envelope_index}(1:Num_S,:);
    X=rand(Num_I,1);
    B=eye(1);
    Bk_ij=local_multiplier_builder(A,X, B);
    aij=pwasys.Abar{x, envelope_index}(1:Num_S, Num_S+1);
    Bkij=[Bkij; zeros(Num_S, Num_S*Num_I)  Bk_ij   aij];
end
end