function [reshaping, M]=LyapContGeneral(pwasys)
% expression AXB with [m,n]=size(A); [n,k]=size(X) and [k,p]=size(B) is
% given. X is variable and A and B are known. we rewrite AXB as S*X_vec
% where X_vec= reshape(X, nk,1) and [mp, nk]=size(S).

% this function considers V(x)=x'px + 2q'x+r function as the Lyapunov
% candiadate
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

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

unrptd_var_length=0;
for x=1:Num_S+1  %% updated from original file for V(x)=x'px
    unrptd_var_length=unrptd_var_length+x;
end


M={};
S=zeros(Num_rltn*Num_S^2, NR*unrptd_var_length);
rltn_index=0;
for x=1:NR
    for y=x+1:NR
        if ~isempty(Fbar{x,y})
            rltn_index=rltn_index+1;
            % synthesizing equivalent of F_{ij}KF_{ij}
            A=Fbar{x,y}';
            X=rand(Num_S+1);
            B=Fbar{x,y};
            [reshaping, F_ijX_iF_ij]=local_multiplier_builder(A,X, B);
            [reshaping, F_ijX_jF_ij]=local_multiplier_builder(A,X, B);
            
            
            % synthesizing S
            row_vec=(rltn_index-1)*Num_S^2+1:rltn_index*Num_S^2;
            col_vec_x=(x-1)*unrptd_var_length+1:x*unrptd_var_length;
            col_vec_y=(y-1)*unrptd_var_length+1:y*unrptd_var_length;
            S(row_vec, col_vec_x)=F_ijX_iF_ij;
            S(row_vec, col_vec_y)=-F_ijX_jF_ij;
            
            
        end
    end
end


%% for i=istar set q=0

xcl=zeros(Num_S, 1);
istar=[];
for i_P=1:NR,
    xcl_is_inside_Ri = all(pwasys.Ebar{i_P}*[xcl;1]>=0-1e-7);
    if xcl_is_inside_Ri,
        istar = union(istar,i_P);                         % Center region(s)
    end
end

istar

[m_S, n_S]=size(S);
for x=istar
    A=[eye(Num_S) zeros(Num_S, 1)];
    X=rand(Num_S+1);
    B=[zeros(Num_S, 1); 1];
    
    [reshaping, F_ijX_iF_ij]=local_multiplier_builder(A,X, B);
    
    % synthesizing S
    row_vec=m_S+1:m_S+Num_S;
    col_vec=(x-1)*unrptd_var_length+1:x*unrptd_var_length;
    S(row_vec, col_vec)=F_ijX_iF_ij;
    m_S=m_S+Num_S;
end

%% analyzing matrix S
S=round(S*10^8)/10^8;
[w,d]=rref(S);
M.S=S;
M.w=w;
M.d=d;
% w is a form of dependency between vectors and d is indices of
% independent vectors
% archving w and d for every calculation


end

%% local function
function [reshaping,S_local]=local_multiplier_builder(A, X, B)
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
        reshaping=[];
        index=0;
        for q=1:n
            for r=q:n
                index=index+1;
                if (r==q)
                    v(index)=alpha(i,j,r,q);
                else
                    v(index)=alpha(i,j,r,q)+alpha(i,j,q,r);
                end
                reshaping=[reshaping; [index r q]];
            end
        end
        S_local=[S_local;v];
    end
end
end