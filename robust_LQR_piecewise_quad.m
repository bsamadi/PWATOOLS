function [K]=robust_LQR_piecewise_quad(istar, pwasys, Param)

%%%% robust LQR computation for R_istar region, based on Jadbababie Theorem
%%%% Theorem 5.1 in [Samadi and Rodrigues, 2008]
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

Q = Param.QLin;
R = Param.RLin;


if L==1
    K=lqr(pwasys.A{istar}, -pwasys.B{istar}, Q, R);
else
    xcl = pwasys.xcl;
    xclbar = [xcl;1];
    
    [NR NS] = size(pwasys.Abar);    % Number of Systems, Number of Regions
    n = size(pwasys.Abar{1},1)-1;   % Number of state variables
    m = size(pwasys.Bbar{1},2);     % Number of inputs
    
    % if isfield(Param,'alpha'),      % VDot <-\alpha V
    %     alpha=Param.alpha;
    % else
    %     alpha = 0;
    % end
    
    Abar = pwasys.Abar;
    Bbar = pwasys.Bbar;
    Ebar = pwasys.Ebar;
    Fbar = pwasys.Fbar;
    Qm=Q^.5;
    Rm=R^.5;
    
    for i=1:NR,
        for j=1:NS,
            A{i,j} = Abar{i,j}(1:end-1,1:end-1);
            a{i,j} = Abar{i,j}(1:end-1,end);
            B{i,j} = Bbar{i,j}(1:end-1,:);
        end
        Echeck{i} = [zeros(1,n) 1; Ebar{i}];
    end
    A1=A{istar,1};
    A2=A{istar,2};
    
    B1=B{istar,1};
    B2=B{istar,2};
    
    S=sdpvar(n,n);
    Y=sdpvar(m,n, 'full');
    
    %%%%%%% matrix S>0 and matrix T_jad_i should be nagative definite
    T_jad_1=[S*A1'+A1*S+Y'*B1'+B1*Y   S*Qm  Y'*R^.5;
        Qm*S           -eye(n)      zeros(n,m);
        Rm*Y       zeros(m, n)        -eye(m)];
    
    T_jad_2=[S*A2'+A2*S+Y'*B2'+B2*Y   S*Qm  Y'*R^.5;
        Qm*S             -eye(n)     zeros(n,m);
        Rm*Y         zeros(m, n)    -eye(m)];
    
    
    constraints= set (S>0);
    constraints=constraints + set(T_jad_1 < 0) + set (T_jad_2 <0);
    obj= -trace(S);
    solvesdp(constraints,obj,sdpsettings('usex0',1, 'shift', 1e-7));
    [p, d]=checkset(constraints);
    S=double(S);
    Y=double(Y);
    K=Y*inv(S);
    
    
    
    
end
