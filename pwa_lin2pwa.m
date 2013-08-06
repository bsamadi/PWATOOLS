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

[NR NS] = size(pwasys.Abar);    % Number of Systems, Number of Regions
n = size(pwasys.Abar{1},1)-1;   % Number of state variables
m = size(pwasys.Bbar{1},2);     % Number of inputs

if isfield(Param,'LimitK'),
    LimitK = Param.LimitK;
end

if isfield(Param,'LimitU'),
    LimitU = Param.LimitU;
end

if isfield(Param,'InSat'),
    InSat = Param.InSat;
end

if isfield(Param,'Lyapunov'),
    Lyapunov = Param.Lyapunov;
else
    Lyapunov = 'PWQ';
end

if isfield(Param,'L2Gain'),
    L2Gain = Param.L2Gain;
    gamma = Param.L2Gain^2;
    Bwbar = pwasys.Bwbar;
    Cbar = pwasys.Cbar;
    Dw = pwasys.Dw;
    if isfield(pwasys,'Du'),
        Du = pwasys.Du;
    else,
        for i=1:NR,
            Du{i} = zeros(size(Cbar{1},1),m);
        end
    end
    for i=1:NR,
        Bw{i} = Bwbar{i}(1:end-1,:);
        C{i} = Cbar{i}(:,1:end-1);
    end
    q = size(pwasys.Bwbar{1},2);    % Number of disturbance inputs
end

if isfield(Param,'Poles'),
    Poles = Param.Poles;
elseif isfield(Param,'KLin'),
    KLin = Param.KLin;
elseif isfield(Param,'QLin'),
    QLin = Param.QLin;
    RLin = Param.RLin;
end

if isfield(Param,'alpha'),
    alpha=Param.alpha;
else
    alpha = 0;
end

istar = [];
Abar = pwasys.Abar;
Bbar = pwasys.Bbar;
Ebar = pwasys.Ebar;
Fbar = pwasys.Fbar;

for i=1:NR,
    for j=1:NS,
        A{i,j} = Abar{i,j}(1:end-1,1:end-1);
        a{i,i} = Abar{i,j}(1:end-1,end);
        B{i,j} = Bbar{i,j}(1:end-1,:);
    end
    Echeck{i} = [zeros(1,n) 1; Ebar{i}];
    xcl_is_inside_Ri = all(Ebar{i}*xclbar>=0);
    if xcl_is_inside_Ri,
        istar = union(istar,i);
    end
end

Lin2PWA = 1;
for i = istar,
    if exist('Poles'),
        K{i} = -acker(A{i},B{i},Poles);
    elseif exist('KLin'),
        K{i} = KLin;
    elseif exist('QLin'),
        K{i} = -lqr(A{i},B{i},QLin,RLin);
    else
        Lin2PWA = 0;
    end
end

yalmip('clear');
constraints=set([]);
p = size(Ebar{1},1);

% Lyapunov function
if strcmp(Lyapunov,'Global'),
    Pg = sdpvar(n,n);
    for i=1:NR,
        P{i} = Pg;
        r{i} = 0;
        % Constrainted Pbar
        Pbar{i} = [P{i} -P{i}*xcl; -xcl'*P{i} xcl'*P{i}*xcl+r{i}];
    end
else,
    for i=1:NR,
        if ismember(i,istar),
            P{i} = sdpvar(n,n);
            Pbar{i} = [P{i} -P{i}*xcl; -xcl'*P{i} xcl'*P{i}*xcl];
        else
            Pbar{i} = sdpvar(n+1,n+1);
        end
    end
end

for i=1:NR,
    Kbar{i} = sdpvar(m,n+1,'full');
end

%%%%%%%%%%%%%%%%%%%%%%%% The center region %%%%%%%%%%%%%%%%%%%
I = eye(n);
One = zeros(n+1);
One(end,end)=1;

if exist('InSat'),
    lambda = sdpvar(1,1);
    constraints = constraints+set(lambda>0,'lambda>0');
end

if Lin2PWA,
    for i=istar,
        k = sdpvar(m,1);
        Kbar{i} = [K{i} k];
    end
end

for i=istar,
    K{i} = Kbar{i}(:,1:n);
    for j=1:NS,
        constraints=constraints+set( (Abar{i,j}+Bbar{i,j}*Kbar{i})*[xcl;1] == 0,['(Abar+Bbar*Kbar)*xclbar=0']);
        if exist('L2Gain'),
            L2Mat{i} = [P{i}*(A{i}+B{i}*K{i})+(A{i}+B{i}*K{i})'*P{i}+(C{i}+Du{i}*K{i})'*(C{i}+Du{i}*K{i}) P{i}*Bw{i}+(C{i}+Du{i}*K{i})'*Dw{i}; (P{i}*Bw{i}+(C{i}+Du{i}*K{i})'*Dw{i})' Dw{i}'*Dw{i}-gamma*eye(q)];
            constraints=constraints+set( L2Mat{i} < 0,['L2Mat' num2str(i) '<0']);
        else,
            if Lin2PWA,
                DV{i,j} = Pbar{i}*(Abar{i,j}+Bbar{i,j}*Kbar{i})+(Abar{i,j}+Bbar{i,j}*Kbar{i})'*Pbar{i};  % Don't use alpha for the center region(s) if you already have the controller
            else
                DV{i,j} = Pbar{i}*(Abar{i,j}+Bbar{i,j}*Kbar{i})+(Abar{i,j}+Bbar{i,j}*Kbar{i})'*Pbar{i}+alpha*Pbar{i};                
            end
            constraints=constraints+set( DV{i,j}(1:n,1:n) < 0,['DV' num2str(i) ',' num2str(j) '<0']);
        end
    end
    if exist('InSat'),
        constraints=constraints+set(P{i} > 0.075*I,['P' num2str(i) '>0']);
        constraints=constraints+set(P{i} < lambda*I,['P' num2str(i) '<lambda*I']);
    else
        constraints=constraints+set(P{i} > 0.01*I,['P' num2str(i) '>0.01*I']);
    end
end



for i=1:NR,
    %%%%%%%%%%%%%%%%%%%%%%%% Inequality Constraints %%%%%%%%%%%%%%%%%%%
    if strcmp(Lyapunov,'PWQ'),
        if i~=istar,
            if ~exist('InSat'),
                Z{i} = sdpvar(p+1,p+1);
                % S procedure
                constraints=constraints+set(0 < Z{i}(:) < 1e5,['Z{' num2str(i) '}(:)>0']);
                constraints=constraints+set(Pbar{i}-Echeck{i}'*Z{i}*Echeck{i} > 0,['P' num2str(i) '>0']);
            end
        end
    end %PWQ%

    if exist('LimitK'),
        if i~=istar | ~Lin2PWA,
            constraints=constraints+set(Kbar{i}(:) < LimitK,'Kbar<LimitK');
            constraints=constraints+set(Kbar{i}(:) > -LimitK,'Kbar>-LimitK');
        end
    end %LimitK%

    if exist('LimitU'),
        LinVar = find(pwasys.NonlinearDomain==0);
        NLinVar = find(pwasys.NonlinearDomain==1);
        if i~=istar | ~Lin2PWA,
            for j=1:size(pwasys.T,2),
                for ik = 0:2^length(LinVar)-1,
                    Index = [];
                    Bin = dec2bin(ik);
                    if length(Bin)<length(LinVar),
                        Bin = ['0'*ones(1,length(LinVar)-length(Bin)) Bin];
                    end
                    for jk = 1:length(LinVar),
                        if jk<=length(Bin),
                            Index(jk) = str2num(Bin(jk));
                        else,
                            Index(jk) = 0;
                        end
                    end
                    Index = Index + 1;
                    V = zeros(size(pwasys.NonlinearDomain));
                    for jk = 1:length(LinVar),
                        V(LinVar(jk)) = pwasys.Domain{LinVar(jk)}(Index(jk));
                    end
                    V(NLinVar) = pwasys.X(pwasys.T(i,j),:);
                    constraints=constraints+set(Kbar{i}*[V; 1] < LimitU(:,2),'U<LimitU(2)');
                    constraints=constraints+set(Kbar{i}*[V; 1] > LimitU(:,1),'U>LimitU(1)');
                end
            end
        end
    end %LimitU%

    if exist('InSat'),
        umax = InSat.Limit;
        if strcmp(InSat.Method,'Ellipsoid'),
            Sp{i} = sdpvar(p+1,p+1);
            ESp{i} = sdpvar(1);
            % S procedure
            constraints = constraints+set( 0 < Sp{i}(:) < 1e5,['Sp{' num2str(i) '}(:)>0']);
            constraints = constraints+set( 0 < ESp{i}(:) < 1e5,['ESp{' num2str(i) '}(:)>0']);
            for j=1:m,
                Kjbar = Kbar{i}(j,:);
                constraints=constraints+set(Kjbar'*Kjbar+ESp{i}*(-Pbar{i}+One)+Echeck{i}'*Sp{i}*Echeck{i}<umax*One,'Sat');
            end
        end
    end %InSat%

    %%%%%%%%%%%%%%%%%%%%%%%%% Negative definite %%%%%%%%%%%%%%%%%%%%%%%%%
    if i~=istar,
        for j=1:NS,
            L{i} = sdpvar(p+1,p+1);
            % S procedure
            constraints=constraints+set(0 < L{i}(:) < 1e5,['L{' num2str(i) '}(:)>0']);
            if exist('L2Gain'),
                L2Mat{i} = [Pbar{i}*(Abar{i}+Bbar{i}*Kbar{i})+(Abar{i}+Bbar{i}*Kbar{i})'*Pbar{i}+(Cbar{i}+Du{i}*Kbar{i})'*(Cbar{i}+Du{i}*Kbar{i})+Echeck{i}'*L{i}*Echeck{i} Pbar{i}*Bwbar{i}+(Cbar{i}+Du{i}*Kbar{i})'*Dw{i}; (Pbar{i}*Bwbar{i}+(Cbar{i}+Du{i}*Kbar{i})'*Dw{i})' Dw{i}'*Dw{i}-gamma*eye(q)];
                constraints=constraints+set( L2Mat{i} < 0,['L2Mat' num2str(i) '<0']);
            else
                DV{i,j} = Pbar{i}*(Abar{i,j}+Bbar{i,j}*Kbar{i})+(Abar{i,j}+Bbar{i,j}*Kbar{i})'*Pbar{i}+alpha*Pbar{i}+Echeck{i}'*L{i}*Echeck{i};
                constraints=constraints+set( DV{i,j} < 0,['DV' num2str(i) '-' num2str(j) '<0']);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%% Continuity constraints %%%%%%%%%%%%%%%%%%%%%%%%
    for j = i:NR,
        if ~isempty(Fbar{i,j}),
            % Continuity of the control input
            constraints=constraints+set((Kbar{i}-Kbar{j})*Fbar{i,j}== 0,['Ki*F' num2str(i) '_' num2str(j) 's=Kis*F' num2str(i) '_' num2str(j) 's']);
            if strcmp(Lyapunov,'PWQ'),
                % Continuity of the Lyapunov function
                constraints=constraints+set(Fbar{i,j}'*(Pbar{i}-Pbar{j})*Fbar{i,j} == 0,['F' num2str(i) '_' num2str(j) 's*Pbari*F' num2str(i) '_' num2str(j) 's=F' num2str(i) '_' num2str(j) 's*Pbaris*F' num2str(i) '_' num2str(j) 's']);
            end
        end
    end
end

if Lin2PWA,
    for i = 1:NR,
        if i~=istar,
            if isfield(pwasys,'Kbar'),
                setsdpvar(Kbar{i},pwasys.Kbar{i});
            else
                setsdpvar(Kbar{i}(:,1:end-1),Kbar{istar(1)}(:,1:end-1));
            end
        end
        if isfield(pwasys,'Pbar'),
            setsdpvar(Pbar{i},pwasys.Pbar{i});
        end
    end
end

if exist('InSat'),
    solution=solvesdp(constraints,lambda,sdpsettings('usex0',1,'shift',1e-7))
else
    solution=solvesdp(constraints,[],sdpsettings('usex0',1,'shift',1e-7))
end
checkset(constraints);
pwactrl.constraints = constraints;

pwactrl.xcl = xcl;
pwactrl.istar = istar;

for i=1:NR,
    pwactrl.Kbar{i} = double(Kbar{i});
    pwactrl.Pbar{i} = double(Pbar{i});
end

pwactrl.Index = pwasys.Index;
pwactrl.X = pwasys.X;
pwactrl.T = pwasys.T;
pwactrl.xcl = pwasys.xcl;
pwactrl.Ebar = pwasys.Ebar;