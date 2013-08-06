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

if isfield(Param,'LimitK'),
    LimitK = Param.LimitK;
end

if isfield(Param,'Poles'),
    Poles = Param.Poles;
elseif isfield(Param,'KLin'),
    KLin = Param.KLin;
elseif isfield(Param,'QLin'),
    QLin = Param.QLin;
    RLin = Param.RLin;
end

if isfield(Param,'LimitU'),
    LimitU = Param.LimitU;
end

alpha=Param.alpha;
xcl = Param.xcl;

M = length(pwasys.Abar);        % Number of regions
xclbar = [Param.xcl;1];
n = size(pwasys.Abar{1},1)-1;   % Number of state variables
m = size(pwasys.Bbar{1},2);     % Number of inputs

istar = [];
Abar = pwasys.Abar;
Bbar = pwasys.Bbar;
Ebar = pwasys.Ebar;
Fbar = pwasys.Fbar;

for i=1:M,
    A{i} = Abar{i}(1:end-1,1:end-1);
    a{i} = Abar{i}(1:end-1,end);
    B{i} = Bbar{i}(1:end-1,:);
    Echeck{i} = [zeros(1,n) 1; Ebar{i}];
    xcl_is_inside_Ri = all(Ebar{i}*xclbar>=0);
    if xcl_is_inside_Ri,
        istar = union(istar,i);
    end
end

istar = istar(end);
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
end

for i=istar
    eig(A{i}+B{i}*K{i})
end

Or = istar;
CR = istar;
i=istar(1);

Kstar = K{istar};
Kbarstar = [Kstar -Kstar*xcl-B{i}\(A{i}*xcl+a{i})];
Kbar{i} = Kbarstar;

yalmip('clear');
constraints=set([]);
obj = sdpvar(1,1);
constraints = constraints+set(obj>0,'Obj>0');
P{i} = sdpvar(n,n);
r{i} = sdpvar(1,1);
Pbar{i} = [P{i} -P{i}*xcl; -xcl'*P{i} xcl'*P{i}*xcl+r{i}];
constraints=constraints+set(P{i} > 0,['P' num2str(i) '>0']);

N=[];
while length(CR)<M,

    % Neighbors of central region
    N=[];
    for i=1:M,
        for j=1:length(CR),
            if ~isempty(pwasys.Fbar{i,CR(j)}) & ~ismember(i,CR),
                N = [N i];
            end
        end
    end

    for i=N,
        p = size(pwasys.Ebar{i},1);
        %%%%%%%%%%%%%%%%%%%%%%%% Declare Variables %%%%%%%%%%%%%%%%%%%%%%%%
        % Constrainted Pbar
        P{i} = sdpvar(n,n);
        r{i} = sdpvar(1,1);
        Pbar{i} = [P{i} -P{i}*xcl; -xcl'*P{i} xcl'*P{i}*xcl+r{i}];
        constraints=constraints+set(P{i} > 0,['P' num2str(i) '>0']);
        Kbar{i} = sdpvar(m,n+1,'full');
        L{i} = sdpvar(p+1,p+1,'symmetric');
        %%%%%%%%%%%%%%%%%%%%%%%% Inequality Constraints %%%%%%%%%%%%%%%%%%%%%%%%
        if exist('LimitK'),
            constraints=constraints+set(Kbar{i}(:) < LimitK,'Kbar<LimitK');
            constraints=constraints+set(Kbar{i}(:) > -LimitK,'Kbar>-LimitK');
        end

        if exist('LimitU'),
            LinVar = find(pwasys.NonlinearDomain==0);
            NLinVar = find(pwasys.NonlinearDomain==1);
            if i~=istar,
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
                            V(LinVar(jk)) = pwasys.Domain{jk}(Index(jk));
                        end
                        V(NLinVar) = pwasys.X(pwasys.T(i,j),:);
                        constraints=constraints+set(Kbar{i}*[V; 1] < LimitU(:,2),'U<LimitU(2)');
                        constraints=constraints+set(Kbar{i}*[V; 1] > LimitU(:,1),'U>LimitU(1)');
                    end
                end
            end
        end

        constraints=constraints+set(L{i}(:) > 0,'L(:)>0');
        %%%%%%%%%%%%%%%%%%%%%%%%% Positive definite %%%%%%%%%%%%%%%%%%%%%%%%%
        DV{i} = Pbar{i}*(Abar{i}*Bbar{i}*Kbar{i})+(Abar{i}*Bbar{i}*Kbar{i})'*Pbar{i}+alpha*Pbar{i}+Echeck{i}'*L{i}*Echeck{i};
        constraints=constraints+set( DV{i} < 0,['DV' num2str(i) '<0']);
        %%%%%%%%%%%%%%%%%%%%%%%% Equality Constraints %%%%%%%%%%%%%%%%%%%%%%%%
        constraints = constraints+set((Abar{i}+Bbar{i}*Kbar{i}) - (Abar{istar}+Bbar{istar}*Kbar{istar}) <  obj ,'A+BK < cte');
        constraints = constraints+set((Abar{i}+Bbar{i}*Kbar{i}) - (Abar{istar}+Bbar{istar}*Kbar{istar}) > -obj ,'A+BK > -cte');
        for j = union(CR,N),
            if ~isempty(Fbar{i,j}),
                constraints=constraints+set((Kbar{i}-Kbar{j})*Fbar{i,j}== 0,'Ki*Fi_is=Kis*Fi_is');
                constraints=constraints+set(Fbar{i,j}'*(Pbar{i}-Pbar{j})*Fbar{i,j} == 0,'Fi_is*Pbari*Fi_is=Fi_is*Pbaris*Fi_is');
            end
        end
    end
    CR = union(CR,N);
end

for i = 1:length(Kbar),
    if i~=istar,
        setsdpvar(Kbar{i}(1:end-1),Kbar{istar}(1:end-1));
    end
end

solution=solvesdp(constraints,obj,sdpsettings('usex0',1,'shift',1e-7));
yalmiperror(solution.problem);
checkset(constraints);

obj = double(obj)

% pwactrl.Q = Q;
pwactrl.xcl = xcl;
pwactrl.istar = istar;

for i=1:M,
    pwactrl.Kbar{i} = double(Kbar{i});
    pwactrl.Pbar{i} = double(Pbar{i});
end
pwactrl.Index = pwasys.Index;
pwactrl.X = pwasys.X;
pwactrl.T = pwasys.T;
pwactrl.xcl = pwasys.xcl;
pwactrl.Ebar = pwasys.Ebar;
pwactrl.obj = obj;