function pwactrl = pwa_sos(pwasys, Param)

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

if isfield(Param,'LimitU'),
    LimitU = Param.LimitU;
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

xcl = pwasys.xcl;
xclbar = [xcl;1];

NR = length(pwasys.Abar);       % Number of regions
n = size(pwasys.Abar{1},1)-1;   % Number of state variables
m = size(pwasys.Bbar{1},2);     % Number of inputs

istar = [];
Abar = pwasys.Abar;
Bbar = pwasys.Bbar;
Ebar = pwasys.Ebar;
Fbar = pwasys.Fbar;

for i=1:NR,
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
Lin2PWA = 1;
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
else
    Lin2PWA = 0;
end

yalmip('clear');
constraints=set([]);
obj = sdpvar(NR,1);
X = sdpvar(n,1);
Xbar = [X;1];
p = size(Ebar{1},1);
for i=1:NR,
    Kbar{i} = sdpvar(m,n+1,'full');
    L{i} = sdpvar(p+1,p+1);

    % Constrainted Pbar
    P{i} = sdpvar(n,n);
    r{i} = sdpvar(1,1);
    Pbar{i} = [P{i} -P{i}*xcl; -xcl'*P{i} xcl'*P{i}*xcl+r{i}];
end

i=istar(1);
I = eye(n+1);

if Lin2PWA,
    k = sdpvar(m,1);
    Kbar{i} = [K{i} k];
end
constraints=constraints+set( (Abar{i}+Bbar{i}*Kbar{i})*[xcl;1] == 0,['(Abar+Bbar*Kbar)*xclbar=0']);

DV{i} = Pbar{i}*(Abar{i}+Bbar{i}*Kbar{i})+(Abar{i}+Bbar{i}*Kbar{i})'*Pbar{i}+alpha*Pbar{i};
constraints=constraints+set( sos(-X'*DV{i}(1:n,1:n)*X) ,['DV' num2str(i) '<0']);

for i=1:NR,
    %%%%%%%%%%%%%%%%%%%%%%%% Declare Variables %%%%%%%%%%%%%%%%%%%%%%%%
    constraints=constraints+set(Pbar{i} > I,['P' num2str(i) '>I']);
    constraints=constraints+set(Pbar{i} < obj(i)*I,['P' num2str(i) '<obj*I']);
    %%%%%%%%%%%%%%%%%%%%%%%% Inequality Constraints %%%%%%%%%%%%%%%%%%%%%%%%
    if exist('LimitK'),
        if i~=istar | ~Lin2PWA,
            constraints=constraints+set(Kbar{i}(:) < LimitK,'Kbar<LimitK');
            constraints=constraints+set(Kbar{i}(:) > -LimitK,'Kbar>-LimitK');
        end
    end

    if exist('LimitU'),
        if i~=istar | ~Lin2PWA,
            for j=1:size(pwasys.T,2),
                constraints=constraints+set(Kbar{i}*[pwasys.X(pwasys.T(i,j),:) 1]' < LimitU(2),'U<LimitU(2)');
                constraints=constraints+set(Kbar{i}*[pwasys.X(pwasys.T(i,j),:) 1]' > LimitU(1),'U>LimitU(1)');
            end
        end
    end

    % S procedure
    constraints=constraints+set(L{i}(:) > 0,['L{' num2str(i) '}(:)>0']);
    %%%%%%%%%%%%%%%%%%%%%%%%% Positive definite %%%%%%%%%%%%%%%%%%%%%%%%%
    if i~=istar,
        DV{i} = Pbar{i}*(Abar{i}+Bbar{i}*Kbar{i})+(Abar{i}+Bbar{i}*Kbar{i})'*Pbar{i}+alpha*Pbar{i}+Echeck{i}'*L{i}*Echeck{i};
        constraints=constraints+set( sos(-Xbar'*DV{i}*Xbar),['DV' num2str(i) '<0']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%% Equality Constraints %%%%%%%%%%%%%%%%%%%%%%%%
    for j = i:NR,
        if ~isempty(Fbar{i,j}),
            constraints=constraints+set((Kbar{i}-Kbar{j})*Fbar{i,j}== 0,['Ki*F' num2str(i) '_' num2str(j) 's=Kis*F' num2str(i) '_' num2str(j) 's']);
            constraints=constraints+set(Fbar{i,j}'*(Pbar{i}-Pbar{j})*Fbar{i,j} == 0,['F' num2str(i) '_' num2str(j) 's*Pbari*F' num2str(i) '_' num2str(j) 's=F' num2str(i) '_' num2str(j) 's*Pbaris*F' num2str(i) '_' num2str(j) 's']);
        end
    end
end

if Lin2PWA,
    for i = 1:NR,
        if i~=istar,
            setsdpvar(Kbar{i}(1:end-1),Kbar{istar}(1:end-1));
        end
    end
end

solution=solvesos(constraints,obj'*obj,sdpsettings('usex0',1,'shift',1e-7))
checkset(constraints);

obj = double(obj);

pwactrl.xcl = xcl;
pwactrl.istar = istar;

% To make V(0) = 0
V0 = double(Pbar{istar}(end,end));
Pbar0 = zeros(n+1,n+1);
Pbar0(end,end) = V0;

for i=1:NR,
    pwactrl.Kbar{i} = double(Kbar{i});
    pwactrl.Pbar{i} = double(Pbar{i})-Pbar0;
end

pwactrl.Index = pwasys.Index;
pwactrl.X = pwasys.X;
pwactrl.T = pwasys.T;
pwactrl.xcl = pwasys.xcl;
pwactrl.Ebar = pwasys.Ebar;
pwactrl.obj = obj;