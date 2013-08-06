function pwactrl = pwai_ctrlsyncvx(pwainc)

% PWA controller synthesis for PWA differential inclusions
%
% Input: pwainc
% xDot \in \conv{ pwainc.A{i}x+pwainc.b{i} } when x\in Region i
% Region i= { x | |q_{ij}| < 1 }
% q_{ij} = Cq{i,j} x+ dq{i,j}
%
% Controller: pwactrl
% u = pwactrl.K{i}*x+pwactrl.k{i}

n = size(pwainc.A{1,1},1);   % Number of state variables
nR = length(pwainc.Domain);  % Number of variables in the domain of nonlinearities
m = size(pwainc.B{1},2);     % Number of inputs
NR = pwainc.NR;              % Number of regions

Iz = cell2mat(pwa_recreg(zeros(n,1),pwainc)); % The regions that include the origin
I0 = zeros(1,length(Iz));
for i=1:length(Iz),
    I0(i) = Iz(i);
end
I1 = setdiff(1:NR,I0);      % Other regions

constraints=set([]);
Q = sdpvar(n,n);
Y = cell(NR,1);
Z = cell(NR,1);
W = cell(NR,1);
mu = cell(NR,1);

constraints = constraints+set(Q>0.5);
constraints = constraints+set(Q<2);
%Q = eye(2);
alpha = 0; % decay rate

for i = I0,
    Y{i} = sdpvar(m,n);
    for j=1:2,
        constraints = constraints+set(pwainc.A{i,j}*Q+Q*pwainc.A{i,j}'+pwainc.B{i}*Y{i}+Y{i}'*pwainc.B{i}'<-alpha*Q);
        constraints = constraints+set(-1000<Y{i}(:)<1000);
    end
end

for i=I1,
    Y{i} = sdpvar(m,n);
    Z{i} = [sdpvar(m,1) zeros(m,nR-1)];
    W{i} = zeros(m,m);
    mu{i} = sdpvar(nR,1);
    Mj{i} = diag(mu{i});
    Bu = pwainc.B{i};
    Bv = [pwainc.B{i} zeros(n,nR-1)];
    Cq = zeros(nR,n);
    Dq = zeros(nR,nR);
    for j=1:size(pwainc.Cq,2),
        Cq(j,:)=pwainc.Cq{i,j};
        Dq(j,j) = pwainc.dq{i,j};
    end
    constraints = constraints+set(mu{i}<0);
    constraints = constraints+set(-1000<Y{i}(:)<1000);
    constraints = constraints+set(-1000<Z{i}(:)<1000);
    for j=1:size(pwainc.A,2),
        Omega11 = pwainc.A{i,j}*Q+Q*pwainc.A{i,j}'+Bu*Y{i}+Y{i}'*Bu'+Bv*Mj{i}*Bv'+Bv*Z{i}'*Bu'+Bu*Z{i}*Bv'+Bu*W{i}*Bu';
        Omega21 = Cq*Q+Dq*Mj{i}*Bv'+Dq*Z{i}'*Bu';
        Omega22 = -Mj{i}+Dq*Mj{i}*Dq';
        constraints = constraints+set([Omega11 Omega21';Omega21 Omega22]<0);
    end
end

solution=solvesdp(constraints,[],sdpsettings('usex0',1,'shift',1e-7));
checkset(constraints);

pwactrl.Q = double(Q);
for j=1:NR,
    pwactrl.Y{j} = double(Y{j});
    pwactrl.Z{j} = double(Z{j});
    pwactrl.mu{j} = double(mu{j});
    pwactrl.K{j} = pwactrl.Y{j}/pwactrl.Q;
end
for j=I1,
    pwactrl.k{j} = double(Mj{j}(1,1))\double(Z{j}(:,1));
end
for j=I0,
    pwactrl.k{j} = zeros(m,1);
end

for i=1:NR,
    pwactrl.Kbar{i} = [pwactrl.K{i} pwactrl.k{i}];
end

pwactrl.Domain = pwainc.Domain;
pwactrl.GR = pwainc.GR;