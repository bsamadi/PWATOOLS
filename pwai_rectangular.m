function pwainc = pwai_rectangular(nlfun)

% pwai_rectangular Piecewise-affine inclusion for a rectangular grid.
%
% pwainc = pwa_rectangular(nlfun)
%
% nlfun.Handle = @function_name;
% nlfun.Domain = [xmin xmax];
% nlfun.GR = Uniform grid resolution;
% nlfun.Resolution = Fine grid resolution; % (Optional)
% nlfun.ObjFun = 'L2'; % or 'Linf' (Optional)
%
% pwainc.Abar = Linear coefficients; % Ybar = Abar{i}*Xbar, Xbar = [X;1]
% pwainc.Err = Approximation error for the fine mesh
% pwainc.W = Fine mesh points
% pwainc.Z = The value of pwainc at W
% pwainc.Obj = The value of the objective function
%
% Copyright: Behzad Samadi, Amirkabir University of Technology, Oct. 2008

% Argument checking
f = nlfun.Handle;                     % Specify the handle of the nonlinear function
n = length(nlfun.Domain);             % Number of variables in the domain of nonlinearities

a=[];
for i=1:n,
    a=[a; nlfun.Domain{i}(1)];        % Set start point of the domain
end

if isfield(nlfun,'Parameters'),
    Param = nlfun.Parameters;         % Set the system parameters
else
    Param = [];
end

if isfield(nlfun,'ObjFun'),
    Option = nlfun.ObjFun;            % Specify the objective (L2 or Linf)
end

%
for i = 1:n,
    if any(diff(nlfun.Domain{i})<0),
        eval(['error(''nlfun.domain{' num2str(i) '} should be increasing.'');'])
    end
end

%
try
    Y = pwa_nlfun(nlfun,a,Param);
catch
    error(['The function can not be evaluated.']);
end

% The Approximation

% Construct the fine mesh
%
% Z = f(W)
% W = [w1; w2; ... ; wm] is m times n
% Z = [z1; z2; ... ; zm] is m times p
%
% zi = f(wi)

% fine mesh
for i=1:n,
    Domain{i} = nlfun.Domain{i}([1 end]);
end
W = pwa_grid(Domain,nlfun.Resolution);
GRP = pwa_grid(Domain,nlfun.GR);
W = union(GRP, W,'rows');    % Add gridding points to the data set
Z = pwa_nlfun(nlfun,W,Param);

% The Approximation
q = size(Z,2);                      % Number of nonlinear functions
Ind = pwa_index(nlfun.GR);
NR = size(Ind,2);                   % Number of regions

% Upper bound

for i=1:NR,
    Abar{i} = sdpvar(q,n+1,'full'); % PWA approximation parameters;
end

p = size(W,1);                      % Number of data points
Index = zeros(size(Ind,1),p);
Err = [];

constraints = set([]);

for i = 1:p,
    Wi = W(i,:)';
    R{i} = pwa_recreg(Wi,nlfun);
    
    for j = 1:length(R{i}),
        Ri = cell2mat(R{i});
        Zhat = (Abar{Ri(j)}*[Wi; 1])';
        Err = [Err;Z(i,:)-Zhat];
        constraints=constraints+set(Z(i,:)<Zhat,'Z<Zhat');
    end
end

% The equilibrium point
Ro = pwa_recreg(nlfun.xcl,nlfun);
Ro = cell2mat(Ro);
Zo = pwa_nlfun(nlfun,nlfun.xcl',Param);

for j = 1:length(Ro),
    Zhat = (Abar{Ro(j)}*[nlfun.xcl; 1])';
    constraints=constraints+set(Zo==Zhat,'Z==Zhat');
end
clear R;

% Objective
if ~isfield(Option,'ObjFun'),
    Option.ObjFun='L2';
end


if strcmp(Option.ObjFun,'L2'),
    %     Obj=sum(sum(Err.*Err));
    Obj = sdpvar(1,1);
    ErrVec = reshape(Err,numel(Err),1);
    Schur = [Obj ErrVec';ErrVec eye(length(ErrVec))];
    constraints=constraints+set(Schur>0,'Shur>0');
elseif strcmp(Option.ObjFun,'Linf'),
    Obj = sdpvar(1,1);
    constraints=constraints+set(Obj>0,'Obj>0');
    constraints=constraints+set(Err>-Obj,'Error>-Obj');
    constraints=constraints+set(Err<Obj,'Error<Obj');
end

% Solve the optimization problem
solution=solvesdp(constraints,Obj,sdpsettings('Shift',1e-6))
checkset(constraints)

for i=1:length(Abar),
    Abar{i} = double(Abar{i});
end
Obj = double(Obj);
Err = double(Err);

GridPoints = cell(n,1);
for i=1:n,
    GridPoints{i}=linspace(Domain{i}(1),Domain{i}(end),nlfun.GR(i)+1);
end

pwainc.Abar = Abar';
pwainc.Err{1} = Err;
pwainc.Obj{1} = Obj;
pwainc.W = W;
pwainc.Z = Z;
pwainc.Domain = Domain;
pwainc.GR = nlfun.GR;
pwainc.Ind = Ind;
pwainc.NR = NR;
pwainc.GridPoints = GridPoints;

% Lower bound
for i=1:NR,
    Abar{i} = sdpvar(q,n+1,'full'); % PWA approximation parameters;
end
Err = [];

constraints = set([]);

for i = 1:p,
    Wi = W(i,:)';
    R{i} = pwa_recreg(Wi,nlfun);
    for j = 1:length(R{i}),
        Ri = cell2mat(R{i});
        Zhat = (Abar{Ri(j)}*[Wi; 1])';
        Err = [Err;Z(i,:)-Zhat];
        constraints=constraints+set(Z(i,:)>Zhat,'Z>Zhat');
    end
end

pwainc.Rw = R';

% The equilibrium point
Ro = pwa_recreg(nlfun.xcl,nlfun);
Ro = cell2mat(Ro);
Zo = pwa_nlfun(nlfun,nlfun.xcl',Param);

for j = 1:length(Ro),
    Zhat = (Abar{Ro(j)}*[nlfun.xcl; 1])';
    constraints=constraints+set(Zo==Zhat,'Z==Zhat');
end



% Objective
if ~isfield(Option,'ObjFun'),
    Option.ObjFun='L2';
end

if strcmp(Option.ObjFun,'L2'),
    %     Obj=sum(sum(Err.*Err));
    Obj = sdpvar(1,1);
    ErrVec = reshape(Err,numel(Err),1);
    Schur = [Obj ErrVec';ErrVec eye(length(ErrVec))];
    constraints=constraints+set(Schur>0,'Shur>0');
elseif strcmp(Option.ObjFun,'Linf'),
    Obj = sdpvar(1,1);
    constraints=constraints+set(Obj>0,'Obj>0');
    constraints=constraints+set(Err>-Obj,'Error>-Obj');
    constraints=constraints+set(Err<Obj,'Error<Obj');
end

% Solve the optimization problem
solution=solvesdp(constraints,Obj,sdpsettings('Shift',1e-6))
checkset(constraints)

for i=1:length(Abar),
    pwainc.Abar{i,2} = double(Abar{i});
end

% Abar = [A a;0 1]

Obj = double(Obj);
Err = double(Err);

pwainc.Err{2} = Err;
pwainc.Obj{2} = Obj;