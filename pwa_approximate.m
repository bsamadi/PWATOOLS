function pwasys = pwa_approximate(nlsys)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Create a piecewise-affine (PWA) system pwasys from the nonlinear      %
%   system nlsys. First, the nonlinear function nlfun is extracted from   %
%   nlsys. Then, nlfun is approximated by a PWA function pwafun. And      %
%   finally, pwafun is used to create the overall PWA system pwasys.      %
%                                                                         %
%   Input:                                                                %
%       nlsys                                                             %
%                                                                         %
%   Output:                                                               %
%       pwasys                                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if  ~isfield(nlsys,'NonlinearFunction'),
    nlsys.NonlinearFunction = inline('0*x;');
end

% Specify the handle of the nonlinear function
nlfun.Handle = nlsys.NonlinearFunction;

% Define the domain of nonlinearity and an index of variables in the domain
nlfun.Domain = [];
nlfun.Index = [];
Domain = nlsys.Domain;

if  isfield(nlsys,'NonlinearDomain'),
    index = find(nlsys.NonlinearDomain==0);
    for i = index',
        nlsys.Domain{i}=[];
    end
end
 
j=1;
for i = 1:length(nlsys.Domain),
    if ~isempty(nlsys.Domain{i}),
        nlfun.Index = [nlfun.Index i];
        nlfun.Domain{j} = nlsys.Domain{i};
        j=j+1;
    end
end

if ~isfield(nlsys,'xcl'),
    error('Closed loop equilibrium point (nlsys.xcl) is not defined.');
else,
     nlfun.xcl = nlsys.xcl(nlfun.Index);
end

% Set the dimension of the nonlinear function
nlfun.Dim = length(nlfun.Domain);
% Obtain indices of state variables described by nonlinear functions
YIndex = find(nlsys.NonlinearStateEquations);
% Obtain indices of state variables present in domain of nonlinearity
XIndex = nlfun.Index;
nlfun.YIndex = YIndex;

% Set the system parameters
if isfield(nlsys,'Parameters')
    nlfun.Parameters = nlsys.Parameters;
end
% Set the sampling resolution
nlfun.Resolution = nlsys.Resolution;
% Set the objective (L2 or Linf)
if isfield(nlsys,'ObjFun')
    nlfun.ObjFun = nlsys.ObjFun;
end

if ~isfield(nlsys,'Method'),
    nlsys.Method = 'Uniform';
end

% Set the resolution of the uniform gird
if isfield(nlsys,'UGR'),
    nlfun.UGR = nlsys.UGR;
else
    nlfun.UGR = 5;
end

    
% Determine approximation method to use and compute the PWA approximation
if strcmpi(nlsys.Method,'uniform'),
    % Set the resolution of the uniform gird
    % Compute the approximation
    pwafun = pwa_uniform(nlfun);

elseif strcmpi(nlsys.Method,'optimaluniform'),
    % Anchored estimation
    if isfield(nlsys,'AbarLin'),
        nlfun.AbarLin = nlsys.AbarLin-[nlsys.A nlsys.a;zeros(1,length(nlsys.a)+1)];
        nlfun.AbarLin = nlfun.AbarLin(YIndex,[XIndex end]);
    end
    % Compute the approximation
    pwafun = pwa_optimal_uniform(nlfun);

elseif strcmpi(nlsys.Method,'multiresolution'),
    % Set the target number of regions
    nlfun.TNR = nlsys.TNR;

    % Anchored estimation
       if isfield(nlsys,'AbarLin'),
        nlfun.AbarLin = nlsys.AbarLin-[nlsys.A nlsys.a;zeros(1,length(nlsys.a)+1)];
        nlfun.AbarLin = nlfun.AbarLin(YIndex,[XIndex end]);
        if isfield(nlsys,'Rstar'),
            nlfun.Rstar = nlsys.Rstar;
        end
         end
    % Compute the approximation
    pwafun = pwa_split(nlfun);
elseif strcmpi(nlsys.Method,'upperbound'),
    % Set the resolution of the uniform gird
    % Compute the approximation
    nlfun.X = nlsys.X;
    Option.Type = 'UB';
    [m n] = size(nlfun.X);
    nlfun.xstar = nlfun.xcl; % The approximation should be exact at xcl
    if n>1,
        pwafun = pwa_bounds(nlfun,Option);
    else,
        pwafun = pwa_bounds_single(nlfun,Option);
    end

elseif strcmpi(nlsys.Method,'lowerbound'),
    % Set the resolution of the uniform gird
    % Compute the approximation
    nlfun.X = nlsys.X;
    Option.Type = 'LB';
    [m n] = size(nlfun.X);
    nlfun.xstar = nlfun.xcl; % The approximation should be exact at xcl
    if n>1,
        pwafun = pwa_bounds(nlfun,Option);
    else,
        pwafun = pwa_bounds_single(nlfun,Option);
    end
end

% Obtain the cell information
CellInfo = pwa_cell(pwafun);

% Obtain dimensions of Triangulation
[mT nT] = size(pwafun.T);

% Define the number of state variables
N = size(nlsys.A,1);

% Determine if B(x) is a matrix
ismatrix = 0;
try
    nlsys.B = 0*nlsys.B+nlsys.B;
    [mB nB] = size(nlsys.B);
    ismatrix = 1;
end

% Obtain indices of state variables not in the domain of nonlinearity
NotXIndex = setdiff(1:N,XIndex);

for i = 1:mT,
    if isempty(YIndex),
        pwasys.Abar{i} = [nlsys.A nlsys.a; zeros(1,N+1)];
    else
        % Create the auxiliary matrices for A and a
        Afun = zeros(N+1,N+1);
        Afun(YIndex,XIndex) = pwafun.Abar{i}(:,1:nlfun.Dim);
        Afun(YIndex,end) = pwafun.Abar{i}(:,end);
        pwasys.Abar{i} = [nlsys.A nlsys.a; zeros(1,N+1)]+Afun;
    end
    % Define augmented cell description matrices
    pwasys.Ebar{i} = zeros(nT,N+1);
    pwasys.Ebar{i}(:,XIndex) = CellInfo.Ebar{i}(:,1:nlfun.Dim);
    pwasys.Ebar{i}(:,end) = CellInfo.Ebar{i}(:,end);
    pwasys.Ebar{i} = [pwasys.Ebar{i};zeros(1,N) 1];
    
    % Define augmented boundary description matrices
    for j = 1:mT,
        if isempty(CellInfo.Fbar{i,j}),
            pwasys.Fbar{i,j} = [];
        else
            pwasys.Fbar{i,j} = zeros(N+1,N);
            pwasys.Fbar{i,j}(end,end) = 1;
            pwasys.Fbar{i,j}(NotXIndex,1:length(NotXIndex))=eye(length(NotXIndex));
            pwasys.Fbar{i,j}([XIndex end],length(NotXIndex)+(1:size(CellInfo.Fbar{i,j},2))) = CellInfo.Fbar{i,j};
        end
    end

    % Create B matrix
    if ismatrix,
        % If B(x) is already a matrix, use "as is"
        pwasys.Bbar{i} = [nlsys.B;zeros(1,size(nlsys.B,2))];
    else
        % If B(x) is not a matrix, find the Chebychedf center of the
        % polytopic region described by P={x| Ebar*xbar>0}
        % [Page 148 Convex optimization book by Boyd]

        yalmip('clear');
        r = sdpvar(1);
        xc = sdpvar(N,1);
        xcbar = [xc;1];
        constraints = set([]);
        for j = 1:nT,
            constraints = constraints + set(pwasys.Ebar{i}*xcbar > r*norm(pwasys.Ebar{i}(j,:)));
        end
        objective = -r;
        solution = solvesdp(constraints,[],objective,sdpsettings('verbose',0,'warning',0));
        xcheb = double(xc);
        NaN = find(isnan(xcheb));
        xcheb(NaN) = 0*NaN;
        if isfield(nlsys,'Parameters')
            B{i} = feval(nlsys.B,xcheb,nlsys.Parameters);
        else,
            B{i} = feval(nlsys.B,xcheb);
        end
        pwasys.Bbar{i} = [B{i}; zeros(1,size(B{i},2))];
    end
end

if size(pwafun.X,2)==1,
    for i = 1:mT,
        pwasys.L{i} = 2*[1 0]/(pwafun.X(pwafun.T(i,1))-pwafun.X(pwafun.T(i,2)));   % Luis Boyd Paper, Systems and Control Letters 2005
        pwasys.l{i} = -(pwafun.X(pwafun.T(i,1))+pwafun.X(pwafun.T(i,2)))/(pwafun.X(pwafun.T(i,1))-pwafun.X(pwafun.T(i,2)));
    end
end

% Pass local values to pwasys
pwasys.Abar = pwasys.Abar';
pwasys.Bbar = pwasys.Bbar';
pwasys.Ebar = pwasys.Ebar';
pwasys.Index = XIndex;
pwasys.X = pwafun.X;
pwasys.Y = pwafun.Y;
pwasys.T = pwafun.T;
pwasys.dt = pwafun.dt;
pwasys.W = pwafun.W;
pwasys.Z = pwafun.Z;
pwasys.Err = pwafun.Err;
pwasys.xcl = nlsys.xcl;
if isfield(nlsys,'Domain'),
    pwasys.Domain = Domain;
end
if isfield(nlsys,'NonlinearDomain'),
    pwasys.NonlinearDomain = nlsys.NonlinearDomain;
end
pwasys.istar = pwa_region(pwasys.xcl,pwasys);
