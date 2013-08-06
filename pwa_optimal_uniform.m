function pwa_fun = pwa_optimal_uniform(nlfun)

% PWA_Optimal_Uniform Piecewise-affine approximation.
%
% pwafun = pwa_optimal_uniform(nlfun)
%
% nlfun.Handle = @function_name;
% nlfun.Domain = [xmin xmax];
% nlfun.UGR = Uniform grid resolution;
% nlfun.Parameters = Parameter vector; % If it is needed
% nlfun.Resolution = Fine grid resolution;
% nlfun.ObjFun = 'L2'; % or 'Linf'
% nlfun.AP = Anchoring point; % Optional
% nlfun.AbarLin = Linear approximation of the function at the anchoring point;  % Optional
%
% pwafun.X = Griding points
% pwafun.Y = The value of the pwafun at griding points
% pwafun.T = Triangles
% pwafun.Abar = Linear coefficient for each triangle; % Y = Abar{i}*X
% pwafun.Err = Approximation error for the fine mesh
% pwafun.W = Fine mesh points
% pwafun.Z = The value of pwafun at W
% pwafun.Obj = The value of the objective function
%
% Copyright: Behzad Samadi, Concordia University August 2005

% Argument checking

f = nlfun.Handle;                     % Specify the handle of the nonlinear function


a=[];
for i=1:length(nlfun.Domain),
    a=[a; nlfun.Domain{i}(1)];    % Set start point of the domain
end

if isfield(nlfun,'Parameters'),
    Param = nlfun.Parameters;             % Set the system parameters
else
    Param = [];
end

if isfield(nlfun,'ObjFun'),
    Option.ObjFun = nlfun.ObjFun;                   % Specify the objective: L2 or Linf
end

if isfield(nlfun,'xcl') && isfield(nlfun,'AbarLin'),
    Option.xcl = nlfun.xcl;              % Specify the anchoring point
    Option.AbarLin = nlfun.AbarLin;      % Specify the linear approximation at the anchoring point
end

n=length(a);

try
    Y = pwa_nlfun(nlfun,a,Param);
catch
    error(['The function can not be evaluated.']);
end

% Construct the fine mesh
%
% Z = f(W)
% W = [w1; w2; ... ; wm] is m times n
% Z = [z1; z2; ... ; zm] is m times p
%
% zi = f(wi)

% fine mesh
for i=1:length(nlfun.Domain),
    Domain{i} = nlfun.Domain{i}([1 end]);
end
W = pwa_grid(Domain,nlfun.Resolution);
Z = pwa_nlfun(nlfun,W,Param);

% The Approximation
X = pwa_grid(nlfun.Domain,nlfun.UGR);
if n == 1,
    X = sort(X);
end

if  ~exist('Option'),
    Option = [];
end

[Y Err Obj T Abar] = pwa_estim(W,Z,X,Option);

pwa_fun.X = X;
pwa_fun.Y = Y;
pwa_fun.T = T;
pwa_fun.Err = Err;
pwa_fun.Obj = Obj;
pwa_fun.W = W;
pwa_fun.Z = Z;
pwa_fun.Abar = Abar;

% pwa_optimal_uniform