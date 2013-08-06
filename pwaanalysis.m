function Y=pwaanalysis(pwasys, setting, Gain)
warning off
% Y=PWAanalysis(pwasys, xcl) analyzes a PWA model or a nonlinear model,
% approximated with PWADIs, to verify if a given point setting.xcl is a stable
% equilibrium point for the system. The reasoning serves as a sufficient
% condition for stability; The code may conclude the (nonlinear) system is
% stable, or end without any conclusion about the stability of the system.
%
% INPUT
% pwasys: PWA/PWADI model
% setting : containg xcl and other parametrs 
% Gain:  (if given) the controller
%
% OUTPUT
% prints a message for the user
% Y     : 1 if stable, 0 otherwise
%
% References:
%
% 1) L. Rodrigues and S. Boyd. Piecewise-affine state feedback for
%    piecewise-affine slab systems using convex optimization. Systems and
%    Control Letters, 54:835-853, 200
%
% 2) A. Hassibi and S. Boyd. Quadratic stabilization and control of
%    piecewise-linear systems. Proceedings of the American Control Conference,
%    6:3659-3664, 1998
%
% 3) B. Samadi and L. Rodrigues. Extension of local linear controllers to
%    global piecewise affine controllers for uncertain non-linear systems.
%    International Journal of Systems Science, 39(9):867-879, 2008.
%
%
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

%% Initialization
Y.open=0;
tol=1e-4; % tolarance for accepting the equilibrium point

[error, cont_index]=IsC1Function(pwasys);
if error==1
    clc
    fprintf('\nThe dynamics should be a C1 function.\n');
    fprintf('It has discontinuity on the border of regions R_%i and R_%i.\n\n', cont_index(1), cont_index(2));
    fprintf('I quit.\n\n');
    return;
end

setting=Initialsetting(pwasys,setting);

try
    if setting.error
        fprintf('\n\n PROGRAM STOPPED BY ERROR! EQUILIBRIUM POINT SHOULD BE DEFINED\n\n')
        setting=rmfield(setting, 'error');
        assignin('base', 'setting', setting);
        return
    end
end

xcl=setting.xcl;
pwasys.xcl=xcl;  % setting the equilibrium point

if nargin==3
    [pwasys, sysloop]=KDefine(pwasys, Gain); % ensuring K and k are defined properly
    % check to see if xcl is an equilibrium point with the given tolarance
else
    [pwasys, sysloop]=KDefine(pwasys); % ensuring K and k are defined properly
end


Isxcl=Isequilibrium(pwasys, xcl, tol);

if ~Isxcl.open & ~Isxcl.closed
    fprintf('\nxcl should be an equilibrium point for the system\n\n');
    fprintf('PWATOOL should quit.\n\n')
    return;
end


fprintf('setting parameters were set.\n')
setting
fprintf('Strike any key to continue\n')

clc

% Not set by user in workspace. These two lines should be set with respect
% to each other and also the for-loop ranging over method_index. If other
% method is added to the code, TotMeth.name, TotMeth.synth and the for-loop
% ranging over method_index should be updated accordingly.
TotMeth.name    ={'ELLIPSOIDAL    ', 'QUADRATIC CURVE', 'ELLIPSOIDAL    ', 'QUADRATIC CURVE'};
TotMeth.Lyapunov={'Global', 'Global', 'PWQ', 'PWQ'};
                  
% checking the Lyapunov function type offered by user
IsGlobal =any(strcmpi(setting.Lyapunov, 'global'));
IsPWQ =any(strcmpi(setting.Lyapunov, 'pwq'));

% checking the methods for approximating offered by user.
Isellipsoid = any(strcmpi(setting.ApxMeth, 'ellipsoidal'));
Isregion    = any(strcmpi(setting.ApxMeth, 'quadratic'));


TotMeth.index=[IsGlobal*Isellipsoid   IsGlobal*Isregion...
    IsPWQ*Isellipsoid       IsPWQ*Isregion];

[I,method_index]=find(TotMeth.index);
NumResult=length(I); % number of successful methods

%% Analysis

clc

%% Closed-loop Analysis (If applicable)
if strcmp(sysloop, 'closed') % if gains are non-zero
    [StabilityClosed, CtrlClosed]=DoAnalysis(pwasys, method_index, setting);
end

%% Open-loop Analysis
pwasys=MekeGainZero(pwasys); % Gains are zero in open-loop case. Kbar does not change
[StabilityOpen, CtrlOpen]=DoAnalysis(pwasys, method_index, setting);
pwasys=KDefine(pwasys); % Updating K and k from Kbar

%% Printing the results
clc
if strcmp(sysloop, 'closed')  % (if applicable)
    Y.closed=PrintResult(pwasys, StabilityClosed, 'closed', method_index, TotMeth);
end

Y.open=PrintResult(pwasys, StabilityOpen, 'open', method_index, TotMeth);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  MAIN CODE ENDS HERE  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Local function 1: Analysis
function  [stability, ctrl]=DoAnalysis(pwasys, method_index, setting)

% output initialization
ctrl=[];

for i=method_index
    if i==1
        % Ellipsoidal approximation: Global
        [pwasys,ctrl{i}]=AnsGlblElp(pwasys, setting);
    elseif i==2
        % Regional approximation: Global
        [pwasys,ctrl{i}]=AnsGlblReg(pwasys, setting);
    elseif i==3
        % Ellipsoidal approximation: PWQ
        [pwasys,ctrl{i}] = AnsPWQElp(pwasys, setting);
    elseif i==4
        % Regional approximation: PWQ
        [pwasys,ctrl{i}] = AnsPWQReg(pwasys, setting);
    end
    stability(i)=CheckAnalysis(ctrl{i});  % 1 is the results are trustable
end

end

%% Local function 2: Making gains zero for open-loop case
function pwasys=MekeGainZero(pwasys)

n=size(pwasys.A{1},1);
m=size(pwasys.B{1},2);
NR=size(pwasys.A,1);

for i=1:NR
    pwasys.K{i}=zeros(m,n);
    pwasys.k{i}=zeros(m,1);
end

end

%% local function 3: checking the results of analysis function

function stability=CheckAnalysis(pwactrl)
% This local function checks the results of to see if the results are
% confirming the stability
stability=0;

Problem=pwactrl.problem;
u1=pwactrl.u1;
u2=pwactrl.u2;

eligib=min(sign(u1))>=0 & norm(u2)<10;

if (eligib & (Problem==0 | Problem==4))
    stability=1; % open loop is stable with i-th method
end

end

%% Local funcrtion 4: printing
function Y=PrintResult(pwasys, stability, sysloop, method_index, TotMeth)
% isolating a linear model from a PWA model so that a definite conclusion
% can be obtained for UNSTABILITY
Y=0;

if strcmp(pwasys.type, 'pwadi')
    modeltype='nonlinear';
elseif strcmp(pwasys.type, 'lower-envelope')
    modeltype='PWA';
end

if strcmp(sysloop, 'closed')
    Looptype='closed-loop';
elseif strcmp(sysloop, 'open')
    Looptype='open-loop';
end

ONEregion=strcmp(pwasys.type, 'lower-envelope') & pwasys.NR==1;
ONEregion=double(ONEregion);

fprintf('\n\n')

if any(stability)
    Y=1;
    fprintf('The following methods verified that the %s %s system is stable at xcl. \n\n', Looptype, modeltype);
    counter=0;
    for i=method_index
        if stability(i)
            counter=counter+1;
            fprintf('(%i)  %s approximation (%s Lyapunov function)\n', counter, TotMeth.name{i}, TotMeth.Lyapunov{i});
        end
    end
else
    if ONEregion
        fprintf('The %s %s system is NOT stable at xcl. \n', Looptype, modeltype);
    else
        fprintf('I could not verify if the %s %s system is stable at xcl. \n', Looptype, modeltype);
    end
end

fprintf('\n\n')
end



%% Local function 1: initializing the setting

function setting=Initialsetting(pwasys,setting)
% This function initializes setting

n = size(pwasys.Abar{1},1)-1;   % order of the system
m = size(pwasys.Bbar{1},2);     % Number of the inputs
clc

if ~isfield(setting, 'Lyapunov')
    setting.Lyapunov= {'global', 'pwq'};
end

if ~isfield(setting, 'ApxMeth')
    setting.ApxMeth={'ellipsoidal', 'quadratic'};
end

if ~isfield(setting, 'alpha')
    setting.alpha=.1;
end

if strcmp(pwasys.type, 'lower-envelope')
    setting.systype='pwa';
    assignin('base', 'nlsys', []);
elseif strcmp(pwasys.type, 'pwadi')
    setting.systype='nonlinear';
end

if isfield(setting, 'xcl')
    if length(setting.xcl)~=n
        fprintf('\n Dimension of setting.xcl does not match the order of the system.\n\n')
        setting.error=1;
    end
else
    if isfield(pwasys, 'xcl')
        setting.xcl=pwasys.xcl;
    else
        fprintf('\n setting.xcl is missing. Equilibrium point should be defined.\n\n')
        setting.error=1;
    end
end


end


