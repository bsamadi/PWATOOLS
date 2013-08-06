function [sys,x0,str,ts] = pwa_sfunc(t,x,u,flag,nlsys,pwasys,SimParam)
switch flag,

    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    case 0,
        [sys,x0,str,ts]=mdlInitializeSizes(pwasys,SimParam);

        %%%%%%%%%%%%%%%
        % Derivatives %
        %%%%%%%%%%%%%%%
    case 1,
        sys=mdlDerivatives(t,x,u,nlsys,pwasys,SimParam);

        %%%%%%%%%%%
        % Outputs %
        %%%%%%%%%%%
    case 3,
        sys=mdlOutputs(t,x,u);

        %%%%%%%%%%%%%%%%%%%
        % Unhandled flags %
        %%%%%%%%%%%%%%%%%%%
    case { 2, 4, 9 },
        sys = [];

        %%%%%%%%%%%%%%%%%%%%
        % Unexpected flags %
        %%%%%%%%%%%%%%%%%%%%
    otherwise
        error(['Unhandled flag = ',num2str(flag)]);

end
% end pwa_sfunc

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================

function [sys,x0,str,ts]=mdlInitializeSizes(pwasys,SimParam)

n = length(pwasys.Abar{1})-1;
m = size(pwasys.Bbar{1},2);

sizes = simsizes;
sizes.NumContStates  = n;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = n;
sizes.NumInputs      = m;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = SimParam.x0;
str = [];
ts  = [0 0];

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u,nlsys,pwasys,SimParam)

xbar = [x;1];

if strcmp(lower(SimParam.sys),'linear'),
    Ri = pwa_region(pwasys.xcl,pwasys);
    sysbar = pwasys.Abar{Ri}*xbar+pwasys.Bbar{Ri}*u;
    sys = sysbar(1:end-1);
elseif strcmp(lower(SimParam.sys),'pwa'),
    Ri = pwa_region(x,pwasys);
    if Ri==0,
        set_param(gcs,'SimulationCommand','stop');
        sys=0*x;
    else
        sysbar = pwasys.Abar{Ri,1}*xbar+pwasys.Bbar{Ri,1}*u;
        sys = sysbar(1:end-1);
    end
elseif strcmp(lower(SimParam.sys),'nonlinear'),
    if isfield(nlsys,'Parameters'),
        Param = nlsys.Parameters;             % Set the system parameters
    else
        Param = [];
    end
    fx = pwa_nlfun(nlsys,x',Param)';
    B = nlsys.B;
    try
        B = B+0*B;
    catch
        Bfun.Handle = nlsys.B;
        Bfun.xcl = nlsys.xcl;
        B = pwa_nlfun(Bfun,x',Param)';
    end
    sys = nlsys.A*x+nlsys.a+fx+B*u;
else
    error('SimParam.sys is undefined.');
end

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
sys = x;
% end mdlOutputs