function [sys,x0,str,ts] = pwa_sfunc(t,x,u,flag,pwadsys,SimParam)

switch flag,

    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    case 0,
        [sys,x0,str,ts]=mdlInitializeSizes(pwadsys,SimParam);

        %%%%%%%%%%%%%%%
        %    Update   %
        %%%%%%%%%%%%%%%
    case 2,
        sys=mdlUpdate(t,x,u,pwadsys,SimParam);

        %%%%%%%%%%%
        % Outputs %
        %%%%%%%%%%%
    case 3,
        sys=mdlOutputs(t,x,u);

        %%%%%%%%%%%%%%%%%%%
        % Unhandled flags %
        %%%%%%%%%%%%%%%%%%%
    case { 1, 4, 9 },
        sys = [];

        %%%%%%%%%%%%%%%%%%%%
        % Unexpected flags %
        %%%%%%%%%%%%%%%%%%%%
    otherwise
        error(['Unhandled flag = ',num2str(flag)]);

end
% end pwa_dsfunc

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================

function [sys,x0,str,ts]=mdlInitializeSizes(pwadsys,SimParam)

n = length(pwadsys.Abar{1})-1;
m = size(pwadsys.Bbar{1},2);

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = n;
sizes.NumOutputs     = n;
sizes.NumInputs      = m;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = SimParam.x0;
str = [];
ts  = [pwadsys.Ts 0]; % Sample period of Ts seconds

% end mdlInitializeSizes
%
%=============================================================================
% mdlUpdate
% Return the update for the discrete states.
%=============================================================================
%
function sys=mdlUpdate(t,x,u,pwadsys,SimParam)

xbar = [x;1];

if strcmp(lower(SimParam.sys),'linear'),
    Ri = pwa_region(pwadsys.xcl,pwadsys);
    sysbar = pwadsys.Abar{Ri}*xbar+pwadsys.Bbar{Ri}*u;
    sys = sysbar(1:end-1);
elseif strcmp(lower(SimParam.sys),'pwa'),
    Ri = pwa_region(x,pwadsys);
    if isnan(Ri),
        set_param(gcs,'SimulationCommand','stop');
        sys=0*x;
    else
        sysbar = pwadsys.Abar{Ri}*xbar+pwadsys.Bbar{Ri}*u;
        sys = sysbar(1:end-1);
    end
else
    error('SimParam.sys is undefined.');
end

% end mdlUpdate
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
sys = x;
% end mdlOutputs