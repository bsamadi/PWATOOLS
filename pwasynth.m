function [Control, DataSim]=pwasynth(pwasys, x0, setting)
warning off
% Control=pwasynth(pwasys, x0, option) synthesizes PWA controllers for pwasys
% and simulates the closed-loop system, starting from x0. option  contains the
% parameters that affect control design and simulation
%
% INPUT
% pwasys  : PWA differential inclusions or PWA model
% x0      : Initial point
% setting : contians the following fields
% |------------------------------------------------------------------------
% |FIELD            STRUCTURE       SIZE             DEFAULT VALUE
% |------------------------------------------------------------------------
% |ApxMeth          Cell           (1 x 2)           {'ellipsoidal', 'quadratic'}
% |Lyapunov         Cell           (1 x 2)           {'global', 'pwq'}
% |SynthMeth        Cell           (1 x 2)           {'lmi', 'bmi'}
% |QLin             Array          (n x n)           Random QLin>0
% |RLin             Array          (m x m)           Random RLin>0
% |RandomQ          scalar         (1 x 1)           1
% |RandomR          scalar         (1 x 1)           1
% |alpha            scalar         (1 x 1)           0.1
% |xcl              Array          (n x 1)           pwasys.xcl (if exists)
% |StopTime         scalar         (1 x 1)           10 sec
% |IterationNumber  scalar         (1 x 1)           5
% |------------------------------------------------------------------------
%
%
% OUTPUT
% Control    : controller gains
% pwacont    : controller gains
% pwasetting : final setting used for design
% DataSim    : Simulation data from simulink
%
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%
%% Is pwasys a C1 function?
[error, cont_index]=IsC1Function(pwasys);
if error==1
    clc
    fprintf('\nThe dynamics should be a C1 function.\n');
    fprintf('It has discontinuity on the border of regions R_%i and R_%i.\n\n', cont_index(1), cont_index(2));
    fprintf('I quit.\n\n');
    Control=[];
    DataSim=[];
    return;
end

%% Initialization
setting=Initialsetting(pwasys, setting);

try
    if setting.error
        fprintf('\n\n PROGRAM STOPPED BY ERROR!\n\n')
        setting=rmfield(setting, 'error');
        assignin('base', 'setting', setting);
        return
    end
end

fprintf('setting parameters were set.\n')
setting
fprintf('Strike any key to continue\n')

if ~setting.failure
    pause;
end


clc
pwactrl = [];
xcl=setting.xcl;
pwasys.x0=x0;

% Not set by user in workspace. These two lines should be set with respect
% to each other and also the for-loop ranging over method_index. If other
% method is added to the code, TotMeth.name, TotMeth.synth and the for-loop
% ranging over method_index should be updated accordingly.
TotMeth.name={'ELLIPSOIDAL', 'QUADRATIC CURVE', 'SLAB REGION', 'ELLIPSOIDAL', 'QUADRATIC CURVE'};
TotMeth.synth={'BMI', 'BMI', 'LMI', 'BMI', 'BMI'};
TotMeth.Lyapunov={'GLOBAL', 'GLOBAL', 'GLOBAL', 'PWQ', 'PWQ'};

% checking the Lyapunov function type offered by user
IsGlobal =any(strcmpi(setting.Lyapunov, 'global'));
IsPWQ =any(strcmpi(setting.Lyapunov, 'pwq'));

% checking the methods for synthesizeing offered by user
IsLMI =any(strcmpi(setting.SynthMeth, 'lmi'));
IsBMI =any(strcmpi(setting.SynthMeth, 'bmi'));

% checking the methods for approximating offered by user.
Isellipsoid = any(strcmpi(setting.ApxMeth, 'ellipsoidal'));
Isregion    = any(strcmpi(setting.ApxMeth, 'quadratic'));


TotMeth.index=[IsGlobal*IsBMI*Isellipsoid   IsGlobal*IsBMI*Isregion...
               IsLMI IsPWQ*IsBMI*Isellipsoid  IsPWQ*IsBMI*Isregion]; 

[I,method_index]=find(TotMeth.index);
NumResult=length(I); % number of successful methods

%% Synthesis
for i=method_index
    if i==1
        % Ellipsoidal approximation: Global
        [pwasys,ctrl{i}] = SyntGlblBMIElp(pwasys, setting);
    elseif i==2
        % Regional approximation: Global
        [pwasys,ctrl{i}] = SyntGlblBMIReg(pwasys, setting);
    elseif i==3
        % Slab Regions: Global
        [pwasys,ctrl{i}] = SyntGlblLMI(pwasys, setting);
    elseif i==4
        % Ellipsoidal approximation: PWQ
        [pwasys,ctrl{i}] = SyntPWQBMIElp(pwasys, setting);
    elseif i==5
        % Regional approximation: PWQ
        [pwasys,ctrl{i}] = SyntPWQBMIReg(pwasys, setting);
    end
    CheckResult(i)=CheckSynthesis(ctrl{i});  % 1 is the results are trustable
end

clc

if ~setting.failure  % initially without any error occurred, number of failure is set to 0
    setting.numberoffailure=0;
end

if any(CheckResult)
    fprintf('\n');
    fprintf('PWATOOL was successful in synthesizing PWA controllers by the following method(s)\n\n');
    counter=0;
    for i=method_index
        if  CheckResult(i)
            counter=counter+1;
            fprintf('(%i)  %s approach with %s approximation (%s Lyapunov)\n',counter, TotMeth.synth{i}, TotMeth.name{i}, TotMeth.Lyapunov{i});
        end
    end
       
    fprintf('\n');
    fprintf('\n\nWaiting for Simulink model to run....\n\n');
    % printing the results and simulation
    [Control, DataSim]=pwasimulation(pwasys, ctrl, CheckResult, TotMeth, method_index, setting); 
else % if no method has convergeed we may try to the maximum number of IterationNumber
    setting.numberoffailure=setting.numberoffailure+1;
    setting.failure=1;
    if setting.numberoffailure<setting.IterationNumber
        if setting.numberoffailure==1
            fprintf('\n\nPWATOOL was unsuccessful in synthesizing PWA controllers in its try\n\n');
            fprintf('PWATOOL will try at most %i more times\n\n', setting.IterationNumber-1);
            setting.confirm=input('Do you want to continue? (Y/N): ', 's');
            if strcmpi('y', setting.confirm)
                setting.RandomQBck=setting.RandomQ; % back-up of the original Q
                setting.RandomRBck=setting.RandomR; % back-up of tghe original R
                setting.RandomQ=1;  % make sure Q and R are randomely set so that
                setting.RandomR=1;  % we can change the initialization
            end
        end
        if strcmpi('y', setting.confirm)
            Control=pwasynth(pwasys, x0, setting);
        else
            Control=[];
            DataSim=[];
            fprintf('\n\n');  % exits automatically in this case
        end
    else
        fprintf('\n\nPWATOOL could not find any PWA controller for the system after\n');
        fprintf('%i number(s) of iteration\n\n', setting.IterationNumber);
        Control=[];
        DataSim=[];
        return;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  MAIN CODE ENDS HERE  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Local function 1: initializing the setting

function setting=Initialsetting(pwasys,setting)
% This function initializes setting

n = size(pwasys.Abar{1},1)-1;   % order of the system
m = size(pwasys.Bbar{1},2);     % Number of the inputs
clc
if ~isfield(setting, 'SynthMeth')
    setting.SynthMeth= {'LMI', 'BMI'};
end

if ~isfield(setting, 'Lyapunov')
    setting.Lyapunov= {'global', 'pwq'};
end

if ~isfield(setting, 'ApxMeth')
    setting.ApxMeth={'ellipsoidal', 'quadratic'};
end

if ~isfield(setting, 'alpha')
    setting.alpha=.1;
end

if ~isfield(setting, 'NormalDirectionOnly')
    setting.NormalDirectionOnly=0;
end


if ~isfield(setting, 'StopTime')
    setting.StopTime=10;
end

if ~isfield(setting, 'IterationNumber')
    setting.IterationNumber=5;
end

if ~isfield(setting, 'RandomQ') || isempty(setting.QLin)
    setting.RandomQ=1;
end

if ~isfield(setting, 'RandomR') || isempty(setting.RLin)
    setting.RandomR=1;
end

if setting.RandomQ
    Q=rand(n); Q=Q'*Q*10^3;
    setting.QLin=Q;
end

if setting.RandomR
    R=rand(m); R=3*R'*R*10^3;
    setting.RLin=R;
end


if strcmp(pwasys.type, 'lower-envelope')
    setting.systype='pwa';
elseif strcmp(pwasys.type, 'pwadi')
    setting.systype='nonlinear';
end

setting.ctrltype='pwa';


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


    
    
if ~isfield(setting, 'failure')
    setting.failure=0;
end


end


%% local function 2: checking the results for stability and reliability

function stability=CheckSynthesis(pwactrl)
% This local function checks the results of both ellipsoidal and regional
% methods on synthesis to see if any of the results are correct

stability=0; % no conclusion

u1=pwactrl.u1;
u2=pwactrl.u2;

eligib=min(sign(u1))>=0 & norm(u2)<10;
problem=pwactrl.problem;

if (eligib & problem==0)
    stability=1; % closed loop is stable with i-th method
end

end


%% local function 3: printing the result in case of success
function [Control, DataSim]=pwasimulation(pwasys, ctrl, CheckResult, TotMeth, method_index, setting)
close all
x0=pwasys.x0;
j=0;
for i=method_index
    if CheckResult(i)
        j=j+1;
        pwactrl=ctrl{i};
        Control.Gain{j}=pwactrl.Kbar;
        Control.AppMethod{j}=TotMeth.name{i};
        pwactrl.AppMethod=TotMeth.name{i};
        Control.SynthesisType{j}=TotMeth.synth{i};
        Control.Lyapunov{j}=TotMeth.Lyapunov{i};
        try
            Control.Sprocedure1{j}=pwactrl.Z;
        end
        try
            Control.Sprocedure1{j}=pwactrl.miu;
        end
        try
            Control.Sprocedure2{j}=pwactrl.W;
        end
        try
            Control.Sprocedure2{j}=pwactrl.bita;
        end
        % assiging the results in BASE workspace for future usage of Simulink
        assignin('base', 'setting', setting);
        assignin('base', 'xcl', setting.xcl);
        evalin('base', 'SimParam.StopTime=setting.StopTime;');
        evalin('base', 'SimParam.sys  = setting.systype;');
        evalin('base', 'SimParam.ctrl = setting.ctrltype;');
        assignin('base', 'x0', x0);
        evalin('base', 'SimParam.x0=x0;');
        assignin('base', 'pwactrl', pwactrl);
        assignin('base', 'pwaModel', pwasys);
        %% in case of convergence run the simulation
        if strcmp(pwasys.type, 'lower-envelope')
            assignin('base', 'NonSys', []);
        elseif strcmp(pwasys.type, 'pwadi')
            evalin('base', 'NonSys=nlsys;');
        end
        
        sim('pwa_sim');
        n = size(pwasys.Abar{1},1)-1;   % order of the system
        m = size(pwasys.Bbar{1},2);     % Number of the inputs
 
        title_text=strcat(TotMeth.synth{i}, '--', TotMeth.name{i}, '--', TotMeth.Lyapunov{i});
        figure(j);
        time=sim_ctrl(:,1);
        subplot(211)
        plot(time, sim_states(:,2:n+1));
        title(strcat(title_text, ' Method: Error Signals'));
        grid
        
        subplot(212)
        plot(time, sim_ctrl(:,2:m+1));
        title(strcat(title_text, ' Method: Controller Output Signals'));
        grid
        figure(j)
        DataSim.states{j}= sim_states;
        DataSim.ctrl{j}= sim_ctrl;
    end
end
assignin('base', 'DataSim', DataSim);


numberoffigures=j;
fprintf('\nStates and control signals transitions are shown in figure(s) ');
for i=1:numberoffigures-1
    fprintf('%i, ', i);
end
fprintf('%i. \n\n', numberoffigures);
fprintf('Controller gains are stored in\n');
assignin('base', 'pwacont', Control);
pwacont=Control;
pwacont

fprintf('\n')

% order of these operation should be kept VVVVVVVVVVVVVVVV
if isfield(setting, 'confirm')
    setting=rmfield(setting, {'confirm'});
    pwasetting=setting;
    setting.RandomQ=setting.RandomQBck;
    setting.RandomR=setting.RandomRBck;
    setting=rmfield(setting, {'RandomQBck', 'RandomRBck'});
    pwasetting=rmfield(pwasetting, {'RandomQBck', 'RandomRBck'});
else
    pwasetting=setting;
end
% order of these operation should be kept ^^^^^^^^^^^^^^^^^

assignin('base', 'pwasetting', pwasetting); % assigning pwasetting
setting=rmfield(setting, {'failure', 'numberoffailure'});
assignin('base', 'setting', setting); % updating setting in the BASE workspace
fprintf('Settings used in synthesizing the PWA controller are stored in\n')
pwasetting


end