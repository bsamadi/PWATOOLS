function pwactrl = pwa_ctrlsynsos(pwasys, pwactrl,Param)

% State Feedback Controller Synthesis
% pwactrl = pwa_ctrlsynsf(pwasys,Param)

% pwasys
% xbarDot = Abar{i}*xbar + Bbar{i}*u
%
% pwasys.Abar = {Abar1, ..., Abarn}
% pwasys.Bbar = {Bbar1, ..., Bbarn}
% pwasys.Cybar = {Cybar1, ..., Cybarn}
% pwasys.Czbar = {Czbar1, ..., Czbarn}
% pwasys.Dzw = {Dzw1, ..., Dzwn}
% pwasys.Dzu = {Dzu1, ..., Dzun}
%
% Polytopic region i: {x| Ebar{i}*xbar > 0 }
%
% Boundary information table:
% Boundary between Ri and Rj : {x| xbar=Fbar{i,j}*sbar }
%
% Param.Lyapunov = 'SOS'
% Param.Order = Order of SOS function (should be even)
% Param.LimitC = Limit for polynomial coefficients

% pwasys.n number of state variables
% pwasys.m number of control inputs
% pwasys.p number of constraints for each region
% pwasys.q number of disturbance inputs
% pwasys.r number of control feedback outputs

if nargin==1,
    pwactrl = [];
end
pwasys = center_region(pwasys);

for i = pwasys.istar,
    A{i} = pwasys.Abar{i}(1:pwasys.n,1:pwasys.n);
    B{i} = pwasys.Bbar{i}(1:pwasys.n,:);
end
pwactrl.Lin2PWA = 1;
if isfield(Param,'Poles'),
    for i = pwasys.istar,
        K{i} = -acker(A{i},B{i},Param.Poles);
    end
elseif isfield(Param,'KLin'),
    for i = pwasys.istar,
        K{i} = Param.KLin;
    end
elseif isfield(Param,'QLin'),
    for i = pwasys.istar,
        K{i} = -lqr(A{i},B{i},Param.QLin,Param.RLin);
    end
else
    pwactrl.Lin2PWA = 0;
end


if ~isfield(Param,'Cont'),
    Param.Cont=1;
end

if ~isfield(Param,'alpha'),
    Param.alpha=0;
end

if ~isfield(Param,'Order'),
    Param.Order=2;
end

x = sdpvar(pwasys.n,1);
xbar = [x;1];

F = set([]);
% Polynomial Lyapunov
if ~isfield(pwactrl,'Pbar'),
    xm = (monolist(x,Param.Order/2));
    P = sdpvar(length(xm));
    V = xm'*P*xm;
    F = set(P(1,1)==0);

    if isfield(Param,'LimitC'),
        F = F+set(-Param.LimitC<P(:)<Param.LimitC);
    else
        F = F+set(-Param.LimitP<P(:)<Param.LimitP);
    end

    % Positive definite
    eps = sdpvar(Param.Order/2,1);
    F = F + set(eps(:)>0);
    LB = 0;
    for i=1:Param.Order/2,
        LB = LB+eps(i)*(x'*x)^i;
    end
    F = F + set(sos(V-LB));    % Lower Bound
else,
    try,
        V = xbar'*pwactrl.Pbar*xbar;
    catch,
        V = xbar'*pwactrl.Pbar{1}*xbar;
    end
end

if isfield(Param,'L2Gain'),
    L2Gain = Param.L2Gain;
    gammaSQ = sdpvar(1);
    F = F + set(gammaSQ > Param.L2Gain^2,['gamma^2>' num2str(Param.L2Gain^2)]);
    for i=1:pwasys.NR,
        if ~iscell(pwasys.Bwbar),
            Bwbar{i}=pwasys.Bwbar;
        end
        if ~iscell(pwasys.Czbar),
            Czbar{i}=pwasys.Czbar;
        end
        if isfield(pwasys,'Cybar'),
            if ~iscell(pwasys.Cybar),
                Cybar{i}=pwasys.Cybar;
            end
        end
        if ~isfield(pwasys,'Dzu'),
            for i=1:pwasys.NR,
                Dzu{i} = zeros(size(Czbar{i},1),pwasys.m);
            end
        end
        pwasys.q = size(Bwbar{1},2);    % Number of disturbance inputs
        if ~isfield(pwasys,'Dzw'),
            for i=1:pwasys.NR,
                Dzu{i} = zeros(size(Czbar{i},1),pwasys.q);
            end
        end
        if ~iscell(pwasys.Dzw),
            Dzw{i}=pwasys.Dzw;
        end
        if ~iscell(pwasys.Dzu),
            Dzu{i}=pwasys.Dzu;
        end
    end
    for i=1:pwasys.NR,
        Bw{i} = Bwbar{i}(1:end-1,:);
        Cz{i} = Czbar{i}(:,1:end-1);
        if isfield(pwasys,'Cybar'),     % Static output feedback
            Cy{i} = Cybar{i}(:,1:end-1);
        end
    end
    if ~iscell(pwasys.Czbar),
        pwasys.Czbar = Czbar;
        if isfield(pwasys,'Cybar'),     % Static output feedback
            pwasys.Cybar = Cybar;
            pwasys.r = size(Cybar{i},1);
        end
        pwasys.Bwbar = Bwbar;
        pwasys.Dzw = Dzw;
        pwasys.Dzu = Dzu;
    end
end

% Controller
Options = [];
if isfield(pwactrl,'Analysis'),
    if pwactrl.Analysis,
        Kbar = pwactrl.Kbar;
    end
elseif pwactrl.Lin2PWA,
    for i=1:pwasys.NR,
        if ismember(i,pwasys.istar),      % For the center region(s)
            k{i} = sdpvar(pwasys.m,1);
            F = F + set(-Param.LimitK<k{i}(:)<Param.LimitK);
            Kbar{i} = [K{i} k{i}];
        else
            if isfield(pwactrl,'Kbar'),
                Kbar{i} = pwactrl.Kbar{i};
            else,
                if isfield(pwasys,'Cybar'),     % Static output feedback
                    Kbar{i} = sdpvar(pwasys.m,pwasys.r+1,'full');
                else,
                    Kbar{i} = sdpvar(pwasys.m,pwasys.n+1,'full');
                end
                F = F + set(-Param.LimitK<Kbar{i}(:)<Param.LimitK);
            end
        end
    end
else
    for i=1:pwasys.NR,
        if isfield(pwasys,'Cybar'),     % Static output feedback
            Kbar{i} = sdpvar(pwasys.m,pwasys.r+1,'full');
        else,
            Kbar{i} = sdpvar(pwasys.m,pwasys.n+1,'full');
        end
        F = F + set(-Param.LimitK<Kbar{i}(:)<Param.LimitK);
        if isfield(pwactrl,'Kbar'),
            assign(Kbar{i},pwactrl.Kbar{i});
            Options = sdpsettings('usex0',1);
        end
    end
end

% Desired equilibrium point
for i=pwasys.istar,
    try,
        if isfield(pwasys,'Cybar'),     % Static output feedback        
            F = F + set( pwasys.Abar{i}*[pwasys.xcl; 1]+pwasys.Bbar{i}*Kbar{i}*[pwasys.Cybar{i}*[pwasys.xcl;1]; 1] == 0,['(Abar+Bbar*Kbar)*xclbar=0']);
        else
            F = F + set( (pwasys.Abar{i}+pwasys.Bbar{i}*Kbar{i})*[pwasys.xcl; 1] == 0,['(Abar+Bbar*Kbar)*xclbar=0']);
        end
    catch,
        disp('Check the constraint for the desired equilibrium point.');
    end
end

if Param.Cont,
    for i=1:pwasys.NR,
        for j = i:pwasys.NR,
            if ~isempty(pwasys.Fbar{i,j}),
                % Continuity of the control input
                try,
                    if isfield(pwasys,'Cybar'),     % Static output feedback
                        F = F + set(Kbar{i}*[pwasys.Cybar{i}*pwasys.Fbar{i,j};zeros(1,pwasys.n-1) 1]== Kbar{j}*[pwasys.Cybar{j}*pwasys.Fbar{i,j};zeros(1,pwasys.n-1) 1],['Ki*F' num2str(i) '_' num2str(j) 's=Kis*F' num2str(i) '_' num2str(j) 's']);
                    else,
                        F = F + set((Kbar{i}-Kbar{j})*pwasys.Fbar{i,j}== 0,['Ki*F' num2str(i) '_' num2str(j) 's=Kis*F' num2str(i) '_' num2str(j) 's']);
                    end
                catch,
                    disp('Check the continuity of the controller.');
                end
            end
        end
    end
end

if exist('L2Gain'),
    w = sdpvar(pwasys.q,1);
end
for i = 1:pwasys.NR,
    % Description of the regions
    S = pwasys.Ebar{i}*xbar;
    S = S(1:end-1);

    % S procedure parameters
        [L(i),cs{i}] = polynomial(S,Param.Order-1,1);
        F = F + set(cs{i}(:)>0);
    if exist('L2Gain'),
        if isfield(pwasys,'Cybar'),     % Static output feedback
            xDot = pwasys.Abar{i}*xbar+pwasys.Bbar{i}*Kbar{i}*[pwasys.Cybar{i}*xbar; 1]+pwasys.Bwbar{i}*w;            
            z = pwasys.Czbar{i}*xbar+pwasys.Dzu{i}*Kbar{i}*[pwasys.Cybar{i}*xbar; 1]+pwasys.Dzw{i}*w;
        else,                           % State feedback
            xDot = (pwasys.Abar{i}+pwasys.Bbar{i}*Kbar{i})*xbar+pwasys.Bwbar{i}*w;
            z = pwasys.Czbar{i}*xbar+pwasys.Dzu{i}*Kbar{i}*xbar+pwasys.Dzw{i}*w;
        end
        xDot = xDot(1:pwasys.n);
        F = F + set(sos(-jacobian(V,x)*xDot-z'*z+gammaSQ*w'*w-L(i)));
    else
        % Stability constraint dV/dt < -alphaV
        xDot = (pwasys.Abar{i}+pwasys.Bbar{i}*Kbar{i})*xbar;
        F = F + set(sos(-jacobian(V,x)*xDot(1:end-1)-Param.alpha*V-L(i)));
    end
end

pwactrl.Kbar = Kbar;

Obj = [];
if isfield(Param,'L2Gain'),
    Obj = gammaSQ;
end

[sol,v,Q] = solvesos(F,Obj,[Options,sdpsettings('sos.model',2,'sos.traceobj',1)]);
checkset(F)
sol

if exist('P'),
    P = double(P);
    P = [P(2:end,2:end) P(2:end,1); P(1,2:end) P(1,1)];
    for i=1:pwasys.NR,
        pwactrl.Pbar{i} = P;
    end
    % pwactrl.xm = sdisplay(xm);
    pwactrl.Note = 'Pbar was rearrenged inside the SOS design process to make it compatible with the xbar notation.';
end


for i=1:pwasys.NR,
    pwactrl.Kbar{i}=double(pwactrl.Kbar{i});
end
if isfield(Param,'L2Gain'),
    pwactrl.gamma = sqrt(double(gammaSQ));
end

pwactrl.Index = pwasys.Index;
pwactrl.X = pwasys.X;
pwactrl.T = pwasys.T;
pwactrl.xcl = pwasys.xcl;
pwactrl.Ebar = pwasys.Ebar;
