function pwactrl = pwa_ctrlsynsf(pwasys, pwactrl, Param)

% State Feedback Controller Synthesis
% pwactrl = pwa_ctrlsynsf(pwasys,pwactrl,Param)

% pwasys
% xbarDot = Abar{i}*xbar + Bbar{i}*u
%
% pwasys.Abar = {Abar1, ..., Abarn}
% pwasys.Bbar = {Bbar1, ..., Bbarn}
%
% Polytopic region i: {x| Ebar{i}*xbar > 0 }
%
% Boundary information table:
% Boundary between Ri and Rj : {x| xbar=Fbar{i,j}*sbar }
%
% Param.Lyapunov = 'Global', 'PWQ'

pwasys = center_region(pwasys);
yalmip('clear');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                pwactrl = Lyapunov_function(pwasys, pwactrl, Param);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pwactrl.constraints = set([]);

if nargin==1,
    pwactrl = [];
end

I = eye(pwasys.n);
xcl = pwasys.xcl;

if isfield(Param,'Lyapunov'),
    Lyapunov = Param.Lyapunov;
else
    Lyapunov = 'PWQ';
end

if ~isfield(pwactrl,'Pbar'),
    % Lyapunov function
    if strcmp(Lyapunov,'Global'),
        Pg = sdpvar(pwasys.n,pwasys.n);
        pwactrl.constraints = pwactrl.constraints+set(Pg > 0*I,['P>0']);
        for i=1:pwasys.NR,
            P{i} = Pg;
            % Constrainted Pbar
            Pbar{i} = [P{i} -P{i}*xcl; -xcl'*P{i} xcl'*P{i}*xcl];
        end
    else,
        for i=1:pwasys.NR,
            if ismember(i,pwasys.istar),
                P{i} = sdpvar(pwasys.n,pwasys.n);
                pwactrl.constraints = pwactrl.constraints+set(P{i} > 0.01*I,['P' num2str(i) '>0.01']);
                if isfield(Param,'LimitP'),
                    pwactrl.constraints=pwactrl.constraints+set(P{i}(:) < Param.LimitP,['P' num2str(i) '<' num2str(Param.LimitP) '*I']);
                end
                Pbar{i} = [P{i} -P{i}*xcl; -xcl'*P{i} xcl'*P{i}*xcl];
            else
                Pbar{i} = sdpvar(pwasys.n+1,pwasys.n+1);
                Z{i} = sdpvar(pwasys.p+1,pwasys.p+1);
                % S procedure
                pwactrl.constraints=pwactrl.constraints+set(0 < Z{i}(:) < 1e5,['Z{' num2str(i) '}(:)>0']);
                pwactrl.constraints=pwactrl.constraints+set(Pbar{i}-pwasys.Ebar{i}'*Z{i}*pwasys.Ebar{i} > 0,['Pbar' num2str(i) '>0']);
                if isfield(Param,'LimitP'),
                    pwactrl.constraints=pwactrl.constraints+set(Pbar{i}(:) < Param.LimitP,['Pbar' num2str(i) '<' num2str(Param.LimitP) '*I']);
                end
            end
        end
    end
    pwactrl.Pbar = Pbar;
end

if strcmp(Lyapunov,'PWQ'),
    for i=1:pwasys.NR,
        for j = i:pwasys.NR,
            if ~isempty(pwasys.Fbar{i,j}),
                % Continuity of the Lyapunov function
                try,
                    pwactrl.constraints=pwactrl.constraints+set(pwasys.Fbar{i,j}'*(Pbar{i}-Pbar{j})*pwasys.Fbar{i,j} == 0,['F' num2str(i) '_' num2str(j) 's*Pbari*F' num2str(i) '_' num2str(j) 's=F' num2str(i) '_' num2str(j) 's*Pbaris*F' num2str(i) '_' num2str(j) 's']);
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 pwactrl = pwa_controller(pwasys, pwactrl, Param);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pwactrl.Lin2PWA = 1;
for i = pwasys.istar,
    A{i} = pwasys.Abar{i}(1:pwasys.n,1:pwasys.n);
    B{i} = pwasys.Bbar{i}(1:pwasys.n,:);
    if isfield(Param,'Poles'),
        K{i} = -acker(A{i},B{i},Param.Poles);
    elseif isfield(Param,'KLin'),
        K{i} = Param.KLin;
    elseif isfield(Param,'QLin'),
        K{i} = -lqr(A{i},B{i},Param.QLin,Param.RLin);
    else
        pwactrl.Lin2PWA = 0;
    end
end

if ~isfield(Param,'alpha'),
    Param.alpha=0;
end

if isfield(pwactrl,'Analysis'),
    if pwactrl.Analysis,
        Kbar = pwactrl.Kbar;
    end
elseif ~isfield(pwactrl,'Kbar'),
    for i=1:pwasys.NR,
        Kbar{i} = sdpvar(pwasys.m,pwasys.n+1,'full');
    end
else,
    if pwactrl.Lin2PWA,
        Kbar = pwactrl.Kbar;
    else,
        for i=1:pwasys.NR,
            Kbar{i} = sdpvar(pwasys.m,pwasys.n+1,'full');
            assign(Kbar{i},pwactrl.Kbar{i})
        end
    end
end

if isfield(Param,'L2Gain'),
    L2Gain = Param.L2Gain;
    gammaSQ = sdpvar(1);
    pwactrl.constraints=pwactrl.constraints+set(gammaSQ > Param.L2Gain^2,['gamma^2>' num2str(Param.L2Gain^2)]);
    for i=1:pwasys.NR,
        if ~iscell(pwasys.Czbar),
            Czbar{i}=pwasys.Czbar;
        end
        if ~iscell(pwasys.Bwbar),
            Bwbar{i}=pwasys.Bwbar;
        end
        if ~iscell(pwasys.Dzw),
            Dzw{i}=pwasys.Dzw;
        end
    end
    if ~isfield(pwasys,'Du'),
        for i=1:pwasys.NR,
            Dzu{i} = zeros(size(Czbar{i},1),pwasys.m);
        end
    end
    for i=1:pwasys.NR,
        Bw{i} = Bwbar{i}(1:end-1,:);
        Cz{i} = Czbar{i}(:,1:end-1);
    end
    if ~iscell(pwasys.Czbar),
        pwasys.Czbar = Czbar;
        pwasys.Bwbar = Bwbar;
        pwasys.Dzw = Dzw;
        pwasys.Dzu = Dzu;
    end
    q = size(Bwbar{1},2);    % Number of disturbance inputs
end

%%%%%%%%%%%%%%%%%%%%%%%% The center region %%%%%%%%%%%%%%%%%%%
I = eye(pwasys.n);
One = zeros(pwasys.n+1);
One(end,end)=1;

if pwactrl.Lin2PWA,
    for i=pwasys.istar,
        k = sdpvar(pwasys.m,1);
        Kbar{i} = [K{i} k];
    end
end

for i=pwasys.istar,
    K{i} = Kbar{i}(:,1:pwasys.n);
    for j=1:pwasys.NS,
        try
             pwactrl.constraints=pwactrl.constraints+set( (pwasys.Abar{i,j}+pwasys.Bbar{i,j}*Kbar{i})*[pwasys.xcl;1] == 0,'(Abar+Bbar*Kbar)*xclbar=0');
        end
        if exist('L2Gain'),
            L2Mat11 = (pwasys.Abar{i}+pwasys.Bbar{i}*Kbar{i})'*pwactrl.Pbar{i}+0.5*(pwasys.Czbar{i}+pwasys.Dzu{i}*Kbar{i})'*(pwasys.Czbar{i}+pwasys.Dzu{i}*Kbar{i});
            L2Mat11 = L2Mat11(1:pwasys.n,1:pwasys.n);
            L2Mat11 = L2Mat11+L2Mat11';
            L2Mat12 = pwactrl.Pbar{i}*pwasys.Bwbar{i}+(pwasys.Czbar{i}+pwasys.Dzu{i}*Kbar{i})'*pwasys.Dzw{i};
            L2Mat22 = pwasys.Dzw{i}'*pwasys.Dzw{i}-gammaSQ*eye(q);
            pwactrl.constraints=pwactrl.constraints+set( [L2Mat11 L2Mat12(1:pwasys.n,:); L2Mat12(1:pwasys.n,:)' L2Mat22] < 0,['L2Mat' num2str(i) '<0']);
%             pwactrl.constraints=pwactrl.constraints+set([(A{i}+B{i}*K{i})'*P{i}+P{i}*(A{i}+B{i}*K{i})+(Cz{i}+Dzu{i}*K{i})'*(Cz{i}+Dzu{i}*K{i}) P{i}*Bw{i}+(Cz{i}+Dzu{i}*K{i})'*Dzw{i}; (P{i}*Bw{i}+(Cz{i}+Dzu{i}*K{i})'*Dzw{i})' Dzw{i}'*Dzw{i}-gammaSQ*eye(size(Dzw{i},2))]<0,'L2Gain<gamma');
        else,
            if pwactrl.Lin2PWA,
                DV{i,j} = pwactrl.Pbar{i}*(pwasys.Abar{i,j}+pwasys.Bbar{i,j}*Kbar{i})+(pwasys.Abar{i,j}+pwasys.Bbar{i,j}*Kbar{i})'*pwactrl.Pbar{i};  % Don't use alpha for the center region(s) if you already have the controller
            else
                DV{i,j} = pwactrl.Pbar{i}*(pwasys.Abar{i,j}+pwasys.Bbar{i,j}*Kbar{i})+(pwasys.Abar{i,j}+pwasys.Bbar{i,j}*Kbar{i})'*pwactrl.Pbar{i}+Param.alpha*pwactrl.Pbar{i};
            end
            pwactrl.constraints=pwactrl.constraints+set( DV{i,j}(1:pwasys.n,1:pwasys.n) < 0,['DV' num2str(i) ',' num2str(j) '<0']);
        end
    end
end

for i=1:pwasys.NR,
    %%%%%%%%%%%%%%%%%%%%%%%% Inequality Constraints %%%%%%%%%%%%%%%%%%%

    if ~isfield(pwactrl,'Kbar'),
        if isfield(Param,'LimitK'),
            if ~ismember(i,pwasys.istar) | ~pwactrl.Lin2PWA,
                pwactrl.constraints=pwactrl.constraints+set(Kbar{i}(:) < Param.LimitK,'Kbar<LimitK');
                pwactrl.constraints=pwactrl.constraints+set(Kbar{i}(:) > -Param.LimitK,'Kbar>-LimitK');
            end
        end %LimitK%
    end

    %%%%%%%%%%%%%%%%%%%%%%%%% Negative definite %%%%%%%%%%%%%%%%%%%%%%%%%
    if i~=pwasys.istar,
        if exist('L2Gain'),
            L{i} = sdpvar(pwasys.p+1,pwasys.p+1);
            pwactrl.constraints=pwactrl.constraints+set(0 < L{i}(:) < 1e5,['L{' num2str(i) '}(:)>0']);
            L2Mat{i} = [(pwasys.Abar{i}+pwasys.Bbar{i}*Kbar{i})'*pwactrl.Pbar{i}+pwactrl.Pbar{i}*(pwasys.Abar{i}+pwasys.Bbar{i}*Kbar{i})+(pwasys.Czbar{i}+pwasys.Dzu{i}*Kbar{i})'*(pwasys.Czbar{i}+pwasys.Dzu{i}*Kbar{i})+pwasys.Ebar{i}'*L{i}*pwasys.Ebar{i} pwactrl.Pbar{i}*pwasys.Bwbar{i}+(pwasys.Czbar{i}+pwasys.Dzu{i}*Kbar{i})'*pwasys.Dzw{i}; (pwactrl.Pbar{i}*pwasys.Bwbar{i}+(pwasys.Czbar{i}+pwasys.Dzu{i}*Kbar{i})'*pwasys.Dzw{i})' pwasys.Dzw{i}'*pwasys.Dzw{i}-gammaSQ*eye(q)];
            pwactrl.constraints=pwactrl.constraints+set( L2Mat{i} < 0,['L2Mat' num2str(i) '<0']);
        else

            for j=1:pwasys.NS,
                L{i,j} = sdpvar(pwasys.p+1,pwasys.p+1);
                % S procedure
                pwactrl.constraints=pwactrl.constraints+set(0 < L{i,j}(:) < 1e5,['L{' num2str(i) '}(:)>0']);
                DV{i,j} = pwactrl.Pbar{i}*(pwasys.Abar{i,j}+pwasys.Bbar{i,j}*Kbar{i})+(pwasys.Abar{i,j}+pwasys.Bbar{i,j}*Kbar{i})'*pwactrl.Pbar{i}+Param.alpha*pwactrl.Pbar{i}+pwasys.Ebar{i}'*L{i,j}*pwasys.Ebar{i};
                pwactrl.constraints=pwactrl.constraints+set( DV{i,j} < 0,['DV' num2str(i) '-' num2str(j) '<0']);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%% Continuity constraints %%%%%%%%%%%%%%%%%%%%%%%%
    for j = i:pwasys.NR,
        if ~isempty(pwasys.Fbar{i,j}),
            % Continuity of the control input
            try
                pwactrl.constraints=pwactrl.constraints+set((Kbar{i}-Kbar{j})*pwasys.Fbar{i,j}== 0,['Ki*F' num2str(i) '_' num2str(j) 's=Kis*F' num2str(i) '_' num2str(j) 's']);
            end
        end
    end
end

if ~isfield(pwactrl,'Kbar'),
    if pwactrl.Lin2PWA,
        for i = 1:pwasys.NR,
            if i~=pwasys.istar,
                assign(Kbar{i}(:,1:end-1),Kbar{pwasys.istar(1)}(:,1:end-1));
            end
        end
    end
end

pwactrl.Kbar = Kbar;

Obj = [];
if isfield(Param,'L2Gain'),
    Obj = gammaSQ;
end
solution=solvesdp(pwactrl.constraints,Obj,sdpsettings('usex0',1));
checkset(pwactrl.constraints);

for i=1:pwasys.NR,
    pwactrl.Kbar{i} = double(pwactrl.Kbar{i});
    pwactrl.Pbar{i} = double(pwactrl.Pbar{i});
end

pwactrl.Index = pwasys.Index;
pwactrl.X = pwasys.X;
pwactrl.T = pwasys.T;
pwactrl.xcl = pwasys.xcl;
pwactrl.Ebar = pwasys.Ebar;
pwactrl.gamma = sqrt(double(gammaSQ));