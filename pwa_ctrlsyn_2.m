 function pwactrl = pwa_ctrlsyn_2(pwasys, pwactrl, Param)

% PWA Controller design
% PWQ stability analysis
%
% pwasys
% xbarDot = Abar{i}*xbar + Bbar{i}*u
%
% pwasys.Abar = {Abar1, ..., Abarn}
% pwasys.Bbar = {Bbar1, ..., Bbarn}
%
% Polytopic region i: {x| Ebar{i}*x > 0 }
%
%
% Boundary information table:
% Boundary between Ri and Rj : {x| xbar=Fbar{i,j}*sbar }
%

xcl = pwasys.xcl;
xclbar = [xcl;1];

[NR NS] = size(pwasys.Abar);    % Number of Systems, Number of Regions
n = size(pwasys.Abar{1},1)-1;   % Number of state variables
m = size(pwasys.Bbar{1},2);     % Number of inputs

if isfield(Param,'LimitK'),     % The upper bound for controller gains
    LimitK = Param.LimitK;
end
if isfield(Param,'LimitP'),     % The upper bound for Lyapunov matrix P entries
    LimitP = Param.LimitP;
end

if isfield(Param,'Lyapunov'),   % Type of the Lyapunov function
    Lyapunov = Param.Lyapunov;
else
    Lyapunov = 'PWQ';           % Piecewise Quadratic Lyapunov functinons are considered by default
end

if isfield(Param,'Poles'),      % Desired poles for the linear controller in the center region
    Poles = Param.Poles;
elseif isfield(Param,'KLin'),   % Gain for the linear controller in the center region
    KLin = Param.KLin;
elseif isfield(Param,'QLin'),   % Q and R for the LQR controller in the center region
    QLin = Param.QLin;
    RLin = Param.RLin;
end

if isfield(Param,'alpha'),      % VDot <-\alpha V
    alpha=Param.alpha;
else
    alpha = 0;
end

istar = [];
Abar = pwasys.Abar;
Bbar = pwasys.Bbar;
Ebar = pwasys.Ebar;
Fbar = pwasys.Fbar;
for i=1:NR,
    for j=1:NS,
        A{i,j} = Abar{i,j}(1:end-1,1:end-1);
        a{i,i} = Abar{i,j}(1:end-1,end);
        B{i,j} = Bbar{i,j}(1:end-1,:);
    end
    Echeck{i} = [Ebar{i};zeros(1,n) 1];
    xcl_is_inside_Ri = all(Ebar{i}*xclbar>=0-1e-7);
    if xcl_is_inside_Ri,
        istar = union(istar,i);                         % Center region(s)
    end
end

Lin2PWA = 1;
for i = istar,
    if exist('Poles'),
        K{i} = -acker(A{i},B{i},Poles);
    elseif e xist('KLin'),
        K{i} = KLin;
    elseif exist('QLin'),
       % K{i} = -lqr(A{i},B{i},QLin,RLin);
        K{i}=robust_LQR(pwasys, pwactrl, Param);
    else
        Lin2PWA = 0;
    end
end

yalmip('clear');
constraints=set([]);
p = size(Ebar{1},1);

if ~isfield(pwactrl,'Pbar'),                        % Is Lyapunov function given?
    % Lyapunov function                             % No
    if strcmp(Lyapunov,'Global'),                   % Global Lyapunov function
        Pg = sdpvar(n,n);
        for i=1:NR,
            P{i} = Pg;
            r{i} = 0;
            % Constrainted Pbar
            Pbar{i} = [P{i} -P{i}*xcl; -xcl'*P{i} xcl'*P{i}*xcl+r{i}];
        end
    else,                                           % Piecewise quadratic Lyapunov function
        for i=1:NR,
            if ismember(i,istar),
                P{i} = sdpvar(n,n);
                Pbar{i} = [P{i} -P{i}*xcl; -xcl'*P{i} xcl'*P{i}*xcl];
            else
                Pbar{i} = sdpvar(n+1,n+1);
            end
        end
    end
else                                                % Lyapunov function is given
    for i=1:NR,
        if ismember(i,istar),
            if Lin2PWA,                             % If there is also a linear controller for the center region, use the controller but not the given Lyapunov function
                P{i} = sdpvar(n,n);
                Pbar{i} = [P{i} -P{i}*xcl; -xcl'*P{i} xcl'*P{i}*xcl];
            else
                if iscell(pwactrl.Pbar),
                    Pbar = pwactrl.Pbar;
                    P{i} = Pbar{i}(n,n);
                else,
                    Pbar{i} = pwactrl.Pbar;
                    P{i} = Pbar{i}(n,n);
                end
            end
        else,
            if iscell(pwactrl.Pbar),
                Pbar{i} = pwactrl.Pbar{i};
            else,
                Pbar{i} = pwactrl.Pbar;
            end
        end
    end
end

if ~isfield(pwactrl,'Kbar'),                        % If a controller is given, use it.
    for i=1:NR,
        Kbar{i} = sdpvar(m,n+1,'full');
    end
else,
    Kbar = pwactrl.Kbar;
end

%%%%%%%%%%%%%%%%%%%%%%%% The center region %%%%%%%%%%%%%%%%%%%
I = eye(n);
One = zeros(n+1);
One(end,end)=1;

K(istar)
if Lin2PWA,
    for i=istar,
        k = sdpvar(m,1);
        %k=[0 0 ]';
        Kbar{i} = [K{i} k];
        
    end
end

for i=istar,
    K{i} = Kbar{i}(:,1:n);
    for j=1:NS,
%         if ~isfield(pwactrl,'Kbar'),
%             constraints=constraints+set( (Abar{i,j}+Bbar{i,j}*Kbar{i})*[xcl;1] == 0,['(Abar+Bbar*Kbar)*xclbar=0']); % The desired equilibrium point
%         end
        if Lin2PWA,
            DV{i,j} = Pbar{i}*(Abar{i,j}+Bbar{i,j}*Kbar{i})+(Abar{i,j}+Bbar{i,j}*Kbar{i})'*Pbar{i};  % Don't use alpha for the center region(s) if you already have the controller
        else
            DV{i,j} = Pbar{i}*(Abar{i,j}+Bbar{i,j}*Kbar{i})+(Abar{i,j}+Bbar{i,j}*Kbar{i})'*Pbar{i}+alpha*Pbar{i};
        end
        constraints=constraints+set( DV{i,j}(1:n,1:n) < 0,['DV' num2str(i) ',' num2str(j) '<0']);    % Vdot(+alpha V)<0
    end
    if ~isfield(pwactrl,'Pbar'),
        constraints=constraints+set(P{i} > 0*I,['P' num2str(i) '>0']);
        if exist('LimitP'),
            constraints=constraints+set(P{i}(:) < LimitP,['P' num2str(i) '<' num2str(LimitP) '*I']);
        end
    end
end

for i=1:NR,
    %%%%%%%%%%%%%%%%%%%%%%%% Inequality Constraints %%%%%%%%%%%%%%%%%%%
    if ~isfield(pwactrl,'Pbar'),
        if strcmp(Lyapunov,'PWQ'),
            if i~=istar,
                Z{i} = sdpvar(p+1,p+1);
                %S procedure
                constraints=constraints+set(1e2 < Z{i}(:),['Z{' num2str(i) '}(:)>0']);
                constraints=constraints+set(Pbar{i} -Echeck{i}'*Z{i}*Echeck{i} > 0,['Pbar' num2str(i) '>0']);
                if exist('LimitP'),
                    constraints=constraints+set(Pbar{i}(:) < LimitP,['Pbar' num2str(i) '<' num2str(LimitP) '*I']);
                end
            end
        end %PWQ%
    end

    if ~isfield(pwactrl,'Kbar'),
        if exist('LimitK'),
            if i~=istar | ~Lin2PWA,
                constraints=constraints+set(Kbar{i}(:) < LimitK,'Kbar<LimitK');
                constraints=constraints+set(Kbar{i}(:) > -LimitK,'Kbar>-LimitK');
            end
        end %LimitK%
    end

    %%%%%%%%%%%%%%%%%%%%%%%%% Negative definite %%%%%%%%%%%%%%%%%%%%%%%%%
    if i~=istar,
        for j=1:NS,
            L{i} = sdpvar(p+1,p+1);
            %S procedure
            constraints=constraints+set(1e2 < L{i}(:),['L{' num2str(i) '}(:)>0']);
            DV{i,j} = Pbar{i}*(Abar{i,j}+Bbar{i,j}*Kbar{i})+(Abar{i,j}+Bbar{i,j}*Kbar{i})'*Pbar{i}+alpha*Pbar{i}+Echeck{i}'*L{i}*Echeck{i};
            constraints=constraints+set( DV{i,j} < 0,['DV' num2str(i) '-' num2str(j) '<0']);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%% Continuity constraints %%%%%%%%%%%%%%%%%%%%%%%%
    for j = i:NR,
        if ~isempty(Fbar{i,j}),
            % Continuity of the control input
            if ~isfield(pwactrl,'Kbar'),
                constraints=constraints+set((Kbar{i}-Kbar{j})*Fbar{i,j}== 0,['Ki*F' num2str(i) '_' num2str(j) 's=Kis*F' num2str(i) '_' num2str(j) 's']);
            end
            if ~isfield(pwactrl,'Pbar'),
                if strcmp(Lyapunov,'PWQ'),
                    % Continuity of the Lyapunov function
                    constraints=constraints+set(Fbar{i,j}'*(Pbar{i}-Pbar{j})*Fbar{i,j} == 0,['F' num2str(i) '_' num2str(j) 's*Pbari*F' num2str(i) '_' num2str(j) 's=F' num2str(i) '_' num2str(j) 's*Pbaris*F' num2str(i) '_' num2str(j) 's']);
                end
            end
        end
    end
end
if ~isfield(pwactrl,'Kbar'),
    if Lin2PWA,
        for i = 1:NR,
            if i~=istar,
                setsdpvar(Kbar{i}(:,1:end-1),Kbar{istar(1)}(:,1:end-1));
            end
        end
    end
end



% Bita=sdpvar(1);
% constraints =constraints + set(Bita > 0);
% for i=1:NR
%     for j=1:NS
%         cl_loop_dif_norm{i,j}= norm(Abar{i,j}+Bbar{i,j}*Kbar{i}-(Abar{istar,j}+Bbar{istar,j}*Kbar{istar}));
%         constraints=constraints+ set(cl_loop_dif_norm{i,j}< Bita);
%     end
% end

obj = [];

% for i=1:NR,
%    obj = obj+norm(Kbar{i});
%  end

solvesdp(constraints,obj,sdpsettings('usex0',1, 'shift', 1e-7))
%solvesdp(constraints)
checkset(constraints)

pwactrl.xcl = xcl;
pwactrl.istar = istar;

for i=1:NR,
    pwactrl.Kbar{i} = double(Kbar{i});
    %pwactrl.L{i}=double(L{i});
    try,
        pwactrl.Pbar{i} = double(Pbar{i});
    end
end

pwactrl.Index = pwasys.Index;
pwactrl.X = pwasys.X;
pwactrl.T = pwasys.T;
pwactrl.xcl = pwasys.xcl;
pwactrl.Ebar = pwasys.Ebar;