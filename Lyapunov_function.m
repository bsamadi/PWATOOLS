function pwactrl = Lyapunov_function(pwasys, pwactrl, Param)

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
        pwactrl.constraints = pwactrl.constraints+set(Pg > 0.01*I,['P>0.01']);
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