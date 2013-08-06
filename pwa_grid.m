function W = pwa_grid(arg1,arg2,arg3)

% Make Grid
% W = makegrid(Domain,N) or W = makegrid(a,b,N)

if nargin==3,
    a = arg1;
    b = arg2;
    N = arg3;
    for i=1:length(a),
        Domain{i} = [a(i) b(i)];
    end
else
    Domain = arg1;
    N = arg2;
end

n = length(Domain);
if length(N)==1,
    N = N*ones(n,1);
end

arg = [];
for i = 1:n,
    if length(Domain{i})==2,
        a = Domain{i}(1);
        b = Domain{i}(2);
        x = a:((b-a)/N(i)):b;
    else
        x = Domain{i};
    end
        eval(['x' num2str(i) '=x;']);
        arg = [arg ',x' num2str(i)];
end
arg(1)=[];

if n >1,
    eval(['[' arg ']=ndgrid(' arg ');']);
end

W = [];
for i = 1:n,
    eval(['W=[W reshape(x' num2str(i) ',prod(size(x' num2str(i) ')),1)];']);
end
