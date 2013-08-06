function  pwasys=PWAExport(model)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

clear pwasys

for i=1:model.NR
    pwasys.A{i,1}=model.A{i};
    pwasys.a{i,1}=model.a{i};
    pwasys.B{i,1}=model.B{i};
end

pwasys.E=model.E;
pwasys.e=model.e;
pwasys.type=model.pwatype;

try
    pwasys.K=model.K;
    pwasys.k=model.k;
end
pwasys.NR=model.NR;
pwasys.xcl=model.xcl;
pwasys.Domain=model.Domain;

%%
pwatype=model.pwatype;
if strcmp(pwatype, 'lower-envelope')
    col_index=1;
elseif strcmp(pwatype, 'null')
    col_index=[];
end
NR=model.NR;

n=model.dim(1);
m=model.dim(2);

for i=1:NR
    try
        pwasys.Kbar{i}=[model.K{i} model.k{i}];
            
    end
    pwasys.Ebar{i}=[model.E{i}  model.e{i}
        zeros(1,n)           1];
    for j=col_index
        pwasys.Abar{i,j}=[model.A{i} model.a{i}
                         zeros(1,n+1)        ];
        
        pwasys.Bbar{i,j}=[model.B{i}
            zeros(1,m)];
    end
end

pwasys.Fbar=FbarBuilder(pwasys);
pwasys.F{1,1}=[];
pwasys.f{1,1}=[];
for i=1:NR
    for j=i+1:NR
        if ~isempty(pwasys.Fbar{i,j})
            pwasys.F{i,j} = pwasys.Fbar{i,j}(1:end-1,1:end-1);
            pwasys.F{j,i} = pwasys.F{i,j};
            pwasys.f{i,j} = pwasys.Fbar{i,j}(1:end-1,end);
            pwasys.f{j,i} = pwasys.f{i,j};
            
        end
    end
end



