function model=DI_identifier(model)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

L=model.NR;

SA_1=[0 0];
Sa_1=[0 0];
SB_1=[0 0];

for i=1:L
    try
        SA_1=SA_1+size(model.A{i});
    catch
        SA_1=zeros(1,2);
    end

    
    %%
    try
        Sa_1=Sa_1+size(model.a{i});
    catch
        Sa_1=zeros(1,2);
    end
 
    %%
    try
        SB_1=SB_1+size(model.B{i});
    catch
        SB_1=zeros(1,2);
    end
    
end

size_1=SA_1+Sa_1+SB_1;

if  size_1~=0
    pwatype='lower-envelope';
else
    pwatype='null';
end

model.pwatype=pwatype;

