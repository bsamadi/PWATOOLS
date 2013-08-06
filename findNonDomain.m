%% Local function 6: nonliner domain finder

function model=findNonDomain(model)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

[nA,mA]=size(model.A);

fxStock=[];
for j=1:10
    for i=1:nA
        a=model.Domain{i}(1,1);
        b=model.Domain{i}(1,2);
        x(i)=a + (b-a).*rand;
    end
    fx=eval(model.fx);
    fxStock=[fxStock fx];
end
model.NonlinearStateEquations=find(any(fxStock'));








