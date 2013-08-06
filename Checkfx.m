%% Local function 2: Checkfx

function [er, model]=Checkfx(model)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

er=0;
[nA,mA]=size(model.A);

if strcmp(class(model.fx), 'function_handle')
    for i=1:nA
        a=model.Domain{i}(1,1);
        b=model.Domain{i}(1,2);
        x(i)=a + (b-a).*rand;
    end
    fx=feval(model.fx, x);
    [nf, mf]=size(fx);
    if nf~=nA
        fprintf('Sizes of A and f(x) do not match.\n');
        er=1;
    end
    if mf~=1
        fprintf('f(x) should be a column vector.\n')
        er=1;
    end
else
    fprintf('model.fx should be a function handle.\n')
    er=1;
end

%% finding the states which have nonlinear dynamics
if er==0
    fxStock=[];
    for j=1:10
        for i=1:nA
            a=model.Domain{i}(1,1);
            b=model.Domain{i}(1,2);
            x(i)=a + (b-a).*rand;
        end
        fx=feval(model.fx, x);
        fxStock=[fxStock fx];
    end
    
    index=find(any(fxStock'));
    temp=zeros(1, nA);
    temp(index)=1;
    model.NonlinearStateEquations=temp;
end

