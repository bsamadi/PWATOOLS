%% Local function 3: CheckBx

function er=CheckBx(model)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

er=0;
x=1; % any number for evaluating B(x)
[nA,mA]=size(model.A);

if strcmp(class(model.Bx), 'function_handle')
    % if model.Bx is a function handle
    for i=1:nA
        a=model.Domain{i}(1,1);
        b=model.Domain{i}(1,2);
        x(i)=a + (b-a).*rand;
    end
    Bx=feval(model.Bx, x);
    [nB, mB]=size(Bx);
else % otherwise 
    [nB, mB]=size(model.Bx);
end

if nB~=nA
    fprintf('Sizes of A and Bx do not match.\n');
    er=1;
end

end

