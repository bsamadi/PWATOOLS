function y=PrintError(u, pwactrl, pwasys)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

n=size(pwasys.A{1}, 1);
% u=[x; error_signal]


x=u(1:n);
Er=u(n+1);

domain=pwasys.Domain;
if Er==1
    clc
    y=1;
    fprintf('\nSIMULATION STOPPED.\n');
    fprintf('\n%s method encountered a problem.\n', pwactrl.AppMethod);
    fprintf('At least one of the states is out of bound in simulation.\n\n')
    fprintf('|---------------------------------------------------------\n');
    fprintf('|  x: states             lower-bound       upper-bound        \n');
    fprintf('|---------------------------------------------------------\n');
    fprintf('|  %9.4f              %9.4f          %9.4f\n', x(1), domain{1}(1), domain{1}(2));
    fprintf('|  %9.4f              %9.4f          %9.4f\n', x(2), domain{2}(1), domain{2}(2));
    fprintf('|  %9.4f              %9.4f          %9.4f\n', x(3), domain{3}(1), domain{3}(2));
    fprintf('|  %9.4f              %9.4f          %9.4f\n', x(4), domain{4}(1), domain{4}(2));
    fprintf('|---------------------------------------------------------\n\n');
elseif Er==0
    y=0;
end



