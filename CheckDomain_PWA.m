%% Local function 5: CheckDomain

function er=CheckDomain_PWA(model)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

er=0;
L=model.NR;
n=model.dim(1);

for i=1:n
    try
        [a, b]=size(model.Domain{i});
        if max(a, b) <2
            fprintf('Domain of all variables should be defined.\n');
            er=1;
        else
            A=strcmp(class(model.Domain{i}(1)), 'double');
            B=strcmp(class(model.Domain{i}(2)), 'double');
            if ~(A & B)
                fprintf('Real variables should be used to define Domains.\n');
                er=1;
                return;
            end
        end
    catch
        fprintf('Domain of all variables should be defined.\n');
        er=1;
        return;
    end
    
end

