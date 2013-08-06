function Ind = pwa_index(GR)

% Creates all the indices for entries of a matrix of size GR
% Usage: for i=IndexMat, do_some_task(i); end

n = length(GR);
Ind = zeros(length(GR),prod(GR));

for i=1:n,

    if 1<i<n,
        m = prod(GR(i+1:end));
        p = prod(GR(1:i-1));
    end
    if i==1,
        m = prod(GR(i+1:end));
        p = 1;
    end
    if i==n,
        m = 1;
        p = prod(GR(1:i-1));
    end
    
    Vec = reshape(ones(m,1)*(1:GR(i)),1,prod(GR(i:end)));
    Ind(i,:) = reshape((ones(p,1)*Vec)',1,prod(GR));
end

% To make it compatible with MATLAB columnwise indexing
Ind = flipud(Ind);