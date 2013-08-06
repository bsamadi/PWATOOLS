function t = pwa_region(x,CellInfo)

% i = pwa_region(x,CellInfo)

% CellInfo.X
% CellInfo.T
% CellInfo.Index Optional
% CellInfo.Ebar  Optional

t=0; % initial output. To prevent a NAN output. 0 would in fact imply an
     % out of bound error

ZeroMinus = -1e-15;

if isfield(CellInfo,'T'),

    n = size(CellInfo.T,2)-1;


    if isfield(CellInfo,'Index')
        N = size(CellInfo.Ebar{1},2);
        if size(x,2)==N-1,
            X = x(:,CellInfo.Index);
        elseif size(x,1)==N-1,
            x = x';
            X = x(:,CellInfo.Index);
        else
            if n>1,
                error(['x should have ' num2str(n) ' columns.']);
            else
                error(['x should have ' num2str(n) ' column.']);
            end
        end
    else
        if size(x,2)==n,
            X = x;
        elseif size(x,1)==n,
            X = x';
        else
            if n>1,
                error(['x should have ' num2str(n) ' columns.']);
            else
                error(['x should have ' num2str(n) ' column.']);
            end
        end
    end

    if size(CellInfo.X,2)~=n,
        if size(CellInfo.X,2)>1,
            error(['x should have ' num2str(size(CellInfo.X,2)) ' columns.']);
        else
            error(['x should have ' num2str(size(CellInfo.X,2)) ' column.']);
        end

    end

%     if n == 2,
%         t = tsearch(CellInfo.X(:,1),CellInfo.X(:,2),CellInfo.T,X(:,1),X(:,2));
%         if isnan(t),
%             t = tsearchn(CellInfo.X,CellInfo.T,X);
%         end
%     else
%         t = tsearchn(CellInfo.X,CellInfo.T,X);
%     end

    t = pointLocation(CellInfo.dt,X);
    
else
    for i=1:length(CellInfo.Ebar),
        if all(CellInfo.Ebar{i}*[x;1]>=ZeroMinus),
            t = i;
            break;
        end
    end
end