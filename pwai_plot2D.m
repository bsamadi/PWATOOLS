function pwai_plot2D(pwainc)

if size(pwainc.W,2)==1,
    figure(20);
    %    plot(pwainc.X,pwainc.Y);
    title('PWA Approximation');
    xlabel('xNL');
    figure(30);
    plot(pwainc.W,pwainc.Z);
    title('Nonlinear Function');
    xlabel('xNL');
    figure(40);
    plot(pwainc.W,pwainc.Err);
    title('Approximation Error');
    xlabel('xNL');
elseif size(pwainc.W,2)==2 && size(pwainc.Z,2)==1,
    figure(10);
    %     triplot(pwainc.T,pwainc.X(:,1),pwainc.X(:,2));
    title('Polytopic Regions');
    xlabel('$x_1$');
    ylabel('$x_2$');
    figure(20);
    %     trimesh(pwainc.T,pwainc.X(:,1),pwainc.X(:,2),pwainc.Y);
    title('PWA Approximation');
    xlabel('$x_1$');
    ylabel('$x_2$');
    
    figure(30);
    nm = sqrt(length(pwainc.W));
    mesh(reshape(pwainc.W(:,1),nm,nm),reshape(pwainc.W(:,2),nm,nm),reshape(pwainc.Z,nm,nm));
    title('Nonlinear Function');
    xlabel('$x_1$');
    ylabel('$x_2$');
    
    figure(40);
    mesh(reshape(pwainc.W(:,1),nm,nm),reshape(pwainc.W(:,2),nm,nm),reshape(pwainc.Err,nm,nm));
    title('Approximation Error');
    xlabel('$x_1$');
    ylabel('$x_2$');
    
elseif size(pwainc.W,2)==2 && size(pwainc.Z,2)==2,

    % Upper bound
    pwainc.vertices = cell(pwainc.NR,1);
    pwainc.patch = cell(pwainc.NR,1);
    for i=1:pwainc.NR,
        GP = pwainc.Ind(:,i);
        pwainc.vertices{i} = [pwainc.GridPoints{1}(GP(1)) pwainc.GridPoints{2}(GP(2));
            pwainc.GridPoints{1}(GP(1)+1) pwainc.GridPoints{2}(GP(2));
            pwainc.GridPoints{1}(GP(1)) pwainc.GridPoints{2}(GP(2)+1);
            pwainc.GridPoints{1}(GP(1)+1) pwainc.GridPoints{2}(GP(2)+1)];
        pwainc.patch{i} = (pwainc.Abar{i,1}*[pwainc.vertices{i} ones(4,1)]')';
    end
    figure(20);
    %     trimesh(pwainc.T,pwainc.X(:,1),pwainc.X(:,2),pwainc.Y(:,1));
    for i=1:pwainc.NR,
        h=patch(pwainc.vertices{i}(:,1),pwainc.vertices{i}(:,2),pwainc.patch{i}(:,1),'c');
        set(h,'Faces',[1 2 4 3]);
    end
    title('PWA Approximation');
    xlabel('$x_1$');
    ylabel('$x_2$');
    hold on;
    nm = sqrt(length(pwainc.W));
    mesh(reshape(pwainc.W(:,1),nm,nm),reshape(pwainc.W(:,2),nm,nm),reshape(pwainc.Z(:,1),nm,nm));
    figure(25);
    %     trimesh(pwainc.T,pwainc.X(:,1),pwainc.X(:,2),pwainc.Y(:,2));
    for i=1:pwainc.NR,
        h=patch(pwainc.vertices{i}(:,1),pwainc.vertices{i}(:,2),pwainc.patch{i}(:,2),'c');
        set(h,'Faces',[1 2 4 3]);
    end
    title('PWA Approximation');
    xlabel('$x_1$');
    ylabel('$x_2$');
    figure(30);
    nm = sqrt(length(pwainc.W));
    mesh(reshape(pwainc.W(:,1),nm,nm),reshape(pwainc.W(:,2),nm,nm),reshape(pwainc.Z(:,1),nm,nm));
    title('Nonlinear Function');
    xlabel('$x_1$');
    ylabel('$x_2$');
    figure(35);
    nm = sqrt(length(pwainc.W));
    mesh(reshape(pwainc.W(:,1),nm,nm),reshape(pwainc.W(:,2),nm,nm),reshape(pwainc.Z(:,2),nm,nm));
    title('Nonlinear Function');
    xlabel('$x_1$');
    ylabel('$x_2$');
    figure(40);
    mesh(reshape(pwainc.W(:,1),nm,nm),reshape(pwainc.W(:,2),nm,nm),reshape(pwainc.Err{1}(:,1),nm,nm));
    title('Approximation Error');
    xlabel('$x_1$');
    ylabel('$x_2$');
    figure(45);
    mesh(reshape(pwainc.W(:,1),nm,nm),reshape(pwainc.W(:,2),nm,nm),reshape(pwainc.Err{1}(:,2),nm,nm));
    title('Approximation Error');
    xlabel('$x_1$');
    ylabel('$x_2$');
    
    % Lower bound
    pwainc.vertices = cell(pwainc.NR,1);
    pwainc.patch = cell(pwainc.NR,1);
    for i=1:pwainc.NR,
        GP = pwainc.Ind(:,i);
        pwainc.vertices{i} = [pwainc.GridPoints{1}(GP(1)) pwainc.GridPoints{2}(GP(2));
            pwainc.GridPoints{1}(GP(1)+1) pwainc.GridPoints{2}(GP(2));
            pwainc.GridPoints{1}(GP(1)) pwainc.GridPoints{2}(GP(2)+1);
            pwainc.GridPoints{1}(GP(1)+1) pwainc.GridPoints{2}(GP(2)+1)];
        pwainc.patch{i} = (pwainc.Abar{i,2}*[pwainc.vertices{i} ones(4,1)]')';
    end
    figure(50);
    %     trimesh(pwainc.T,pwainc.X(:,1),pwainc.X(:,2),pwainc.Y(:,1));
    for i=1:pwainc.NR,
        h=patch(pwainc.vertices{i}(:,1),pwainc.vertices{i}(:,2),pwainc.patch{i}(:,1),'c');
        set(h,'Faces',[1 2 4 3]);
    end
    title('PWA Approximation');
    xlabel('$x_1$');
    ylabel('$x_2$');
    hold on;
    nm = sqrt(length(pwainc.W));
    mesh(reshape(pwainc.W(:,1),nm,nm),reshape(pwainc.W(:,2),nm,nm),reshape(pwainc.Z(:,1),nm,nm));
    figure(55);
    %     trimesh(pwainc.T,pwainc.X(:,1),pwainc.X(:,2),pwainc.Y(:,2));
    for i=1:pwainc.NR,
        h=patch(pwainc.vertices{i}(:,1),pwainc.vertices{i}(:,2),pwainc.patch{i}(:,2),'c');
        set(h,'Faces',[1 2 4 3]);
    end
    title('PWA Approximation');
    xlabel('$x_1$');
    ylabel('$x_2$');
    figure(60);
    nm = sqrt(length(pwainc.W));
    mesh(reshape(pwainc.W(:,1),nm,nm),reshape(pwainc.W(:,2),nm,nm),reshape(pwainc.Z(:,1),nm,nm));
    title('Nonlinear Function');
    xlabel('$x_1$');
    ylabel('$x_2$');
    figure(65);
    nm = sqrt(length(pwainc.W));
    mesh(reshape(pwainc.W(:,1),nm,nm),reshape(pwainc.W(:,2),nm,nm),reshape(pwainc.Z(:,2),nm,nm));
    title('Nonlinear Function');
    xlabel('$x_1$');
    ylabel('$x_2$');
    figure(70);
    mesh(reshape(pwainc.W(:,1),nm,nm),reshape(pwainc.W(:,2),nm,nm),reshape(pwainc.Err{2}(:,1),nm,nm));
    title('Approximation Error');
    xlabel('$x_1$');
    ylabel('$x_2$');
    figure(75);
    mesh(reshape(pwainc.W(:,1),nm,nm),reshape(pwainc.W(:,2),nm,nm),reshape(pwainc.Err{2}(:,2),nm,nm));
    title('Approximation Error');
    xlabel('$x_1$');
    ylabel('$x_2$');
end
