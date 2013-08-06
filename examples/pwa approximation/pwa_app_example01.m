%% Function definition

nlfun.Handle = @targetfun01;
nlfun.xstar = 0;
nlfun.Domain = {[-pi pi]};
nlfun.Resolution = 100;
nlfun.ObjFun = 'L2';

%% nlfun.Method = 'Uniform';
nlfun.UGR = 7;
pwa_app01 = pwa_uniform(nlfun);

%% nlfun.Method = 'optimaluniform';
nlfun.UGR = 7;
pwa_app02 = pwa_optimal_uniform(nlfun);

%% nlfun.Method = 'MultiResolution';
nlfun.TNR = 7;
pwa_app03 = pwa_split(nlfun);

%% nlfun.Method = 'upperbound';
Option.Type = 'UB';
nlfun.X = pwa_app03.X;
pwa_app04 = pwa_bounds_single(nlfun,Option);
    
%% nlfun.Method = 'lowerbound';
Option.Type = 'LB';
nlfun.X = pwa_app03.X;
pwa_app05 = pwa_bounds_single(nlfun,Option);
    
%% Plot the results

figure(100);
plot(pwa_app01.W, pwa_app01.Z, pwa_app01.X, pwa_app01.Y);
set(gca,'XTick',pwa_app01.X, 'XGrid','on','XLim',[-pi pi]);
xlabel('x');
ylabel('f(x)');
title('Uniform Grid, Obj. Fun. = 0.2611');
make_tex('uniform');

figure(200);
plot(pwa_app02.W, pwa_app02.Z, pwa_app02.X, pwa_app02.Y);
set(gca,'XTick',pwa_app02.X, 'XGrid','on','XLim',[-pi pi]);
xlabel('x');
ylabel('f(x)');
title('Optimal Grid, Obj. Fun. = 0.0519');
make_tex('optimal'); 

figure(300);
plot(pwa_app03.W, pwa_app03.Z, pwa_app03.X, pwa_app03.Y);
set(gca,'XTick',pwa_app03.X, 'XGrid','on','XLim',[-pi pi]);
xlabel('x');
ylabel('f(x)');
title('MultiResolution, Obj. Fun. = 0.0287');
make_tex('multi');

figure(400);
plot(pwa_app04.W, pwa_app04.Z, pwa_app04.X, pwa_app04.Y, pwa_app05.X, pwa_app05.Y);
set(gca,'XTick',pwa_app04.X, 'XGrid','on','XLim',[-pi pi]);
xlabel('x');
ylabel('f(x)');
title('MultiResolution Envelope');
make_tex('env');
