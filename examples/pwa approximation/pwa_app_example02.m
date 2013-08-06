%% Function definition

nlfun.Handle = @targetfun02;
nlfun.xstar = [0; 0];
nlfun.Domain = {[-pi pi],[-1 1]};
nlfun.Resolution = 10;
nlfun.ObjFun = 'L2';

%% nlfun.Method = 'Uniform';
nlfun.UGR = 7;
pwa_app01 = pwa_uniform(nlfun);

pwa_plot(pwa_app01);

%% nlfun.Method = 'optimaluniform';
nlfun.UGR = 7;
pwa_app02 = pwa_optimal_uniform(nlfun);

pwa_plot(pwa_app02);

figure(40);
zlim([-0.4,0.4]);

%% nlfun.Method = 'MultiResolution';
nlfun.TNR = 45;
pwa_app03 = pwa_split(nlfun);

pwa_plot(pwa_app03);
figure(40);
zlim([-0.4,0.4]);


%% nlfun.Method = 'upperbound';
Option.Type = 'UB';
nlfun.X = pwa_app03.X;
pwa_app04 = pwa_bounds(nlfun,Option);

pwa_plot(pwa_app04);
    
%% nlfun.Method = 'lowerbound';
Option.Type = 'LB';
nlfun.X = pwa_app03.X;
pwa_app05 = pwa_bounds(nlfun,Option);

pwa_plot(pwa_app05);