%        ------------------------------------------------------- 
%       |                                                       |
%       |       Section 1: THIS PART IS FILLED BY THE USER      |
%       |                                                       |
%        ------------------------------------------------------- 


% % % Nonlinear model dynamics ***  dx/dt=A x + a + f(x) + B(x) u  *** 
% % % The model is saved in a cell called 'model'.
% % 
%                        START OF MODEL DEFINITION

clear
clc

clear model; 

%This line states that the model is nonlinear. You should not change it. 
model.type='nonlinear';

model.dim=[2, 0];  % Order of the system and number of the inputs 
%---------- Enter matrice A (2 x 2) here
model.A=[
                         0                         1 
                         -1                        1 ];

%---------- Enter the affine term a (2 x 1) here
model.aff=[0;0];

%---------- Enter the matrix B(x) (2 x 2) as a function of x here
model.Bx=[
                         0                         0
                         0                         0];

%---------- Enter the domains of 2 variables here
model.Domain={[0 10],[0 10]};

%---------- Enter the function handler here to call the function that generates f(x) (2 x 1)
model.fx=@vdp_fun;

%---------- Change the elements whose corresponding variable appear in f(x) to '1'. For example for
%           a second order system if only x(1) appears in f(x), then model NonlinearDomain=[1;0]

model.NonlinearDomain=[1;1]; 

%---------- Enter the desired equilibrium point here
model.xcl=[0;0];

%---------- Enter the griddig type here. Choose from {'Uniform','Optimaluniform','Multiresolution'}
model.mtd='Uniform';

%---------- Enter the desired number of evaluation points for each variable
model.Resolution = 20;

%---------- Enter the desired number of partitions for each variable
model.NR = 5;

%---------- Enter Rstar region here. Rstar region is OPTIONAL and will be in effect
%           only if you choose 'Multiresolution' mode. Example: Rstar=[-.1 .1]
model.Rstar=[];

% Save the model and run it by pressing (F5) to obtian the PWA pproximation

%                        END OF MODEL DEFINITION







%          ------------------------------------------------------- 
%         |                                                       |
%         |Section 2: THIS PART SHOULD NOT BE MODIFIED BY THE USER|
%         |                                                       |
%          ------------------------------------------------------- 


rehash;
[Err, model]=NonCheckModel(model);
if Err==0
    [pwainc, pwasys, nlsys]=PWAComp(model);
    fprintf('\nNonlinear system saved in nlsys.\n\n');
    fprintf('PWA Differential Inclusion approximation saved in pwainc.\n\n');
    fprintf('PWA approximation saved in pwasys.\n\n');
elseif Err==1
    fprintf('\nI quit because model has error.\n\n');
end

pwa_plot(pwasys);