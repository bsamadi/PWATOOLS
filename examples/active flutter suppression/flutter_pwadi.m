%        ------------------------------------------------------- 
%       |                                                       |
%       |       Section 1: THIS PART IS FILLED BY THE USER      |
%       |                                                       |
%        ------------------------------------------------------- 


% Nonlinear model dynamics ***  dx/dt=A x + a + f(x) + B(x) u  *** 
% The model is saved in a cell called 'model'.

%                        START OF MODEL DEFINITION


clear model; 

%This line states that the model is nonlinear. You should not change it. 
model.type='nonlinear';

model.dim=[4, 2];  % Order of the system and number of the inputs 
%---------- Enter matrice A (4 x 4) here
model.A=[
                         0                         0                         1                         0
                         0                         0                         0                         1
         -293.269207891109          -100.58975526267         -5.90266024354767        -0.405417880783147
          1885.94659843821          743.792579996859          34.7278766431219          2.46868729488239];

%---------- Enter the affine term a (4 x 1) here
model.aff=[0;0;0;0];

%---------- Enter the matrix B(x) (4 x 2) as a function of x here
model.Bx=[
                         0                         0
                         0                         0
         -7606.77789655331         -7642.55122635784
          14250.4450732197          9021.91847103964];

%---------- Enter the domains of 4 variables here
model.Domain={[-1 1],[-1 1],[-1 1],[-5 5]};

%---------- Enter the function handler here to call the function that generates f(x) (4 x 1)
model.fx=@Flutter_Nonlinearity;

%---------- Change the elements whose corresponding variable appear in f(x) to '1'. For example for
%           a second order system if only x(1) appears in f(x), then model NonlinearDomain=[1;0]
model.NonlinearDomain=[0;1;0;0];

%---------- Enter the desired equilibrium point here
model.xcl=[0;0;0;0];

%---------- Enter the griddig type here. Choose from {'Uniform','Optimaluniform','Multiresolution'}
model.mtd=['Uniform'];

%---------- Enter the desired number of partitions
model.NR=[6];

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
