%          -------------------------------------------------------
%          |                                                       |
%          |       Section 1: THIS PART IS FILLED BY THE USER      |
%          |                                                       |
%           -------------------------------------------------------


% PWA model dynamics ***  dy/dt=A{i} y + a{i} + B{i} u  ***
% The model is saved in a cell called 'model'.

%                        START OF MODEL DEFINITION


clear model pwasys;

%This line states that the model is PWA. You should not change it.

%This line states that the model is PWA. You should not change it. 
model.type='PWA';
model.dim=[2, 1];  % Order of the system and number of the inputs 
model.NR=3;        % Number of the regions 

% -------- LINEAR DYNAMICS (2 x 2)
%
model.A{1}=[-30.00      -20.00
                .05       -0.25];


model.A{2}=[-30.00      -20.00
                .05       -0.20];


model.A{3}=[-30.00      -20.00
                .05        0.10];

%

%--------- AFFINE DYDNAMICS (2 x 1)
%
%
model.a{1}=[24
    0];


model.a{2}=[24
    0.11];


model.a{3}=[24
    -0.07];

%
%--------- INPUT GAIN (2 x 1)
%
model.B{1}=[20
    0];


model.B{2}=[20
    0];


model.B{3}=[20
    0];


%---------- EQUILIBRIUM POINT
%
model.xcl=-model.A{2}\model.a{2};
%
model.Domain={[-20 20],[-2e4 2e4]};
%
%

%
            model.E{1}=[0    .5*1e-4
                        0   -1];
            model.e{1}=[ 1
                       0.2];
 
                   
            model.E{2}=[0    -.5*1e-4  
                        0    1];
            model.e{2}=[1
                        -0.6];
                   
                    
            model.E{3}=[0    1
                        0   -1];
            model.e{3}=[-.2
                        0.6];

model.K{1}=zeros(1,2);
model.K{2}=zeros(1,2);
model.K{3}=zeros(1,2);

%---------- Enter control gains k{i} (1 x 1) here, if applicable
%---------- Controller output is construted as u=K{i}x+k{i}
model.k{1}=zeros(1,1);
model.k{2}=zeros(1,1);
model.k{3}=zeros(1,1);


% Save the model and run it by pressing (F5) to export the PWA model


%                        END OF MODEL DEFINITION








%        -------------------------------------------------------
%       |                                                       |
%       |Section 2: THIS PART SHOULD NOT BE MODIFIED BY THE USER|
%       |                                                       |
%        -------------------------------------------------------
clc
rehash;
[Err, model]=PWACheckModel(model);
if Err==0
    pwasys=PWAExport(model);
    fprintf('\n\nPWA approximation saved in pwasys.\n\n');
elseif Err==1
    fprintf('\nI quit because model has error.\n\n');
end