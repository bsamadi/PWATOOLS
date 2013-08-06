%          -------------------------------------------------------
%          |                                                       |
%          |       Section 1: THIS PART IS FILLED BY THE USER      |
%          |                                                       |
%           -------------------------------------------------------


% PWA model dynamics ***  dy/dt=A{i} y + a{i} + B{i} u  ***
% The model is saved in a cell called 'model'.

%                        START OF MODEL DEFINITION


clear model;

%This line states that the model is PWA. You should not change it.
model.type='PWA';
model.dim=[3, 1];  % Order of the system and number of the inputs
model.NR=4;        % Number of the regions

%---------- Enter matrices A{i} (3 x 3) here
model.A{1}=[0    1.0000         0
            0   -0.0100         0
       0.1507         0         0];
model.A{2}=[ 0    1.0000         0
             0   -0.0100         0
        0.8584         0         0];
model.A{3}=[ 0    1.0000         0
             0   -0.0100         0
        0.1507         0         0];
model.A{4}=[ 0    1.0000         0
             0   -0.0100         0
        0.8584         0         0];

%---------- Enter the affine term a{i} (3 x 1) here
model.a{1}=[0
            0
           -0.6670];
model.a{2}=[ 0
             0
             0];
model.a{3}=[0
            0
            0.6670];
model.a{4}=[0
            0
            0];

%---------- Enter matrices B{i} (3 x 1) here
model.B{1}=[ 0
             1
             0];
model.B{2}=[ 0
             1
             0];
model.B{3}=[ 0
             1
             0];
model.B{4}=[ 0
             1
             0];

%---------- Enter control gains K{i} (1 x 3) here, if applicable
%---------- Controller output is construted as u=K{i}x+k{i}
model.K{1}=zeros(1,3);
model.K{2}=zeros(1,3);
model.K{3}=zeros(1,3);
model.K{4}=zeros(1,3);

%---------- Enter control gains k{i} (1 x 1) here, if applicable
%---------- Controller output is construted as u=K{i}x+k{i}
model.k{1}=zeros(1,1);
model.k{2}=zeros(1,1);
model.k{3}=zeros(1,1);
model.k{4}=zeros(1,1);

%---------- Enter region equation matrices E{i} and e{i} here
model.E{1}=[0.4686         0         0
           -0.7277         0         0];
model.e{1}=[ 0.8834
            -0.6859];

model.E{2}=[-1.0000         0         0
             0.7277         0         0];
model.e{2}=[ 0
             0.6859];

model.E{3}=[-0.4686         0         0
             0.7277         0         0];
model.e{3}=[ 0.8834
            -0.6859];

model.E{4}=[1.0000         0         0
           -0.7277         0         0];
model.e{4}=[0
           0.6859];


%---------- Enter the domains of 3 variables here
model.Domain={[-1.8850    1.8850],[ -20    20],[ -10    10]};

%---------- Enter the desired equilibrium point here
model.xcl=[0
           0
           0];

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
