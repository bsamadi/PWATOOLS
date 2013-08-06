function pwacreate(n,m, myname, NR)
%
% pwacreate writes an m-file for approximating a nonlinear funcion with a
% PWA model if the number of the inputs to the function is 3. If pwacreate
% is called with 4 inputs, then an m-file is created for inserting a PWA
% model.
%
% INPUTS
% n      : order of the system
% m      : number of the inputs to the system
% myname : file name which should be created. Example: 'myfilename.m'
% NR     : Number of regions in a PWA model (only if a PWA model is directly used)
%
% OUTPUT
% An m-file called myname which conatins the format of the function that
% should be filled and executed by the user.
%
% run pwacrtnon to see the details of approximating a nonlinear model
%
% run pwacrtpwa to see the details of importing a PWA model.
%
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%

if nargin==3
    ImportNonlinear(n,m, myname);
elseif nargin==4
    ImportPWA(n,m,myname, NR);
else
    fprintf('Too few parameters used for calling pwacreate function\n\n');
end

open(myname)
end
%% local function: ImportNonlinear

function ImportNonlinear(n,m, myname);
fid = fopen(myname,'wt');
fprintf(fid, '%%        ------------------------------------------------------- \n');
fprintf(fid, '%%       |                                                       |\n');
fprintf(fid, '%%       |       Section 1: THIS PART IS FILLED BY THE USER      |\n');
fprintf(fid, '%%       |                                                       |\n');
fprintf(fid, '%%        ------------------------------------------------------- \n\n\n');

fprintf(fid, '%% Nonlinear model dynamics ***  dx/dt=A x + a + f(x) + B(x) u  *** \n');
fprintf(fid, '%% The model is saved in a cell called ''model''.\n\n');
fprintf(fid, '%%                        START OF MODEL DEFINITION\n\n\n');

fprintf(fid,'clear model; \n\n');

fprintf(fid, '%%This line states that the model is nonlinear. You should not change it. \n');
fprintf(fid, 'model.type=''nonlinear'';\n\n');
fprintf(fid, 'model.dim=[%i, %i];  %% Order of the system and number of the inputs \n\n', n, m);

fprintf(fid, '%%---------- Enter matrice A (%i x %i) here\n',n,n);
fprintf(fid, 'model.A=[];\n\n');

fprintf(fid, '%%---------- Enter the affine term a (%i x 1) here\n',n);
fprintf(fid, 'model.aff=[];\n\n');

fprintf(fid, '%%---------- Enter the matrix B(x) (%i x %i) as a function of x here\n',n,m);
fprintf(fid, 'model.Bx=[];\n\n');

fprintf(fid, '%%---------- Enter the domains of %i variables here\n', n);
fprintf(fid, 'model.Domain={');
for i=1:n-1
    fprintf(fid, '[],');
end
fprintf(fid, '[]};\n\n');

fprintf(fid, '%%---------- Enter the function handler here to call the function that generates f(x) (%i x 1)\n', n);
fprintf(fid, 'model.fx=@\n\n');

fprintf(fid, '%%---------- Change the elements whose corresponding variable appear in f(x) to ''1''. For example for\n');
fprintf(fid, '%%           a second order system if only x(1) appears in f(x), then model NonlinearDomain=[1;0]\n');
fprintf(fid, 'model.NonlinearDomain=[');
for i=1:n-1
    fprintf(fid, '0;');
end
fprintf(fid, '0];\n\n');

fprintf(fid, '%%---------- Enter the desired equilibrium point here\n');
fprintf(fid, 'model.xcl=[];\n\n');

fprintf(fid, '%%---------- Enter the griddig type here. Choose from {''Uniform'',''Optimaluniform'',''Multiresolution''}\n');
fprintf(fid, 'model.mtd=[];\n\n');

fprintf(fid, '%%---------- Enter the desired number of partitions\n');
fprintf(fid, 'model.NR=[];\n\n');

fprintf(fid, '%%---------- Enter Rstar region here. Rstar region is OPTIONAL and will be in effect\n');
fprintf(fid, '%%           only if you choose ''Multiresolution'' mode. Example: Rstar=[-.1 .1]\n');
fprintf(fid, 'model.Rstar=[];\n\n');

fprintf(fid, '%% Save the model and run it by pressing (F5) to obtian the PWA pproximation\n\n\n');

fprintf(fid, '%%                        END OF MODEL DEFINITION\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');

fprintf(fid, '%%          ------------------------------------------------------- \n');
fprintf(fid, '%%         |                                                       |\n');
fprintf(fid, '%%         |Section 2: THIS PART SHOULD NOT BE MODIFIED BY THE USER|\n');
fprintf(fid, '%%         |                                                       |\n');
fprintf(fid, '%%          ------------------------------------------------------- \n\n\n');
fprintf(fid, 'rehash;\n');
fprintf(fid, '[Err, model]=NonCheckModel(model);\n');
fprintf(fid, 'if Err==0\n');
fprintf(fid, '    pwasys=PWAComp(model);\n');
fprintf(fid, '    fprintf(''%sNonlinear system saved in nlsys.%s'');\n', '\n','\n\n');
fprintf(fid, '    fprintf(''PWA Differential Inclusion approximation saved in pwainc.%s'');\n', '\n\n');
fprintf(fid, '    fprintf(''PWA approximation saved in pwasys.%s'');\n', '\n\n');
fprintf(fid, 'elseif Err==1\n');
fprintf(fid, '    fprintf(''%sI quit because model has error.%s'');\n', '\n','\n\n');
fprintf(fid, 'end\n');

fclose(fid);

end

%% local function: ImportPWA

function ImportPWA(n,m, myname, NR)

fid = fopen(myname,'wt');
fprintf(fid, '%%          ------------------------------------------------------- \n');
fprintf(fid, '%%          |                                                       |\n');
fprintf(fid, '%%          |       Section 1: THIS PART IS FILLED BY THE USER      |\n');
fprintf(fid, '%%          |                                                       |\n');
fprintf(fid, '%%           ------------------------------------------------------- \n\n\n');

fprintf(fid, '%% PWA model dynamics ***  dy/dt=A{i} y + a{i} + B{i} u  *** \n');
fprintf(fid, '%% The model is saved in a cell called ''model''.\n\n');

fprintf(fid, '%%                        START OF MODEL DEFINITION\n\n\n');
fprintf(fid,'clear model; \n\n');

fprintf(fid, '%%This line states that the model is PWA. You should not change it. \n');
fprintf(fid, 'model.type=''PWA'';\n');
fprintf(fid, 'model.dim=[%i, %i];  %% Order of the system and number of the inputs \n', n, m);
fprintf(fid, 'model.NR=%i;        %% Number of the regions \n\n', NR);

fprintf(fid, '%%---------- Enter matrices A{i} (%i x %i) here\n',n,n);

for i=1:NR
    fprintf(fid, 'model.A{%i}=[];\n', i);
end

fprintf(fid, '\n');
fprintf(fid, '%%---------- Enter the affine term a{i} (%i x 1) here\n',n);

for i=1:NR
    fprintf(fid, 'model.a{%i}=[];\n', i);
end

fprintf(fid, '\n');
fprintf(fid, '%%---------- Enter matrices B{i} (%i x %i) here\n',n, m);

for i=1:NR
    fprintf(fid, 'model.B{%i}=[];\n', i);
end

fprintf(fid, '\n');
fprintf(fid, '%%---------- Enter control gains K{i} (%i x %i) here, if applicable\n',m, n) ;
fprintf(fid, '%%---------- Controller output is construted as u=K{i}x+k{i}\n') ;

for i=1:NR
    fprintf(fid, 'model.K{%i}=zeros(%i,%i);\n', i, m,n);
end
fprintf(fid, '\n');
fprintf(fid, '%%---------- Enter control gains k{i} (%i x 1) here, if applicable\n',m) ;
fprintf(fid, '%%---------- Controller output is construted as u=K{i}x+k{i}\n') ;
for i=1:NR
    fprintf(fid, 'model.k{%i}=zeros(%i,1);\n', i, m);
end

fprintf(fid, '\n');
fprintf(fid, '%%---------- Enter region equation matrices E{i} and e{i} here\n');
for i=1:NR
    fprintf(fid, 'model.E{%i}=[];\n', i);
    fprintf(fid, 'model.e{%i}=[];\n\n', i);
end
% fprintf(fid, '\n');
% fprintf(fid, '%%---------- Enter boundary equations matrices F{i,j} (%i x %i) and f{i,j} (%i x %i) here\n',n, n-1, n,1);
% fprintf(fid, '%%           Leave matrices which refer to non-existing boundaries blank. \n');
% for i=1:NR
%     for j=i+1:NR
%         fprintf(fid, 'model.F{%i,%i}=[];\n', i,j);
%         fprintf(fid, 'model.f{%i,%i}=[];\n\n', i,j);
%      end
% end

fprintf(fid, '\n');
fprintf(fid, '%%---------- Enter the domains of %i variables here\n', n);
fprintf(fid, 'model.Domain={');
for i=1:n-1
    fprintf(fid, '[],');
end
fprintf(fid, '[]};\n', n);

fprintf(fid, '\n');
fprintf(fid, '%%---------- Enter the desired equilibrium point here\n');
fprintf(fid, 'model.xcl=[];\n\n');

fprintf(fid, '%% Save the model and run it by pressing (F5) to export the PWA model\n\n\n');

fprintf(fid, '%%                        END OF MODEL DEFINITION\n\n\n\n\n\n\n\n\n\n\n\n\n');

fprintf(fid, '%%        ------------------------------------------------------- \n');
fprintf(fid, '%%       |                                                       |\n');
fprintf(fid, '%%       |Section 2: THIS PART SHOULD NOT BE MODIFIED BY THE USER|\n');
fprintf(fid, '%%       |                                                       |\n');
fprintf(fid, '%%        ------------------------------------------------------- \n\n\n');
fprintf(fid, 'clc\n');
fprintf(fid, 'rehash;\n');
fprintf(fid, '[Err, model]=PWACheckModel(model);\n');
fprintf(fid, 'if Err==0\n');
fprintf(fid, '    pwasys=PWAExport(model);\n');
fprintf(fid, '    fprintf(''%sPWA approximation saved in pwasys.%s'');\n','\n\n', '\n\n');
fprintf(fid, 'elseif Err==1\n');
fprintf(fid, '    fprintf(''%sI quit because model has error.%s'');\n', '\n','\n\n');
fprintf(fid, 'end\n');


fclose(fid);

end


