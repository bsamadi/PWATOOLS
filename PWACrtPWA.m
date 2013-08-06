clear model pwasys
echo on;
clc
% 
% *************************************************************************
%  CREATING THE MODEL (PWA)
% *************************************************************************
% 
% Reference:
%
% (1) L. Rodrigues and J. P. How, "Automated Control Design for a
%     Piecewise-Affine Approximation of a Class of Nonlinear Systems," IEEE
%     American Control Conference, pp.3189-3194, Arlington, Virginia,
%     USA, June 2001 
%
% Consider a PWA model in regions R_i (i=1,2,...,NR):
%
%               dx/dt = A{i} x+ a{i}+ B{i} u    
%
%                   u = K{i} x+ k{i}                  (IF APPLICABLE)
%
pause; %strike any key to continue
%
% In order to show how we store the dynamics data of a PWA model , we use 
% an example from (L. Rodrigues and S. Boyd, 2005); we give a PWA modeling
% of a nonlinear resistor.
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% | Creating a m-file for writing the data of the PWA model               |
% |-----------------------------------------------------------------------|
%
% The nonlinear resistor example from (L. Rodrigues and S. Boyd, 2005) is 
% of the order two (n=2) and has one input (m=1). The system domain is
% split into three regions (NR=3). 
%
% Therefore, to start creating the PWA model we, first, type
%
%            pwacreate(2, 1, 'nonlinear_resistor_pwa.m', 3)
%
pause; %strike any key to continue
%
% or in general 
%
%            pwacreate(n, m, 'myfilename.m', NR)
%
% where myfilename.m is the m-file which is to store the PWA model. By 
% doing so, 'myfilename.m' will be created (if does not exist before)
% or rewritten (if existed before).
%
% Now open 'nonlinear_resistor_pwa.m' which we created previously.
%
pause; %strike any key to continue
clc
            open nonlinear_resistor_pwa.m
% 
% We explain how 'nonlinear_resistor_pwa.m' has been filled.
%
pause; %strike any key to continue
%
% nonlinear_resistor_pwa.m has two sections of which, only, Section 1, 
% should be filled by the user. Section 2 should remain as is to let the
% code be executed after completion.
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model dimension and type                                       |
% |-----------------------------------------------------------------------|
%
% The following parameters are set autmatically when the file is created,
% based on the data that user has passed to 'pwacreate'
%
        model.type='PWA';
        model.NR=3;
        model.dim=[2 1];
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.A : linear dynamics                                            |
% |-----------------------------------------------------------------------|
%
% model.A is a cell in which each matrix model.A{i} (n x n) represents the 
% linear dynamics of the PWA model in region R_i  (i=1,2,...,NR). For the
% nonlinear resistor example we will have
%
%
           model.A{1}=[-30.00      -20.00 
                          .05       -0.25];
                       
                       
           model.A{2}=[-30.00      -20.00 
                          .05       -0.20];
 
                       
           model.A{3}=[-30.00      -20.00 
                          .05        0.10];
 
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.a : affine term                                                |
% |-----------------------------------------------------------------------|
%
% model.a is a cell in which each vector model.a{i} (n x 1) represents
% the affine term of the PWA model in region R_i (i=1,2,...,NR). For the
% nonlinear resistor example we will have
%
           model.a{1}=[24 
                        0];
                       
                       
           model.a{2}=[24
                     0.11];
 
                       
           model.a{3}=[24 
                    -0.07];
%
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.B : Input gain
% |-----------------------------------------------------------------------|
%
% model.B is a cell in which each matrix model.B{i} represents the input
% gain of the PWA model in region R_i (i=1,2,...,NR). For the nonlinear
% resistor example we will have
%
           model.B{1}=[20 
                        0];
                       
                       
           model.B{2}=[20
                     0.11];
 
                       
           model.B{3}=[20 
                    -0.07];
%
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.E & model.e: Regions 
% |-----------------------------------------------------------------------|
%
% model.E and model.e define the region equations. In each region R_i let 
% E{i}=model.E{i} and e{i}=model.e{i}. Then, region R_i is 
%
%                R_i={x | E{i} x + e{i} >= 0 }.    i=1,2,...,NR
%
% For the nonlinear resistor example to mean that
%
%                R_1: x(2)<= 0.2
%                R_2: x(2)>= 0.2 & x(2)<=0.6
%                R_3: x(2)>= 0.6
% we write
%
            model.E{1}=[0    0 
                        0   -1];
            model.e{1}=[ 1
                       0.2];
 
                   
            model.E{2}=[0    0  
                        0    1];
            model.e{2}=[1
                        -0.6];
                   
                    
            model.E{3}=[0    1
                        0   -1];
            model.e{3}=[-.2
                        0.6];
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.F & model.f: Regions boundaries (COMPUTED BY PWATOOL)
% |-----------------------------------------------------------------------|
%
% Region boundaries, if applicable, for two regions R_i and R_j is
% described by an expression  F{i,j} s + f{i,j} where F{i,j} is a matrix of
% size (n x n-1) and f{i,j} is a vector of size (n x 1).
%
%    _
%   |   F{i,j}=[];            If R_i and R_j don't have a common boundary
%   |   f{i,j}=[];            with each other.              
%   |
%  /    otherwise
%  \
%   |   
%   |   F{i,j} s + f{i,j}      For parameter 's' in vector space R^(n-1),
%   |_                         contains the common boundary of R_i and R_j. 
%                               
%
% PWATOOL computes the matrices F and vectors f from the region equations
% matrices E and e. Therefore, user need not enter F and f separately.
%
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.xcl: Equilibrium point
% |-----------------------------------------------------------------------|
%
% Further analysis or synthesis processes for a PWA model requires us to
% define an equilibrium point for the system. For the nonlinear resistor
% example we use the open-loop equilibrium point:
%
            model.xcl=[0.371428571428570
                       0.642857142857146];
%
pause; %strike any key to continue
%
% |-----------------------------------------------------------------------|
% |  model.Domain: Domain of the variables
% |-----------------------------------------------------------------------|
%
% The domain of the variables is cell of size (1 x n). For the nonlinear
% resistor example with -20<= x(1) <=20  and  -20,000<= x(2) <=-20,000  
% we will have:
%
           model.Domain={[-20 20],[-2e4 2e4]};
%
%
% Since at this time the system is open-loop we do not enter any value for
% matrices K{i} and k{i} and we stop at this point. 
%
%
%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                %       MODELING IS COMPLETE        % 
%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  Checking and saving the data
% |-----------------------------------------------------------------------|
%
% Now, save the changes in the file 'nonlinear_resistor_pwa.m' and run it
% by pressing F5 or typing 'run nonlinear_resistor_pwa'.
%
% PWATOOL imports the data that user has entered and checks for dimension
% consistency between the model parameters by calling the function
%
% [Err, model]=PWACheckModel(model);
%
pause; %strike any key to continue
 
              [Err, model]=PWACheckModel(model);
%
%
pause; %strike any key to continue
%
% If there is no error (Err=0), function 'PWAExport' is, then, called to
% build the PWA model.
%
%             pwasys=PWAExport(model);
%
%
pause; %strike any key to continue
 
              pwasys=PWAExport(model);
  
% 
% |-----------------------------------------------------------------------|
% |  Outputs and results
% |-----------------------------------------------------------------------|
%
% PWATOOL stores the variables that user has entered, in a variable called 
% 'pwasys'.
%
pwasys
%
%
pause; %strike any key to continue
clc
%
% The table below shows a summary of the fields of 'pwasys'
%      
% |------------------------------------------------------------------------
% |PARAMETER    STRUCTURE          REPRESENTRS             Size of Each 
% |                                                        Matrix / Vector
% |------------------------------------------------------------------------
% |A            NR x 1  cell       Linear dynamics         n x n
% |a            NR x 1  cell       affine term             n x 1
% |B            NR x 1  cell       Input gain              n x m
% |E            NR x 1  cell       Regions equation        * x n
% |e            NR x 1  cell       Regions equations       * x 1 
% |F            NR x NR cell       Boundary equations      n x (n-1) 
% |f            NR x NR cell       Boundary equations      n x 1 
% |K            NR x 1            (IF APPLICABLE)          m x n
% |k            NR x 1            (IF APPLICABLE)          m x 1 
% |Abar         NR x 1  cell       Linear dynamics         (n+1) x (n+1)
% |Bbar         NR x 1  cell       Input gain              (n+1) x m
% |Ebar         NR x 1  cell       Regions equation          *   x (n+1)
% |Fbar         NR x NR cell       Boundary equations      (n+1) x (n+1) 
% |NR           scalar             Number of regions       1 x 1 
% |type         string             Type of the model       'pwadi' or 'lower-envelope' 
% |------------------------------------------------------------------------
%
% *: The number of rows in matrix E and vector e depend on the gridding 
% type. If the system is slab, then [2, n]=size(E) and [2, 1]=size(e); if 
% the gridding type is uniform and polytopic, then [n+1, n]=size(E) and
% [n+1, 1]=size(e); In general, a  polytopic region with 'h' vertices
% needs h inequality to be described.  
 
pause; %strike any key to return to the main menu
echo off
 


