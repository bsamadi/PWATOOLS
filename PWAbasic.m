clear model pwasys
echo on;
clc
%                          OVERVIEW OF PWATOOL                       
%
% PWATOOL starts with importing a nonlinear or piecewise-affine (PWA) model.
% It also analyzes the PWA model and synthesizes PWA controllers for the PWA
% model. 
%
% If  the model is nonlinear, PWATOOL approximates it with a PWA model and
% uses PWA Differential Inclusions (PWADI) for analyzing the nonlinear
% model and synthesizing PWA controllers for it.
%
%
pause; %strike any key to continue
clc
%                   SECTION A: DEFINING A PWA MODEL
%
% A piecewise-affine model is described by NR set of linear dynamics
% (including an affine term) in NR regions R_i (i=1,2,...,NR) as follows:
%
%       dx/dt = A{i} x+ a{i}+ B{i} u          if x is in R_i
%
% Each triple [A{i}, a{i}, B{i}] indicates the dynamics of the system in
% region R_i (i=1,2,...,NR). 
%
pause; %strike any key to continue
%
% Technically, we store the dynamics data in array cells. For example if
% the model is called pwasys, then we will have three fields pwasys.A,
% pwasys.a and pwasys.B all with sizes of NR x 1 so that at the i-th row,
% pwasys.A{i}, pwasys,a{i} and pwasys.B{i} would denote the dynamics of 
% the system in region R_i. 
%
%
pause; %strike any key to continue
clc
%                       SECTION B: DEFINING PWADI 
%
% Piecewise-affine differential inclusions builds on PWA modeling. Here,
% instead of only one dynamics in each region R_i, we have a pair of
% dynamics which we may refer to them as lower- and upper-envelopes. A
% PWADI is described by envelopes sigma_k for k=1,2 as follows:
%
%      sigma_k: dx/dt=A{i,k} x+ a{i,k}+ B{i,k} u       if x is in R_i
%                                                      (i=1,2,..., NR) 
%
pause; %strike any key to continue
%
% In fact, each sigma_k is itself a PWA model. However, when we refer to
% PWADI we mainly mean that sigma_1 and sigma_2 together show a similar
% behavior such that they can approximate a nonlinear dynamics (this is
% elaborated more in the following sections).
%
pause; %strike any key to continue
%
% Technically and similar to PWA modeling, we store PWADI dynamics matrices
% in cell arrays. However, here for a PWADI model called pwainc, the fields
% pwainc.A, pwainc.a and pwainc.B are of size NR x 2 so that, at the i-th row
% pwasys.A{i,1} denotes the linear dynamics of sigma_1 and pwasys.A{i,2}
% denotes the linear dynamics of sigma_2.
%
%
pause; %strike any key to continue
clc
%     SECTION C: REGION CONSTRUCTION IN PWA AND PWADI MODELING
%
%
% One related and important issue in PWA and PWADI modeling is how the
% domain of definition is split into regions R_i (i=1,2,...,NR). Each
% region R_i is described by a convex "polyhedron" which, if bounded, it
% may also be called a convex "polytope". 
% see http://www.wordiq.com/definition/Polyhedron for more details.
%
pause; %strike any key to continue
%
% A polyhedron is bounded by several "hyperplanes" which are defined over a
% subset of state space. A hyperplane is described by a linear-affine 
% equation
%
%            E_1*x(1)+E_2*x(2)+...+E_n*x(n)+ e = 0,  
%
% where E=[E_1, E_2, . . .,E_n] is a row vector, e is a scalar and one or
% more coefficients E_i are zero in the above definition.
%
pause; %strike any key to continue
%
% A polyhedron is constructed as the intersection of a finite number of 
% half-spaces 
%
%          E_{i,1}*x(1)+E_{i,2}*x(2)+...+E_{i,n}*x(n)+ e{i} >= 0
%
% where E=[E_{i,j}] is a matrix whose element at the i-th row and j-th
% column is E_{i,j} and e=[e{1}; e{2};...;e{n}] is column vector.
%
pause; %strike any key to continue
clc
% A "slab" region is a polyhedron which is described by two hyperplanes,
% parallel to each other. Therefore, a slab region is bounded only in one
% direction.
%
pause; %strike any key to continue
%
% We refer to region synthesis as "gridding" in our notation. Furthermore, 
% we use the term 'Uniform' to indicate that a gridding has been done
% uniformly on the subset of the state space which appear in hyperplanes.  
% 
pause; %strike any key to continue
%
% We sometime need to approximate the polytopic regions with ellipsoids. An
% ellipsoid equation which approximates a polytopic region is described as
%
%        Ellipsoid Equation:    ||A x + b|| <= 1
%
% where A is a matrix and b is a column vector. 
%
pause; %strike any key to continue
%
% If the regions are slab, a degenerate ellipsoid will approximate the
% region:
%
%        Degenerate Ellipsoid Equation:    ||C x + d|| <= 1
%
% where C is a row vector and d is a scalar.
%
%
pause; %strike any key to continue
clc
%     SECTION D:  APPROXIMATING A NONLINEAR SYSTEM WITH A PWA MODEL
% 
% We explain how a nonlinear model should be defined for PWATOOL so that it
% can be later approximated by a PWADI model.
%
pause; %strike any key to continue
%
% A nonlinear model which is expressed in terms of x is of the form below.
%
%                       dx/dt=A x + a + f(x) + B(x) u
%
% PWATOOL gets the parameters A and a , functions f(x) and B(x) and 
% approximates the nonlinear dynamics with a PWADI model. Here, A 
% represents the linear dynamics, a the affine term, f(x) the
% nonlinearity and B(x) the input gain of the system.
%
% Let us explain the structure of the parameters above and show how PWATOOL
% gets them from the user.
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% | Defining 'A' and 'a'                                                  |
% |-----------------------------------------------------------------------|
% 
% We assume the order of the system is n and there are m inputs to the
% system.
%
% Linear dynamics 'A' (n x n) and affine term 'a' (n x 1) can be described
% by any MATLAB variable of class "double" or any MATLAB expression which
% generates a variable of class "double". The following formats are valid.
% 
%
% A=[0 1; 2 4];           A=rand(3);                A=[0    1
%                                                      2    4];
%
% a=[2; 1];                                          a=[2
%                                                      1];
%                                                   
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% | Defining 'f(x)'                                                       |
% |-----------------------------------------------------------------------|
%                                                  
% f(x), when evaluated, should generate a column vector of size (n x 1). 
% f(x) should be coded separately as a MATLAB function, say 'myfunc_fx',  
% and be called by a function handle in the main code. In the following
% "NonLinearity" would be the variable which contains the function handle
% for generating f(x):
% 
% NonLinearity=@myfunc_fx;
% 
% where 'myfunc_fx.m' is defined as a function script:
% 
% function F=myfunc_fx(x)
%
% F=[f1(x); f2(x); ...; fn(x)];   % fn(.) are scalar functions
% 
pause; %strike any key to continue
clc
%                                          
% |-----------------------------------------------------------------------|
% | Defining 'B(x)'                                                       |
% |-----------------------------------------------------------------------|
% 
% If input gain B(x) is a constant matrix of size (n x m) it can easily be
% defined in a similar way to A and a. The following formats are valid:
%
% B=[2 0; 1 3];          B=[20 ; 1];
% 
% However, in general, B(x) is a function of x which should be defined
% in a similar way to f(x) by calling a function handle. It is important
% that the function which is called generates a matrix of size (n x m). For
% example "InputGain", in the following example, would be the variable which
% contains the function handle for generating B(x):
%
% InputGain=@myfunc_Bx(x)
% 
% where 'myfunc_Bx.m' is defined as follows.
%
% function Bx=myfunc_Bx(x)
%
% Bx=[B11(x) B12(x). . .B1m(x)
%     B21(x) B22(x). . .B2m(x)
%     .      .          .
%     .      .          . 
%     .      .          .
%     Bn1(x) Bn2(x). . .Bnm(x)];   %Bij(x) are scalar functions
%      
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% | Defining The Domain 'R'                                               |
% |-----------------------------------------------------------------------|
% 
% In addition to the nonlinear dynamics parameters, we should also specify
% the domain of the variables and the type of the gridding which we like to
% apply to the domain to obtain the regions in PWA modeling.
%
% A domain is specified by a cell of 'n' elements.
% 
% For example a domian
 
% R={[x1_min x1_max], [x2_min  x2_max], . . . , [xn_min  xn_max]}
% 
% means     x1_min <=  x(1)<=   x1_max, 
%           x2_min <=  x(2)<=   x2_max,
%           
% and so on. 
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% | Domain of the nonlinear dynamics 'NonlinearDomain'                    |
% |-----------------------------------------------------------------------|
%
% In approximating the nonlinear model with a PWA model, PWATOOL needs to
% know which variables appear in the nonlinear function f(x). We store the
% indices of these variables in a variable.
% 
% We place a '1' at the j-th position of a vector of length 'n' for every
% variable x(j) that appear in f(x). For example, we write
%
% NonlinearDomain=[0 1 0 0 1]
%
% to mean that only x(2) and x(5) appear in f(x). In other words, f(x) can
% be written as f(x)=g(x(2), x(5)) where g(.) is a function.
% 
% 
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% | Gridding Type and The Number of Regions                               |
% |-----------------------------------------------------------------------|
% 
% We, next, need to tell PWATOOL how it should grid the domain. The current
% options are
%  
% 1) 'Uniform',
% 2) 'OptimalUniform',
% 3) 'Multiresolution'.
% 
% In 'Uniform' and 'OptimalUniform' gridding, the domain is gridded
% uniformly based on the number of regions (NR) which user chooses. 
% 
% In 'Multiresolution' gridding, the gridding is more flexible in the sense
% that user can set certain points in the regions a priori.
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% | Defining 'Rstar'  ('Multiresolution' mode)                            |
% |-----------------------------------------------------------------------|
% 
% In 'Multiresolution' mode, user can specify a region 'Rstar' which
% encompasses the equilibrium point 'xcl' and will not be split by PWATOOL.  
% We write, for example, 
% 
% Rstar=[-0.1 ; 0.5]
% 
% to state that  -.1 <= xcl <=.5  and that to make sure 'xcl' is within the 
% region 'Rstar' above.
%
% In general 'Rstar' is a matrix of size (* x ND) where ND is the number of
% variables which appear in f(x) and * denotes any number greater than or
% equal to 1. 
%
% Each row of Rstar indicates a point in the domain of f(x) which will 
% later be an edge in the final gridding. Therefore, 
%
%  model.Rstar=[x_1   y_1
%               x_2   y_2
%               x_3   y_3
%               x_4   y_4];
%                 
% shows the points P1=[x_1 y_1], P2=[x_2 y_2], P3=[x_3 y_3] and
% P4=[x_4 y_4] that define a region in the final gridding. This ensures us
% that PWATOOL doesn't split the domain that covers P1 to P4. 
%  
pause; %strike any key to continue
clc
% In brief, if x is defined over the vector space R^n and there are m
% inputs to the system, then the table below gives a summary of the
% structure of the system parameters:
% 
% |------------------------------------------------------------------------
% |PARAMETER   REPRESENTRS          STRUCTURE          Type of data    
% |------------------------------------------------------------------------
% |A           Linear dynamics      n x n  matrix      double    
% |a           Affine term          n x 1  vector      double    
% |f(x)        Nonlinear dynamics   n x 1  vector      function handle    
% |B(x)        Input gain           n x m  matrix      double/function handle    
% |------------------------------------------------------------------------
%
% |------------------------------------------------------------------------
% |PARAMETER                            STRUCTURE
% |------------------------------------------------------------------------
% |Domain                               cell  (1 x n)
% |Domain of the nonlinear dynamics     binary array (1 x n) 
% |Gridding                             {'Uniform', 'OptimalUniform', 'Multiresolution'} 
% |Rstar (OPTIOANL)                     array (1 x 2)   
% |------------------------------------------------------------------------
%
pause % strike any key to continue
clc
%                  SECTION E: IMPORTING A PWA MODEL
% 
% 
% In addition to approximating a nonlinear system with a PWADI model, 
% PWATOOL can directly import a PWA model. 
%
pause % strike any key to continue
%
% |-----------------------------------------------------------------------|
% | Defining a PWA model
% |-----------------------------------------------------------------------|
%
% Consider a PWA model
%
%              dx/dt = A{i} x+ a{i}+ B{i} u  
%
%                  u = K{i} x+ k{i}            (OPTIONAL)
%
% If the controllers K{i} and k{i} are zero or are not defined, the system
% is open-loop.
%
pause; %strike any key to continue
%
% All variables should be of class double. For example for a PWA model of
% order 3 and with 1 input which is defined over a domain with 2 regions
% we will have:
%
% A{1}=[2 1 0; 0 -1 0; -1 1 2];   A{2}=[-1 0 0; 0 0 1; 3 5 0];
% a{1}=[1;1;2];                   a{2}=[0;1;0];
% B{1}=[0;0;1];                   B{2}=[1;0;1];
%
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% | Defining the domain
% |-----------------------------------------------------------------------|
%
% Domains of the variables are defined in a similar way to that of
% nonlinear systems; we use a cell of size (1 x n) to store the lower and
% upper limits of each variable:
%
% R={[x1_min x1_max], [x2_min  x2_max], . . . , [xn_min  xn_max]}
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% | Defining the regions R_i
% |-----------------------------------------------------------------------|
%
% Other than the dynamics of the system, a PWA model definition requires
% the equations of the regions. 
%
% Each region R_i should be defined in terms of linear inequalities as
% follows
%
%                R_i={x | E{i} x + e{i} >= 0 }.    i=1,2,...,NR.
%
% where E{i} is a matrix and e{i} is a column vector. We store matrices
% E{i} and e{i}, similar to the dynamics parameters, in cell arrays.
%
pause; %strike any key to continue
%
% Example: We write 
% 
%                       E{1}=[1 0; -2 0]; and e{1}=[0; 3]
% 
% to mean that          x(1)>=0   and   -2*x(1)+3>=0
% or
% equivalently          0 <=  x(1)  <= 1.5
%
pause; %strike any key to continue
clc
% Table below shows a summary of the system and the regions parameters for
% a PWA model of order n and with m inputs for a gridding with NR regions:
%
% |------------------------------------------------------------------------
% |PARAMETER     REPRESENTRS             STRUCTURE         Size of Each 
% |                                                        Matrix / Vector
% |------------------------------------------------------------------------
% |A             Linear dynamics         NR x 1  cell      n x n
% |a             affine term             NR x 1  cell      n x 1
% |B             Input gain              NR x 1  cell      n x m
% |K             Controller (OPTIONAL)   NR x 1  cell      m x n
% |k             Controller (OPTIONAL)   NR x 1  cell      m x 1 
% |E             Regions equation        NR x 1  cell      * x n
% |e             Regions equations       NR x 1  cell      * x 1 
% |------------------------------------------------------------------------
%
%
pause; %strike any key to continue
clc
%       SECTION F: AGGREGATED VECTOR SPACE Y=[X;1] 
%
%
% For consistency purposes, PWATOOL, also saves the system data in an
% aggregated vector space y=[x;1] in R^(n+1). In producing the results the
% fields 'Abar', 'Bbar' and 'Ebar' denote this notation.
%
% We, briefly, explain the details of this notation.
%
pause; %strike any key to continue
%
% |-----------------------------------------------------------------------|
% |  Abar: dynamics of the system
% |-----------------------------------------------------------------------|
%
% Abar contains the linear dynamics of the system. Let A{i} be the linear
% dynamics of the PWA (or PWADI) system, a{i} the affine term and 'n' the
% order of the system. Then, Abar is constructed as follows:
%
%             Abar{i}=[A{i}       a{i}
%                     zeros(1,n)       1]; 
%
% 
pause; %strike any key to continue
clc
%
% |-----------------------------------------------------------------------|
% |  Bbar: Input gain
% |-----------------------------------------------------------------------|
%
% Bbar contains the input gain of the system. Let B{i} be the input gain
% and 'm' the number of the inputs to the system. Then, Bbar is constructed
% as follows:
%
%          Bbar{i}=[B{i}
%                   zeros(1,m)]; 
%
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  Ebar: Region equations
% |-----------------------------------------------------------------------|
%
% Ebar contains the region equations matrices. Let E{i} and e{i} describe
% the region equation matrices. Then, Ebar is constructed as follows:
%
%           Ebar{i}=[E{i}        e{i}
%                    zeros(1,n)     1]; 
%
%
pause; %strike any key to continue
clc
%
% |-----------------------------------------------------------------------|
% |  Fbar: Region boundaries
% |-----------------------------------------------------------------------|
%
% Fbar contains the boundary equations matrices. Let F{i,j} and f{i,j}
% describe the boundary equation matrices. Then, Fbar is constructed as
% follows:
%    _
%   |
%   |      Fbar{i,j}=[];                      if F{i,j}=[];
%  /
%  \ 
%   |
%   |      Fbar{i,j}=[F{i,j}    f{i,j}        if F{i,j} is not empty
%   |                 zeros(1,n)    1]; 
%   |_
%
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  Kbar: Control gains
% |-----------------------------------------------------------------------|
%
% Kbar contains the control gains. Let K{i} and k{i} be the linear and
% affine control gain in each region. Then Kbar is constructed as follows:
%
%          Kbar{i}=[K{i}  k{i}]   i=1,2,...,NR
% 
% 
% 
pause; %strike any key to continue
clc
% Table below shows a summary of the system parameters in the aggregated
% vector space y=[x;1].
%
% |------------------------------------------------------------------------
% |PARAMETER    STRUCTURE          REPRESENTRS             Size of Each 
% |                                                        Matrix / Vector
% |------------------------------------------------------------------------
% |Abar         NR x 1  cell       Linear dynamics         (n+1) x (n+1)
% |Bbar         NR x 1  cell       Input gain              (n+1) x m
% |Ebar         NR x 1  cell       Regions equation          *   x (n+1)
% |Fbar         NR x NR cell       Boundary equations      (n+1) x (n+1) 
% |------------------------------------------------------------------------
%
pause; %strike any key to return to the main menu
echo off
 


