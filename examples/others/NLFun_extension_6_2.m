function F = NLFun_extension_6_2(X)

x=X(1);

if (-4 <= x & x< -2)
    F=[0;-1.8/2*(x+2)-1];
elseif (-2 <= x & x<-1) 
    F=[0;-1];
elseif (-1 <= x & x < 1)
    F=[0;x];
elseif (1 <= x  & x < 2)
    F=[0;1];
elseif  (2 <= x  & x <= 4)
    F=[0;-1.8/2*(x-2)+1];
end
    
