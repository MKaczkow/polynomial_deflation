clc
clear all

eps = 5e-12;

x =[1; -1; 0; 1; 1];
xprim = zeros(5, 1);

A = [42 -50 -160 -4 378; 
    -44 46 154 20 -390; 
    -37 25 114 26 -297; 
    -43 25 120 38 -333; 
    -25 21 82 14 -209];

b = A*x;
deltaX = zeros(1, 1000000);
Aprim = zeros (5);

parfor n=1:1000000
    dist = (2*rand(5)-1)*eps;
    Aprim = A + dist;
    xprim = Aprim\b;
    deltaX(n) = (norm(xprim - x, 2))/(norm(x, 2));
end

deltaXmax = max(deltaX);
display(deltaXmax);