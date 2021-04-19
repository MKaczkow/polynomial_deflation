clc
clear all

%deklaracja zmiennych i niektórych funkcji
syms x v w;
eps = 5e-12;
y = cos(x.^2+2).*exp(x.^3+2);
p = linspace(0,1,1000);

%---------------------------------------
%------------PODPUNKT 1-----------------
%---------------------------------------

%wspolczynnik transmisji obliczony analitycznie
Tx = x/y*diff(y,x);

%wspolczynnik transmisji obliczony za pomoca rachunku epsilonow "na kartce"
Tx2 = (3*x.^3)-((2*x.^2)*tan(x.^2+2));

figure(1);
title('T(x)');
fplot (Tx, [0,1], '-g');
hold on;
fplot (Tx2, [0,1], 'x r');
hold off;

fTx = @(x) (3*x.^3)-((2*x.^2).*tan(x.^2+2));   %maks w 1

%pozostale wspolczynniki obliczone za pomoca rachunku epsilonow "na kartce"
Kcos = 1;
Kexp = 1;

Km2 = @(x) 1+x.^3;                        %maks w 1
Ks2 = @(x) x.^3+2-(x.^2+2).*tan(x.^2+2);  %maks w 0
Kp2 = @(x) x.^3 - (x.^2).*tan(x.^2+2);    %maks w 1


%pozostale wspolczynniki obliczone analitycznie

%sumowanie
y = cos(x.^2+2).*exp(x.^3+2);
ys = subs(y, x^2+2, v);
Ksprim = v/ys*diff(ys,v);
Ksprim = subs(Ksprim, v, x^2+2);
ys = subs(y, x^3+2, w);
Ksbis = w/ys*diff(ys,w);
Ksbis = subs(Ksbis, w, x^3+2);
Ks = Ksprim + Ksbis;

figure(2);
title('Ks(x)');
fplot (Ks, [0,1], '-k');
hold on;
fplot (Ks2, [0,1], 'x r');
hold off;


%---------------------------------------
%------------PODPUNKT 2-----------------
%---------------------------------------

delta1 = (fTx(p)+Kcos+Kexp+Km2(p)+Ks2(p)+Kp2(p)).*eps;
dMAX1 = max(abs(delta1));
display (dMAX1);


%---------------------------------------
%------------PODPUNKT 3-----------------
%---------------------------------------

epsy = de2bi(0:4096);
epsy(epsy==0) = -1;
epsy = epsy*eps;

delta2 = zeros(1, 4096);
y = cos(p.^2+2).*exp(p.^3+2);

parfor n = 1:4096
    ep = epsy(n,:);
    yd = cos((((p.*(1+ep(1))).^2).*(1+ep(2))+2).*(1+ep(3))).*(1+ep(4)).*(1+ep(5)).*exp((((p.*(1+ep(6))).^2).*(p.*(1+ep(7))).*(1+ep(8)).*(1+ep(9))+2).*(1+ep(10)))*(1+ep(11));
    delta2(1, n) = abs((yd-y)/y);
end

dMAX2 = max(delta2);
display (dMAX2);

%---------------------------------------
%------------PODPUNKT 4-----------------
%---------------------------------------

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

dXMAX = max(deltaX);
display(dXMAX);
