clc
clear all

coefficients0 = [1, 9, -59, -1155, 1316, 44308, -162720];
p_roots = transpose(roots(coefficients0)); %pierwiastki uzyskane z roots, w kolejnoœci sa to -8+7j, -8-7j, -9, 8, 4+2j, 4-2j
c_roots = zeros(1,6); %miejsce na pierwiastki obliczone metodami iteracyjnymi

%przyblizenie calkowite wartosci p_roots(przy u¿yciu format long widaæ, ¿e
%procedura roots nie daje dok³adnie wartoœci ca³kowitych pierwiastków
for n = 1:6

    if isreal(p_roots(n))
        p_roots(n) = round(p_roots(n));
    else
        a = round(real(p_roots(n)));
        b = round(imag(p_roots(n)));
        p_roots(n) = a + b*1i;
    end  
end

%---------------------------------------------------------------------
%-----------------------------PUNKT 1---------------------------------
%---------------------------------------------------------------------

Deltax = 10^-3;

% krok 1 - metoda siecznych
c_roots(4) = sieczne(10, Deltax, coefficients0);

%krok 2 - deflacja liniowa za pomoca algorymtu Hornera
coefficients1 = deflacja_lin(c_roots(4), coefficients0);

%krok 3 - metoda stycznych (Newtona)
c_roots(3) = styczne(-10, Deltax, coefficients1);

%krok 4 - deflacja liniowa za pomoca algorymtu Hornera

coefficients2 = deflacja_lin(c_roots(3), coefficients1);

%krok 5 - metoda Mullera w wersji I

c_roots(2) = muller_I(-10, Deltax, coefficients2);
c_roots(1) = conj(c_roots(2));

%kork 6 - deflacja kwadratowa za pomoca algorytmu Hornera

coefficients3 = deflacja_kw(c_roots(2), coefficients2);

%krok 7 - metoda Mullera w wersji II

c_roots(5) = muller_II(10, Deltax, coefficients3);
c_roots(6) = conj(c_roots(5));

%---------------------------------------------------------------------
%-----------------------------PUNKT 2---------------------------------
%---------------------------------------------------------------------

Deltax = zeros(1, 13);
Agreg_sqr = zeros(1, 13);
Agreg_inf = zeros(1,13);

for n = 1:14
    
    Deltax(n) = 10^-(17-n);
    
    c_roots(4) = sieczne(10, Deltax(n), coefficients0);
    
    coefficients1 = deflacja_lin(c_roots(4), coefficients0);
    c_roots(3) = sieczne(-10, Deltax(n), coefficients1);
    
    coefficients2 = deflacja_lin(c_roots(3), coefficients1);
    c_roots(2) = muller_I(-10, Deltax(n), coefficients2);
    c_roots(1) = conj(c_roots(2));
    
    coefficients3 = deflacja_kw(c_roots(2), coefficients2);
    c_roots(5) = muller_II(10, Deltax(n),coefficients3);
    c_roots(6) = conj(c_roots(5));
    
    Agreg_sqr(n) = norm(c_roots-p_roots, 2)/norm(p_roots,2);
    Agreg_inf(n) = norm(c_roots-p_roots, Inf)/norm(p_roots,Inf);
    
end

figure(1);
loglog(Deltax, Agreg_sqr, 'b:*');
title('Zagregowany blad sredniokwadratowy');

figure(2);
loglog(Deltax, Agreg_inf, 'r:*');
title('Zagregowany blad maksymalny');



%funkcje pomocnicze

function [root] = sieczne(start, Deltax, coeff)

x_next = start + 1;
x_i = start;
x_prev = start + 1;
delta = abs(x_next - x_i);
count = 0;

while ((delta >= Deltax) && (count<1000))
    
    x_next = x_i - ((x_i - x_prev)/(polyval(coeff,x_i)-polyval(coeff,x_prev)))*polyval(coeff,x_i);
    delta = abs(x_next - x_i);
    x_prev = x_i;
    x_i = x_next;
    count = count+1;
    
end

root = x_i;

end


function [root] = styczne(start, Deltax, coeff)

x_next = start + 1;
x_i = start;
delta = abs(x_next - x_i);
count = 0;

while ((delta >= Deltax) && (count<1000))
    
    x_next = x_i - (polyval(coeff,x_i))/(polyval(polyder(coeff),x_i));
    delta = abs(x_next - x_i);
    x_i = x_next;
    count = count+1;
    
end

root = x_i;

end


function [root] = muller_I(spoint, Deltax, coeff)

x_next = 0;
x_i = spoint;
delta = abs(x_next - x_i);
count = 0;

while ((delta >= Deltax) && (count<1000))
    
    a = (1/2)*polyval(polyder(polyder(coeff)),x_i);
    b = polyval(polyder(coeff),x_i);
    c = polyval(coeff,x_i);
    x_next = x_i - (2*c)/(b + sign(b)*sqrt(b^2 - 4*a*c));
    delta = abs(x_next - x_i);
    x_i = x_next;
    count = count+1;
    
end

root = x_i;

end


function [root] = muller_II(start, Deltax, coeff)

x_i = start;
x_prev_prev = start + 1;
x_prev = start + 0.5;
x_next = start + 1;

delta = abs(x_next - x_i);
count = 0;

while ((delta > Deltax) && (count<1000))
    
    M = [-(x_i-x_prev)^2, x_i-x_prev; -(x_i-x_prev_prev)^2, x_i-x_prev_prev];
    v = [polyval(coeff,x_i) - polyval(coeff,x_prev); polyval(coeff,x_i) - polyval(coeff,x_prev_prev)];
    
    if det(M)==0
        break;
    end
    
    %zamiast funkcji Inv() u¿yto dzielenia macierzy za pomoc¹ operatora '\'
    wektor = M\v;
    
    a = wektor(1);
    b = wektor(2);
    c = polyval(coeff,x_i);
    x_next = x_i - (2*c)/(b + sign(b)*sqrt(b^2 - 4*a*c));
    delta = abs(x_next - x_i);
    x_prev_prev = x_prev;
    x_prev = x_i;
    x_i = x_next;
    count = count+1;
    
end

root = x_i;

end


function [deflated] = deflacja_lin(root, coeff)

deflated = zeros(1,length(coeff)-1);
deflated(1) = coeff(1);

for n = 2 : length(coeff)-1
    deflated(n) = coeff(n) + (deflated(n-1)*root);
end

end


function [deflated] = deflacja_kw(root, coeff)

deflated = zeros(1,length(coeff)-2);
p = 2*real(root);
r = -(abs(root))^2;

deflated(1) = coeff(1);
deflated(2) = coeff(2) + deflated(1)*p;

for n = 3 : length(coeff)-2
    deflated(n) = coeff(n) + deflated(n-1)*p + deflated(n-2)*r;
end    
   
end


