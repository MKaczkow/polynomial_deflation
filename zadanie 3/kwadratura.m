function [wynik] = kwadratura(coeff, a, b, N)
H=(b-a)/N;
if (N==2)
    wynik = 1/6*polyval(coeff,a) + 4/6*polyval(coeff,a+H) + 1/6*polyval(coeff,b);
end
if (N==3)
    wynik =1/8*polyval(coeff,a) + 3/8*polyval(coeff,a+H) + 3/8*polyval(coeff,a+2*H) + 1/8*polyval(coeff,b);
end
if (N==4)
    wynik =7/90*polyval(coeff,a) + 32/90*polyval(coeff,a+H) + 12/90*polyval(coeff,a+2*H) + 32/90*polyval(coeff,a+3*H) + 7/90*polyval(coeff,b);
end
if (N==5)
    wynik =19/288*polyval(coeff,a) + 75/288*polyval(coeff,a+H) + 50/288*polyval(coeff,a+2*H) + 50/288*polyval(coeff,a+3*H) + 75/288*polyval(coeff,a+4*H) + 19/288*polyval(coeff,b);
end
if (N==6)
    wynik =41/840*polyval(coeff,a) + 216/840*polyval(coeff,a+H) + 27/840*polyval(coeff,a+2*H) + 272/840*polyval(coeff,a+3*H) + 27/840*polyval(coeff,a+4*H) + 216/840*polyval(coeff,a+5*H) + 41/840*polyval(coeff,b);
end
wynik=wynik*(b-a);
end

    