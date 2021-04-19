function [value] = wielomian(coeff, point)

%funkcja oblicza wartosc danego wielomianu (za pomoca wektora wspolczynnikow) w danym punkcie 

value = 0;
degree = length(coeff)-1;
co_degree = degree;

for n = 0:co_degree

    value = value + coeff(n+1).*(point^degree);
    degree = degree-1;
    
end
    
end