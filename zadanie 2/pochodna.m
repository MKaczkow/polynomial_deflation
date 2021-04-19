function [derr] = pochodna(coeff)

degree = length(coeff)-1;
derr = zeros(1, degree-1);
co_degree = degree;

for n = 0:co_degree-1
    
    derr(n+1) = degree*coeff(n+1);
    degree = degree - 1;

end