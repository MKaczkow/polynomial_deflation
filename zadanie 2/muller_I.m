function [root] = muller_I(spoint, Deltax, coeff)

x_next = 0;
x_i = spoint;
delta = abs(x_next - x_i);
count = 0;

while ((delta >= Deltax) && (count<1000))
    
    a = (1/2)*wielomian(pochodna(pochodna(coeff)),x_i);
    b = wielomian(pochodna(coeff),x_i);
    c = wielomian(coeff,x_i);
    x_next = x_i - (2*c)/(b + sign(b)*sqrt(b^2 - 4*a*c));
    delta = abs(x_next - x_i);
    x_i = x_next;
    count = count+1;
    
end

root = x_i;

end