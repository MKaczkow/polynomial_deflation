function [root] = muller_II(start, Deltax, coeff)

%w = @(x) x.^2 - 8*x + 20;

x_i = start;
x_prev_prev = start + 1;
x_prev = start + 0.5;
x_next = start + 1;

delta = abs(x_next - x_i);
count = 0;

while ((delta > Deltax) && (count<1000))
    
    M = [-(x_i-x_prev)^2, x_i-x_prev; -(x_i-x_prev_prev)^2, x_i-x_prev_prev];
    v = [wielomian(coeff,x_i) - wielomian(coeff,x_prev); wielomian(coeff,x_i) - wielomian(coeff,x_prev_prev)];
    
    if det(M)==0
        break;
    end
    M = inv(M);
    wektor = M*v;
    
    a = wektor(1);
    b = wektor(2);
    c = wielomian(coeff,x_i);
    x_next = x_i - (2*c)/(b + sign(b)*sqrt(b^2 - 4*a*c));
    delta = abs(x_next - x_i);
    x_prev_prev = x_prev;
    x_prev = x_i;
    x_i = x_next;
    count = count+1;
    
end

root = x_i;

end