function [root] = styczne(start, Deltax, coeff)

x_next = start + 1;
x_i = start;
delta = abs(x_next - x_i);
count = 0;

while ((delta >= Deltax) && (count<1000))
    
    x_next = x_i - (wielomian(coeff,x_i))/(wielomian(pochodna(coeff),x_i));
    delta = abs(x_next - x_i);
    x_i = x_next;
    count = count+1;
    
end

root = x_i;

end