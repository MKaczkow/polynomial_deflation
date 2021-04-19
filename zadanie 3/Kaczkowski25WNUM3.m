clc
clear all

%---------------------------------------------------------------------
%-----------------------------PUNKT 1---------------------------------
%---------------------------------------------------------------------

ru = 8;
rl = -9;
coeff = [1, 9, -59, -1155, 1316, 44308, -162720];

%obliczenie ca³ki za pomoc¹ wyra¿eñ symbolicznych
syms x;
%y = symfun(x.^6 + 9.*x.^5 - 59.*x.^4 -1155.*x.^3 + 1316.*x.^2 + 44308.*x - 162720, [x]);
w = x.^6 + 9.*x.^5 - 59.*x.^4 -1155.*x.^3 + 1316.*x.^2 + 44308.*x - 162720;

%rzeczywiste pierwiastki wielomianu (uzyskane z zadania 2) to: -9 i 8
p_int = int(w, [rl ru]);
display(abs(p_int));
display(1004310547/420);

%---------------------------------------------------------------------
%-----------------------------PUNKT 2---------------------------------
%---------------------------------------------------------------------

DELTA = 5*10^-7;

% %liczba podprzedzialow
% subs_num = zeros(1,5);
% 
% %dlugosc kroku calkowania
% h = zeros(1,5);
% 
% for N=2:6
%    I=0;
%    subs=1;
%    
%    while(abs(I-p_int) > DELTA)
%        I=0;
%        %podzial przedzialu na podprzedzialy
%        x=(ru-rl)/subs;
%        
%        for i = 1:subs
%            a = rl + (i-1) * x;
%            b = rl + i * x;
%            %wykonanie kwadratury
%            I = I + kwadratura(coeff, a, b, N);
%        end
%        
%    subs=subs+3;
%    
%    end
%    
%    h(N-1) = x/N;
%    subs_num(N-1) = subs;
%      
%    
% end
% 
% %wykonanie wykresow
% figure(1); 
% semilogy((2:6),h, 'b-o'); 
% xlabel("N"); 
% ylabel("D³ugoœæ kroku ca³kowania"); 
% title('Kwadratura - d³ugoœæ kroku ca³kowania');
% grid on;
% 
% figure(2); 
% semilogy((2:6),subs_num, 'k-x') ; 
% xlabel("N"); 
% ylabel("liczba podprzedzia³ow"); 
% title('Kwadratura - liczba podprzedzia³ów');
% grid on;

%---------------------------------------------------------------------
%-----------------------------PUNKT 3---------------------------------
%---------------------------------------------------------------------

%UWAGA UWAGA UWAGA UWAGA UWAGA UWAGA UWAGA UWAGA UWAGA UWAGA
%wykonywanie kodu z punktu 3 zajmuje du¿o czasu, oko³o 20 min

N = [linspace(10,100,10), linspace(200,1000,9), linspace(2000,10000,9)];
c_eq = 0;
c_rand = 0;
estimate_eq = zeros(1,length(N));
estimate_rand =zeros(1, length(N));
delta_abs_equal = zeros(1, length(N));
delta_abs_random = zeros(1, length(N));

%wyznaczenie obszaru Omega:
%x nalezy do [-9, 8]
w_int = @(x) x.^7./7 + (3.*x.^6)./2 - (59.*x.^5)./5 - (1155.*x.^4)./4 + (1316.*x.^3)./3 + 22154.*x.^2 - 162720.*x;
num = linspace(rl, ru, 1000);
%oszacowanie zbioru wartosci funkcji w przedziale [-9,8]
MIN = min(w_int(num));
MAX = max(w_int(num));
%y nalezy do (MIN, MAX) oraz wiemy ze MIN < 0 i MAX > 0 zatem pole:
Omega = (abs(ru) + abs(rl))*(abs(MIN) + abs(MAX));

for k = 1:length(N)
    
%losowanie punktów

%wybor punktow
x_eq = linspace(-9, 8, N(k));
y_eq = (linspace(MIN, MAX, N(k))).';

%sprawdzenie czy punkty nale¿¹ do pola ograniczonego wykresem (eq)
for i = 1:N(k)
    for j = 1:N(k)          
        if (abs(w_int(x_eq(1,i))) > abs(y_eq(j,1)))
            c_eq=c_eq+1;
        end
    end
end

%wykonanie losowania N(k)^2 par punktów
for i = 1:((N(k))^2)
    x_rand = rand()*(abs(ru) + abs(rl)) + rl;
    y_rand = rand()*(abs(MIN) + abs(MAX)) + MIN;
    if (abs(w_int(x_rand))) >= abs(y_rand)
            c_rand=c_rand+1;
    end
end

%oszacowanie wartosci calki
estimate_rand(k) = (c_rand/(N(k)^2))*Omega;
delta_abs_random(k) = abs(p_int-estimate_rand(k));
estimate_eq(k) = (c_eq/(N(k)^2))*Omega;
delta_abs_equal(k) = abs(p_int-estimate_eq(k));
c_eq = 0;
c_rand = 0;

end
 
%wykonanie wykresu
figure(2)
loglog(N, delta_abs_random, 'r-x');
hold on;
loglog(N, delta_abs_equal, 'g-*');
hold off;
title('Monte-Carlo - modu³ b³êdu');
xlabel('Iloœæ punktów')
ylabel('Wartoœæ b³êdu bezwglêdnego');
legend('rozk³ad losowy', 'rozk³ad równomierny', 'Location', 'NorthEast');

figure(3)
semilogx(N, estimate_rand, 'r-x');
hold on;
semilogx(N, estimate_eq, 'g-*');
hold off;
title('Monte-Carlo - estymata pola');
xlabel('Iloœæ punktów')
ylabel('Estymata pola F');
legend('rozk³ad losowy', 'rozk³ad równomierny', 'Location', 'NorthEast');

%-------------------------------------------------------------------------
%-------------------------------FUNKCJE-----------------------------------
%-------------------------------------------------------------------------

function [wynik] = kwadratura(coeff, a, b, N)
H=(b-a)/N;
if (N==2) %kwadratura Simpsona
    wynik = 1/6*polyval(coeff,a) + 4/6*polyval(coeff,a+H) + 1/6*polyval(coeff,b);
end
if (N==3) %kwadratura trzech ósmych
    wynik =1/8*polyval(coeff,a) + 3/8*polyval(coeff,a+H) + 3/8*polyval(coeff,a+2*H) + 1/8*polyval(coeff,b);
end
if (N==4) %kwadratura Milne'a
    wynik =7/90*polyval(coeff,a) + 32/90*polyval(coeff,a+H) + 12/90*polyval(coeff,a+2*H) + 32/90*polyval(coeff,a+3*H) + 7/90*polyval(coeff,b);
end
if (N==5) %kwadratura Bode'a
    wynik =19/288*polyval(coeff,a) + 75/288*polyval(coeff,a+H) + 50/288*polyval(coeff,a+2*H) + 50/288*polyval(coeff,a+3*H) + 75/288*polyval(coeff,a+4*H) + 19/288*polyval(coeff,b);
end
if (N==6) %kwadratura Weddle'a
    wynik =41/840*polyval(coeff,a) + 216/840*polyval(coeff,a+H) + 27/840*polyval(coeff,a+2*H) + 272/840*polyval(coeff,a+3*H) + 27/840*polyval(coeff,a+4*H) + 216/840*polyval(coeff,a+5*H) + 41/840*polyval(coeff,b);
end
wynik=wynik*(b-a);
end

    