%% Quest�o 1: Andamento  da  pot�ncia  recebida,  em  fun��o  da  altura  da  antena  de recep��o

clearvars
close all
clc

% Raio equivalente da Terra para atmosfera padr�o (em metros)
re = 8500e3;
% Altura da antena emissora (em metros)
he = 80;
% Varia��o da altura da antena recetora (em metros)
hr = 0 : 0.25 : 150;
% Dist�ncia (em metros)
d = 50e3;
% Frequ�ncia (em hertz)
f = [1e9 5e9 10e9 20e9 60e9];

for i = 1 : length(f)
       
figure;
    
% Comprimento de onda (em metros)
lambda = 3e8 / f(i);

% Modelo: Terra Plana
Prp = 20 * log10(abs(2 * sin( (2*pi*he.*hr) / (lambda*d))));

subplot(1, 2, 1);
plot(hr, Prp);

title('Andamento da pot�ncia (Modelo: Terra Plana)');
xlabel('Altura da Antena de Recep��o (em metros)');
ylabel('Pot�ncia Recebida - Pot�ncia associada ao raio direto (em decib�is)');
legend( sprintf('%dGhz', f(i) / 1e9) );
grid on;

% Modelo: Terra Esf�rica, Sistema: Europeu

% C�lculo da localiza��o do ponto de reflex�o
xe = roots( [1 -1.5*d (0.5*(d^2) - re*(he+hr)) re*he*d] );
for j = 1 : length(xe)
    if xe(j) >= 0 && xe(j) <= d
        xe = xe(j);
        break;
    end
end

% Altura equivalente do emissor 
het = he - ( (xe^2) / (2 * re) );
% Altura equivalente do recetor
hrt = hr - ( ((d - xe)^2) / (2 * re) );

% C�lculo do �ngulo de incid�ncia
psi = atan(het / xe);

% C�lculo do fator de diverg�ncia
Dv = 1 / ( sqrt( 1 + ( (2*xe*(d-xe)) / (re*d*sin(psi)) ) ) );

% Modelo: Terra esf�rica
Pre = 10 * log10( abs( 1 + Dv^2 - 2 * Dv * cos( (4*pi*het.*hrt) / (lambda*d) ) ) );

subplot(1, 2, 2);
plot(hr, Pre);

title('Andamento da pot�ncia (Modelo: Terra Esf�rica, Sistema: Europeu)')
xlabel('Altura da Antena de Recep��o (em metros)');
ylabel('Pot�ncia Recebida - Pot�ncia associada ao raio direto (em decib�is)');
legend( sprintf('%dGhz', f(i) / 1e9) );

hold on;
grid on;

end

hold off;

%% Quest�o 2: Varia��o dos �ngulos de fogo de emiss�o e recep��o, em fun��o da dist�ncia da liga��o

clearvars
close all
clc

% Altura da antena emissora
he = 30;
% Altura da antena recetora
hr = 150;

% Varia��o da dist�ncia entre as antenas
dx = 0 : 0.125e3 : 50e3;

% Raio equivalente da Terra para atmosfera padr�o (em metros)
re = 8500e3;

%% Modelo: Terra Plana

% C�lculo do declive da equa��o da reta que descreve o raio direto
md = (hr - he) ./ dx;

% �ngulo de fogo na antena emissora
afe1 = abs(atand(md));

% �ngulo de fogo na antena recetora
afr1 = - afe1; 

%% Modelo: Terra Esf�rica, Sistema: Europeu

% C�lculo do declive da superf�cie da terra na antena emissora
dydxe = dx ./ (2 * re);
% C�lculo do declive da superf�cie da terra na antena recetora
dydxr = -dx ./ (2 * re);

% C�lculo do decaimento na antena emissora
ye = - ( ( (- dx/2).^2 ) / (2*re) );
% C�lculo do decaimento na antena recetora
yr = - ( (dx/2).^2 ) / (2*re);

% Altura da antena emissora (em rela��o � superf�cie terrestre)
heref = ye + he;
% Altura da antena recetora (em rela��o � superf�cie terrestre)
hrref = yr + hr;

% C�lculo do declive da equa��o da reta que descreve o raio direto
md = (hrref - heref) ./ dx;

% �ngulo de fogo na antena emissora
afe2 = ( abs(atand(dydxe)) - abs(atand(md)) );

% �ngulo de fogo na antena recetora
afr2 = - ( abs(atand(dydxr)) + abs(atand(md)) ); 

%% Gr�ficos

figure;
subplot(2, 1, 1);
plot(dx, afe1, dx, afe2);
grid on;
axis tight;

title('�ngulo de fogo no emissor')
xlabel('Dist�ncia de liga��o [em Km]');
ylabel('�ngulo de fogo [em graus]');
legend('Modelo: Terra Plana', 'Modelo: Terra Esf�rica, Sistema: Europeu');

subplot(2, 1, 2);
plot(dx, afr1, dx, afr2);
grid on;
axis tight;

title('�ngulo de fogo no recetor')
xlabel('Dist�ncia de liga��o [em Km]');
ylabel('�ngulo de fogo [em graus]');
legend('Modelo: Terra Plana', 'Modelo: Terra Esf�rica, Sistema: Europeu');

%% Quest�o 3: Compara��o Terra Plana - Terra Esf�rica e An�lise da Diferen�a de fase

clearvars
close all
clc

% Altura da antena emissora (em metros)
he = 150;
% Altura da altura da antena recetora (em metros)
hr = 150;
% Varia��o da dist�ncia (em metros)
dx = 0 : 50e3;
% Raio equivalente da Terra para atmosfera padr�o (em metros)
re = 8500e3;

% Frequ�ncia (em hertz)
f = 100e6;
% Comprimento de onda (em metros)
lambda = 3e8 / f;

for i = 1 : length(dx)
    d = dx(i);
    
    % C�lculo da localiza��o do ponto de reflex�o
    xe = roots( [1 -1.5.*d (0.5*(d.^2) - re*(he+hr)) re*he.*d] );
    for j = 1 : length(xe)
        if xe(j) >= 0 && xe(j) <= d
            xe = xe(j);
            break;
        end
    end

    % Altura equivalente do emissor 
    het = he - ( (xe^2) / (2 * re) );
    % Altura equivalente do recetor
    hrt = hr - ( ((d - xe)^2) / (2 * re) );
    
    deltaRte = (2.*het.*hrt)./d;
    deltaRtp = (2*he*hr)./d;
    
    phite = (-2*pi/lambda) * deltaRte;
    phitp = (-2*pi/lambda) * deltaRtp;
    
    dphi(i) = phitp - phite;
end

figure;
plot(dx, abs(dphi));
grid on;
title('Compara��o Terra Plana - Terra Esf�rica');
xlabel('Dist�ncia de liga��o (em metros)');
ylabel('Diferen�a de fase (em m�dulo e em graus)');

%% Quest�o 4

clearvars
close all
clc

% Frequ�ncia
fGhz = 20;
% Dist�ncia
d = 50;
% Press�o atmosf�rica (em hectopascais)
P = 1013;

%% Variando a Temperatura

% Temperatura (em graus cent�grados)
T1 = 0:100;

% Humidade relativa (em percentagem)
H1 = 0.80;

% C�lculo da atenua��o provocada pelo vapor de �gua

% C�lculo da press�o relativa � atmosfera padr�o (em atmosferas)
% 1013: Press�o m�dia da Atmosfera pr�xima da superf�cie terrestre (em hectopascais)
rp1 = P/1013;

% C�lculo da temperatura relativa � atmosfera padr�o
% 288: Temperatura m�dia da Atmosfera (em kelvins) 
rt1 = 288./(273 + T1);

% C�lculo da press�o parcial do vapor de �gua saturado (em hectopascais)
ppvas1 = 6.1121.*exp(17.502.*T1./(240.97+T1));

% C�lculo da press�o parcial do vapor de �gua (em hectopascais)
ppva1 = H1.*ppvas1;

% C�lculo da densidade do vapor de �gua
rho1 = 216.7 * (ppva1./(T1+273.3));

% C�clulo do coeficiente de atenua��o provocado pelo vapor de �gua (em dB/km)
gammaw1 = ((3.27e-2).*rt1 + (1.67e-3) .* rho1.*((rt1.^7)./rp1) + (7.7e-4) .* (fGhz.^0.5) + (3.79./(((fGhz-22.235).^2) + 9.81.*(rp1.^2) .* rt1))+((11.73 .* rt1)./(((fGhz-183.31).^2) + 11.85.*(rp1.^2).*rt1))+((4.01.*rt1)./(((fGhz-325.153).^2) + 10.44.*(rp1.^2) .* rt1)))*(fGhz.^2).*rho1.*rp1.*rt1.*1e-4; 

% C�lculo da atenua��o provocada pelo oxig�nio

% C�lculo do coeficiente de atenua��o provocado pelo oxig�nio (em dB/km)

% Para freq. menores ou iguais a 57GHz
if fGhz <= 57
    gammao1 = (((7.27.*rt1)./((fGhz.^2)+0.351.*(rp1.^2).*(rt1.^2))) + (7.5./(((fGhz - 57).^2)+2.44.*(rp1.^2).*(rt1.^5)))).*(fGhz.^2).*(rp1.^2).*(rt1.^2).*1e-3;
% Para freq. superiores a 57GHz e inferiores a 63GHz
elseif fGhz > 57 && fGhz < 63
    gammao1 = (((fGhz-60).*(fGhz-63))./18).*((((7.27.*rt1)./((57.^2)+0.351.*(rp1.^2).*(rt1.^2))) + (7.5./(((57 - 57).^2)+2.44.*(rp1.^2).*(rt1.^5)))).*(57.^2).*(rp1.^2).*(rt1.^2).*1e-3)-1.66.*(rp1.^2).*(rt1.^8.5).*(fGhz-57).*(fGhz-63) + (((fGhz-57).*(fGhz-60))./18).*((2.*1e-4.*(rt1.^1.5).*(1-1.2e-5.*(63.^1.5))+(4./(((63- 63).^2)+1.5.*(rp1.^2).*(rt1.^5)))+((0.28.*(rt1.^2))./(((63 - 118.75).^2)+2.84.*(rp1.^2).*(rt1.^2)))).*(63.^2).*(rp1.^2).*(rt1.^2).*1e-3);
% Para freq. superiores ou iguais a 350 GHz
elseif fGhz >= 350 
    gammao1 = (2.*1e-4.*(rt1.^1.5).*(1-1.2e-5.*(fGhz.^1.5))+(4./(((fGhz- 63).^2)+1.5.*(rp1.^2).*(rt1.^5)))+((0.28.*(rt1.^2))./(((fGhz - 118.75).^2)+2.84.*(rp1.^2).*(rt1.^2)))).*(fGhz.^2).*(rp1.^2).*(rt1.^2).*1e-3;
end

% C�lculo da atenua��o suplementar devida � presen�a da atmosfera (em dB)
Atg1 = (gammaw1 + gammao1).*d;

%% Variando a humidade relativa

% Temperatura (em graus cent�grados)
T2 = 25;

% Humidade relativa (em percentagem)
H2 = 0:0.01:1;

% C�lculo da atenua��o provocada pelo vapor de �gua

% C�lculo da press�o relativa � atmosfera padr�o (em atmosferas)
% 1013: Press�o m�dia da Atmosfera pr�xima da superf�cie terrestre (em hectopascais)
rp2 = P/1013;

% C�lculo da temperatura relativa � atmosfera padr�o
% 288: Temperatura m�dia da Atmosfera (em kelvins) 
rt2 = 288/(273 + T2);

% C�lculo da press�o parcial do vapor de �gua saturado (em hectopascais)
ppvas2 = 6.1121*exp(17.502*T2/(240.97+T2));

% C�lculo da press�o parcial do vapor de �gua (em hectopascais)
ppva2 = H2.*ppvas2;

% C�lculo da densidade do vapor de �gua
rho2 = 216.7 .* (ppva2./(T2+273.3));

% C�clulo do coeficiente de atenua��o provocado pelo vapor de �gua (em dB/km)
gammaw2 = ((3.27e-2).*rt2 + (1.67e-3) .* rho2.*((rt2.^7)./rp2) + (7.7e-4) .* (fGhz.^0.5) + (3.79./(((fGhz-22.235).^2) + 9.81.*(rp2.^2) .* rt2))+((11.73 .* rt2)./(((fGhz-183.31).^2) + 11.85.*(rp2.^2).*rt2))+((4.01.*rt2)./(((fGhz-325.153).^2) + 10.44.*(rp2.^2) .* rt2)))*(fGhz.^2).*rho2.*rp2.*rt2.*1e-4; 

% C�lculo da atenua��o provocada pelo oxig�nio

% C�clulo do coeficiente de atenua��o provocado pelo oxig�nio (em dB/km)

% Para freq. menores ou iguais a 57GHz
if fGhz <= 57
    gammao2 = (((7.27.*rt2)./((fGhz.^2)+0.351.*(rp2.^2).*(rt2.^2))) + (7.5./(((fGhz - 57).^2)+2.44.*(rp2.^2).*(rt2.^5)))).*(fGhz.^2).*(rp2.^2).*(rt2.^2).*1e-3;
% Para freq. superiores a 57GHz e inferiores a 63GHz
elseif fGhz > 57 && fGhz < 63
    gammao2 = (((fGhz-60).*(fGhz-63))./18).*((((7.27.*rt2)./((57.^2)+0.351.*(rp2.^2).*(rt2.^2))) + (7.5./(((57 - 57).^2)+2.44.*(rp2.^2).*(rt2.^5)))).*(57.^2).*(rp2.^2).*(rt2.^2).*1e-3)-1.66.*(rp2.^2).*(rt2.^8.5).*(fGhz-57).*(fGhz-63) + (((fGhz-57).*(fGhz-60))./18).*((2.*1e-4.*(rt2.^1.5).*(1-1.2e-5.*(63.^1.5))+(4./(((63- 63).^2)+1.5.*(rp2.^2).*(rt2.^5)))+((0.28.*(rt2.^2))./(((63 - 118.75).^2)+2.84.*(rp2.^2).*(rt2.^2)))).*(63.^2).*(rp2.^2).*(rt2.^2).*1e-3);
% Para freq. superiores ou iguais a 350 GHz
elseif fGhz >= 350 
    gammao2 = ((2.*1e-4.*(rt2.^1.5).*(1-1.2e-5.*(fGhz.^1.5))+(4./(((fGhz- 63).^2)+1.5.*(rp2.^2).*(rt2.^5)))+((0.28.*(rt2.^2))./(((fGhz - 118.75).^2)+2.84.*(rp2.^2).*(rt2.^2)))).*(fGhz.^2).*(rp2.^2).*(rt2.^2).*1e-3);
end

% C�lculo da atenua��o suplementar devida � presen�a da atmosfera (em dB)
Atg2 = (gammaw2 + gammao2).*d;

%% Gr�ficos
figure;
plot(T1, Atg1);
grid on;
axis tight;
title('Atenua��o provocada pelo vapor de �gua e pelo oxig�nio');
xlabel('Temperatura [�C]')
ylabel('Atenua��o [dB/Km]')

figure;
plot(H2, Atg2);
grid on;
axis tight;
title('Atenua��o provocada pelo vapor de �gua e pelo oxig�nio');
xlabel('Humidade Relativa [%]')
ylabel('Atenua��o [dB/Km]')

%% Quest�o 5: Representa��oda  atenua��o  provocada  pela  chuva,  n�o  excedida  numa  determinada percentagem do tempo, relativamente � m�dia anual ou ao pior m�s

clearvars
close all
clc

% Percentagem de tempo, relativamente � m�dia anual, em que o valor da intensidade da precipita��o � excedido (em percentagem)
p = 0.001:0.001:1;
% Frequ�ncia
fGhz = 20;
% Dist�ncia
d = 50;
% Regi�o (1 - H, 2 - K)
region = 1;
% Polariza��o (1 - Horizontal, 2 - Vertical)
polarizacao = 1;

% Tabela 2: k e alpha

% 1� linha: 1 gigahertz
% 5� linha: 5 gigahertz
% 10� linha: 10 gigahertz
% 20� linha: 20 gigahertz
% 60� linha: 60 gigahertz

% 1� coluna: kH
% 2� coluna: kV
% 3� coluna: alphaH
% 4� coluna: alphaV

kAlphaTable(1, 1:4) = [0.0000387 0.0000352 0.912 0.880];
kAlphaTable(4, 1:4) = [0.000650 0.000591 1.121 1.075];
kAlphaTable(6, 1:4) = [0.00175 0.00155 1.308 1.265];
kAlphaTable(10, 1:4) = [0.0101 0.00887 1.276 1.264];
kAlphaTable(20, 1:4) = [0.0751 0.0691 1.099 1.065];
kAlphaTable(60, 1:4) = [0.707 0.642 0.826 0.824];

% C�lculo por interpola��o do valor de k e alpha para a frequ�ncia de 5Ghz

alphaH5Ghz = kAlphaTable(4, 3) + ( ( ( kAlphaTable(6, 3) - kAlphaTable(4, 3) ) / ( log10(6) - log10(4) ) ) * ( log10(5) - log10(4) ) );
alphaV5Ghz = kAlphaTable(4, 4) + ( ( ( kAlphaTable(6, 4) - kAlphaTable(4, 4) ) / ( log10(6) - log10(4) ) ) * ( log10(5) - log10(4) ) );

kH5Ghz = 10^( log10(kAlphaTable(4,1)) + ( ( ( log10(kAlphaTable(6,1)) - log10(kAlphaTable(4,1)) ) / ( log10(6) - log10(4) ) ) * ( log10(5) - log10(4) ) ) );
kV5Ghz = 10^( log10(kAlphaTable(4,2)) + ( ( ( log10(kAlphaTable(6,2)) - log10(kAlphaTable(4,2)) ) / ( log10(6) - log10(4) ) ) * ( log10(5) - log10(4) ) ) );

kAlphaTable(5, 1:4) = [kH5Ghz kV5Ghz alphaH5Ghz alphaV5Ghz];

% C�lculo da intensidade de precipita��o para uma percentagem de tempo de 0.01 (em mil�metros por hora)
if region == 1 % Regi�o: H
    R = 32;
elseif region == 2 % Regi�o: K
    R = 42;
end

% C�lculo do k e do alpha
if polarizacao == 1 % Polariza��o: Horizontal
    k = kAlphaTable(fGhz, 1); 
    a = kAlphaTable(fGhz, 3);
elseif polarizacao == 2 % Polariza��o: Vertical    
    k = kAlphaTable(fGhz, 2); 
    a = kAlphaTable(fGhz, 4);
end

% C�culo do coeficiente de atenua��o provocada pela precipita��o (em dB/km)
gammar = k .* R.^a;

% C�lculo do comprimento eficaz do percurso 
def = d / (1 + (d / 35 *( exp(-0.015 * R) ) ) );

% C�lculo da atenua��o provocada pela precipita��o n�o excedida em 0.01% do tempo (relativamente � m�dia anual)
Ar1 = gammar * def;

% C�lculo da atenua��o provocada pela precipita��o n�o excedida numa determinada percentagem do tempo (relativamente � m�dia anual)
Arp = Ar1 .* 0.12 .* p.^(-(0.546 + 0.043 .* log10(p))); 

% C�lculo da atenua��o provocada pela precipita��o n�o excedida numa determinada percentagem do tempo (relativamente ao pior m�s)
pm = 0.3 .* p.^1.15;
Arpm = Ar1 .* 0.12 .* pm.^(-(0.546 + 0.043 .* log10(pm)));

figure;
plot(p, Arp, p, Arpm);
grid on;
axis tight;
title('Atenua��o provocada pela chuva');
xlabel('Percentagem de Tempo [%]')
ylabel('Atenua��o [dB]')
legend('relativamente � m�dia anual', 'relativamente ao pior m�s');

%% Quest�o 6: Assumindo uma varia��o linear do �ndice de refrac��o, efectue a representa��o do perfilda liga��o  em  sistema  Europeu a  partir  do  conceito  de  raio  equivalente

clearvars
close all
clc

% Temperatura (em graus cent�grados)
T = 18;
% Press�o atmosf�rica (em hectopascais)
P = 1017;
% C�lculo da press�o parcial do vapor de �gua (em hectopascais)
ppva = 10;
% Altura da antena de emiss�o
he = 30;
% Altura da antena de rece��o
hr = 150;
% Frequ�ncia (em Hz)
f = 10e9;
% Raio da terra
r0 = 6370;
% dndh
dndh = [-500e-6 -157e-6 -43e-6 50e-6];
% Dist�ncia (em metros)
d = 50000;

lambda = 3e8/f;

% C�lculo da refratividade
N = (77.6/(T+273))*(P + 4810*(ppva/(T+273)));

% C�lculo do �ndice de refrac��o
n = 1 + N * 1e-6;

% Varia��o do eixo das abcissas (em metros)
dx = 0 : d/1e6 : d;

for i = 1 : 4

dndhi = dndh(i);
    
% C�lculo do raio equivalente a partir do fator K e da varia��o do indice de refra��o
Kre = 1/(1+((r0/n)*dndhi));
req = Kre * r0 * 1e3;

% C�lculo do decaimento da terra para cada valor de dx
y = -( (dx - d/2).^2 )./(2*req);

% C�lculo das novas alturas relativamente ao decaimento da Terra(alturas equivalentes)
het = he + y(1);
hrt = hr + y(length(dx));

% C�lculo do declive da equa��o da reta que descreve o raio direto
md = (hrt - het) / d;

% C�lculo da ordenada na origem
bd = het;

% C�lculo da equa��o da reta que descreve o raio direto (yd = md*x + bd)
rd = md * dx + bd;

% Raio de cada sec��o do elips�ide para cada valor de dx (em metros)
r1 = sqrt( (dx .* (d-dx)) * lambda / d ); 

% Equa��es que descrevem o primeiro elips�ide de Fresnel
yup = rd + r1;
ydwn = rd - r1;

figure;

% Desenho da superf�cie da Terra
plot(dx, y, 'k');
hold on;
% Desenho do raio direto
plot(dx, rd, ['--' 'r']);
hold on;
% Desenho do primeiro elips�ide de Fresnel
plot(dx, ydwn, 'b');
plot(dx, yup, 'b');
grid on;
axis tight;
title('Perfil da liga��o em sistema Europeu');
ylabel('Altura [metros]')
xlabel('Dist�ncia [metros]')
legend('Superf�cie da Terra', 'Raio direto', 'Elips�ide de Fresnel');

end

%% Quest�o 7: Represente num gr�fico os v�rios valores do raio de curvatura do raio �ptico normalizado ao raio da Terra

clearvars
close all
clc

dndh = (-500:10:100).*10^-6;
r0 = 6370e3;
d = 50e3;
type = 0;
dg=0;
ctE = -(((dg-(d./2)).^2)./(2.*r0));

hEf = 100 + ctE;
phi = atan(d/(2*r0));
x = 0:1:d;

if (type == 0)% Modelo Europeu
    rS = sqrt(r0^2 - (x - d/2).^2)-r0;
elseif (type == 1)% Modelo Americano
    rS = sqrt(r0^2 - x.^2)-r0;
end

figure;
plot(x, rS, '-k');
hold on
plot([0 0], [ctE hEf], '-g', 'LineWidth', 2);
hmax = get(gca, 'ylim');
hold on

for i=1:length(dndh)
    % Raio �tico
    n = 1+ 315*10^-6*exp(-0.136*hEf*10^-3); %Pelo slide 28 (Influ�ncia da Atmosfera)
    p = -n./(dndh(i).*sin(phi)); %Pelo slide 37 (Influ�ncia da Atmosfera)
    x = 0:1:d;
    y = tan(phi).*x + hEf;

    if (dndh(i) == 0 || p == Inf)
            h = plot(x, phi+hEf,'-y');
            hold on
    end
    
    if (dndh(i) > 0) 
        ro = (-sqrt(p^2 - x.^2)- p) + y;
        h = plot(x, ro); '-y';
        hold on
    elseif (dndh(i) < 0)
        ro = -((-sqrt(p^2 - x.^2) + p)) + y;
        h = plot(x, ro); '-y';
        hold on
    end
end

xlabel('Dist�ncia [em metros]'); 
ylabel('Altura [em metros]');
title('Raios �ticos');
grid on;

%% Quest�o 8: Representa��o da diferen�a de atenua��o de propaga��o dada pelos m�todos de c�lculo das perdas por difra��o em Terra Esf�rica

clearvars
close all
clc

% Frequ�ncia (em GHz)
fGhz = 10;
% Polariza��o (1 - Horizontal, 2 - Vertical)
polarizacao = 1;
% Altura da antena de emiss�o
he = 100;
% Altura da antena de rece��o
hr= 100;
% Frequ�ncia (em Hz)
f = 10e9;
% Dist�ncia
d = 100e3:50e3:1000e3;
% Raio equivalente da Terra para atmosfera padr�o (em metros)
re = 6370*1e3*(4/3);
% Tipo de Terreno
terrain = 3;

lambda = 3e8/f;

% Tabela 1: Tipo de terreno

% 1� linha: Mar
% 2� linha: Muito H�mido
% 3� linha: Normal
% 4� linha: Seco

% 1� coluna: Permitividade relativa para cada tipo de terreno
% 2� coluna: Condutividade el�trica para cada tipo de terreno (em siemens por metro)

terrainTable = [
    81 5 ; 
    25 0.02 ;
    15 0.005 ;
    4 0.001 
];

y = -( (d - d/2).^2 )/(2*re);

% Altura equivalente da antena de emiss�o
het = he + y(1);
% Altura equivalente da antena de rece��o
hrt = hr + y(length(d));
% Efeito do Solo e Polariza��o
beta = 1;

% Fun��o x
x = beta * d *((pi/(lambda*(re^2)))^(1/3));
% Fun��o F(x)
Fx = 11 + 10 * log10(x) - 17.6 * x;

% Fun��o y1
y1 = 2 * beta * het *(((pi^2)/((lambda^2)*re))^(1/3));
% Fun��o y2
y2 = 2 * beta * hrt *(((pi^2)/((lambda^2)*re))^(1/3));

% C�lculo do K
KH = (((2*pi*re)/lambda)^(-1/3))*(((terrainTable(terrain, 1)-1)^2)+(60*lambda*terrainTable(terrain, 2))^2)^(-1/4);
if polarizacao == 0
    K = KH;
else
    K = KH*sqrt((terrainTable(terrain, 1)^2)+(60*lambda*terrainTable(terrain, 2))^2);
end

% Fun��o G(y2)
if (y2 < 2)
    if 10*K < y2
        Gy2 = 20*log10(y2 + 0.1 * (y2^3));
    elseif 10/K < y2 <10*K
        Gy2 = 2 + 20*log10(K) + 9*log10(y2/K)*(log10(y2/K)+1);
    elseif y2 < 10/K
        Gy2 = 2 + 20*log10(K);
    end
else
    Gy2 = 17.6 * sqrt(y2 - 1.1) - 5 * log10(y2 - 1.1)-8;
end

% Fun��o G(y1)
if (y1 < 2)
    if 10*K < y1
        Gy1 = 20*log10(y1 + 0.1 * (y1^3));
    elseif 10/K < y1 <10*K
        Gy1 = 2 + 20*log10(K) + 9*log10(y1/K)*(log10(y1/K)+1);
    elseif y1 < 10/K
        Gy1 = 2 + 20*log10(K);
    end
else
    Gy1 = 17.6 * sqrt(y1 - 1.1) - 5 * log10(y1 - 1.1)-8;
end

% Difrac��o
Ad = -(Fx + Gy1 + Gy2);

figure;
plot(d, Ad);
grid on;
axis tight;
ylabel('Atenua��o [dB]');
xlabel('Dist�ncia [metros]');
title('Perdas por Difrac��o em Terra Esf�rica');

%% Quest�o 9: Leitura de um determinado perfil de terreno a partir de dados existentes num ficheiro Excel e adicione � representa��o da terra a representa��o desse mesmo perfil

clearvars
close all
clc

% Dist�ncia (em Km)
dKm = 45;
% Dist�ncia (em metros)
d = dKm*1e3;
% Raio equivalente da Terra para atmosfera padr�o (em metros)
re = 8500e3;

perfilTerreno = (xlsread('terrainProfile.xlsx'))';

if length(perfilTerreno) ~= dKm
    error('N�mero de pontos definidos no ficheiro excel n�o pode ser diferente da dist�ncia (em quil�metros)')
end

% Varia��o do eixo das abcissas (em metros)
dx = 0 : d/(length(perfilTerreno)-1) : d;

% C�lculo do decaimento da terra para cada valor de dx
y = -( (dx - d/2).^2 )./(2*re);

% Gr�fico
figure;
plot(dx/1e3, perfilTerreno+y);
grid on;
axis tight;
title('Representa��o da terra e do perfil de terreno');
xlabel('Dist�ncia [em Km]');
ylabel('Altura [em metros]');

%% Quest�o 10: C�lculo da atenua��o introduzida por difrac��o

% Frequ�ncia (em Hz)
f = 13e9;
% Altura da antena de emiss�o
he = 100;
% Altura da antena de rece��o
hr = 100;
% Frequ�ncia (em megahertz)
fMhz = f / 1e6;

lambda = 3e8/f;

% Atenua��o em espa�o livre (em decib�is)
A0 = 32.4 + 20*log10(dKm) + 20*log10(fMhz);

% C�lculo das novas alturas relativamente ao decaimento da Terra (alturas equivalentes)
het = he + y(1);
hrt = hr + y(length(dx));

AtenuacaoObstaculos = Deygout(d, het, hrt, lambda, dKm, perfilTerreno, re, 0, 0) + A0;
