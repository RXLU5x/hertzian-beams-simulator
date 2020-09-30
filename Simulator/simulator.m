%% Inicializa��o da aplica��o

clearvars
close all
clc

%% Defini��o das constantes

% Raio da terra (em metros)
r0 = 6370e3;

% Velocidade da luz no vazio (em metros por segundo)
c = 3e8;

% Permitividade el�trica do v�cuo (em farad por metro)
e0 = 8.854e-12;

% Ind�ce de refra��o a uma altura de 0km
n0 = 1 + 315e-6;

%% Defini��o das tabelas

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

% Tabela 3: 
% Ficheiro excel

perfilTerreno = (xlsread('terrainProfile.xlsx'))';

%% Configura��o do modelo (e sistema)

fprintf('CONFIGURA��O DO MODELO\n\n');

% Modelo
fprintf('Selecione o modelo que pretende simular:\n\n1: Terra plana\n2: Terra esf�rica\n\n');
model = abs(input('Modelo: '));
disp(' ');

switch model
    case 1
        fprintf('Selecionou como modelo: Terra plana\n\n');
    case 2
        fprintf('Selecionou como modelo: Terra esf�rica\n\n');
        
        fprintf('CONFIGURA��O DO SISTEMA\n\n');
        
        % Sistema
        fprintf('Selecione o sistema que pretende simular:\n\n1: Europeu\n2: Americano');
        system = abs(input('\n\nSistema: '));
        disp(' ');
        
        switch system
            case 1
                fprintf('Selecionou como sistema: Europeu\n\n');
            case 2
                fprintf('Selecionou como sistema: Americano\n\n');
            otherwise
                error('Sistema selecionado n�o existe');
        end
    otherwise
        error('Modelo selecionado n�o existe');
end

%% Configura��o das antenas

fprintf('CONFIGURA��O DAS ANTENAS\n\n');

% Polariza��o
fprintf('Selecione a polariza��o:\n\n1: Horizontal\n2: Vertical\n\n');
polarizacao = abs(input('Polariza��o: '));
disp(' ');

switch polarizacao
    case 1
        fprintf('Selecionou como polariza��o: Horizontal\n\n');
    case 2
        fprintf('Selecionou como polariza��o: Vertical\n\n');
    otherwise
        error('Polariza��o selecionada n�o existe');
end

% Dist�ncia entre as antenas (em metros)
d = abs(input('Dist�ncia entre as antenas (em metros): '));
disp(' ');

%% Configura��o da antena emissora

fprintf('CONFIGURA��O DA ANTENA EMISSORA\n\n');

% Frequ�ncia (em hertz)
% 1: 1e9 Hz
% 2: 5e9 Hz
% 3: 10e9 Hz
% 4: 20e9 Hz
% 5: 60e9 Hz
fprintf('Selecione a frequ�ncia:\n\n1: 1Ghz\n2: 5Ghz\n3: 10Ghz\n4: 20Ghz\n5: 60Ghz\n\n');
f = abs(input('Frequ�ncia: '));
disp(' ');

switch f
    case 1
        fprintf('Selecionou como frequ�ncia: 1Ghz\n\n');
        f = 1e9;
    case 2
        fprintf('Selecionou como frequ�ncia: 5Ghz\n\n');
        f = 5e9;
    case 3
        fprintf('Selecionou como frequ�ncia: 10Ghz\n\n');
        f = 10e9;
    case 4
        fprintf('Selecionou como frequ�ncia: 20Ghz\n\n');
        f = 20e9;
    case 5
        fprintf('Selecionou como frequ�ncia: 60Ghz\n\n');
        f = 60e9;
    otherwise
        error('Frequ�ncia inv�lida');
end

% Altura (em metros)
he = abs(input('Altura (em metros): '));

% Di�metro do prato (em metros)
le = abs(input('Di�metro do prato (em metros): ')); 

% Pot�ncia (em decib�is)
Pe = 10*log10(abs(input('Pot�ncia (em watts): ')));

% Rendimento (em percentagem)
ne = abs(input('Rendimento (em percentagem): '));

% Comprimento do guia de onda (em metros)
cge = abs(input('Comprimento do guia de onda (em metros): '));

% Atenua��o no guia de onda (em decib�is por quil�metro)
Aekm = abs(input('Atenua��o no guia de onda (em decib�is por quil�metro): '));

disp(' ');

%% Configura��o da antena recetora

fprintf('CONFIGURA��O DA ANTENA RECETORA\n\n');

% Altura (em metros)
hr = abs(input('Altura (em metros): '));

% Di�metro do prato (em metros)
lr = abs(input('Di�metro do prato (em metros): ')); 

% Rendimento (em percentagem)
nr = abs(input('Rendimento (em percentagem): '));

% Comprimento do guia de onda (em metros)
cgr = abs(input('Comprimento do guia de onda (em metros): '));

% Atenua��o no guia de onda (em decib�is por quil�metro)
Arkm = abs(input('Atenua��o no guia de onda (em decib�is por quil�metro): '));

disp(' ');

%% Configura��o do meio ambiente

fprintf('CONFIGURA��O DO MEIO AMBIENTE\n\n');

% Terreno:
% 1: Mar
% 2: Muito H�mido
% 3: Normal
% 4: Seco
fprintf('Selecione o tipo de terreno:\n\n1: Mar\n2: Muito H�mido\n3: Normal\n4: Seco\n\n');
terrain = abs(input('Terreno: '));
disp(' ');

switch terrain
    case 1
        fprintf('Selecionou como tipo de terreno: Mar\n\n');
    case 2
        fprintf('Selecionou como tipo de terreno: Muito H�mido\n\n');
    case 3
        fprintf('Selecionou como tipo de terreno: Normal\n\n');
    case 4
        fprintf('Selecionou como tipo de terreno: Seco\n\n');
    otherwise
        error('Terreno selecionado n�o existe');
end

% Regi�o:
% 1: H
% 2: K
fprintf('Selecione a regi�o:\n\n1: H\n2: K\n\n');
region = abs(input('Regi�o: '));
disp(' ');

switch region
    case 1
        fprintf('Selecionou como regi�o: H\n\n');
    case 2
        fprintf('Selecionou como regi�o: K\n\n');
    otherwise
        error('Regi�o selecionada n�o existe\n\n');
end

% Press�o atmosf�rica (em hectopascais)
P = abs(input('Press�o atmosf�rica (em hectopascais): '));

% Temperatura (em graus cent�grados)
T = abs(input('Temperatura (em graus cent�grados): '));

% Humidade relativa (em percentagem)
H = abs(input('Humidade relativa (em percentagem): '));

% Varia��o do �ndice de refra��o a uma altura de 0km (dn/dh) (em km^-1)
% Valor para atmosfera padr�o: -43e-6
dndh0 = abs(input('Varia��o do �ndice de refra��o (em km^-1): ')) * 1e-3;

% Percentagem de tempo, relativamente � m�dia anual, em que o valor da intensidade da precipita��o � excedido (em percentagem)
p = abs(input('Percentagem de tempo, relativamente � m�dia anual, em que o valor da intensidade da precipita��o � excedido (em percentagem): '));

disp(' ');

figure;
title('Simulador');
xlabel('Dist�ncia (m)');
ylabel('Altura (m)');
grid on;
hold on;

%% C�lculos auxiliares

% Dist�ncia (em quil�metros)
dKm = d/1e3;

% Frequ�ncia (em gigahertz)
fGhz = f/1e9;

% Frequ�ncia (em megahertz)
fMhz = f/1e6;

% Comprimento de onda (em metros)
lambda = c/f;

% Varia��o do eixo das abcissas (em metros)
dx = 0 : d/1e6 : d;

%% Desenho do referencial

if model == 1 % Modelo: Terra plana
    % C�lculo do decaimento da terra para cada valor de dx
    y = zeros(1, length(dx));
    
    % C�lculo do declive da reta tangente � superficie terrestre no emissor
    dydxe = 0;

    % C�lculo do declive da reta tangente � superficie terrestre no recetor
    dydxr = 0;
elseif model == 2 && system == 1 % Modelo: Terra esf�rica, Sistema: Europeu
    % C�lculo do raio equivalente a partir do fator K e da varia��o do indice de refra��o
    Kre = 1/(1+((r0/n0)*dndh0));
    r0 = Kre * r0;
    
    % C�lculo do decaimento da terra para cada valor de dx
    y = -( (dx - d/2).^2 )/(2*r0);
    
    % Altura equivalente do emissor 
    het = he + y(0);
    % Altura equivalente do recetor
    hrt = hr + y(length(dx));
    
    % C�lculo do declive da reta tangente � superficie terrestre no emissor
    dydxe = d/2 / r0;

    % C�lculo do declive da reta tangente � superficie terrestre no recetor
    dydxr = - (d - d/2) / r0;
    
elseif model == 2 && system == 2 % Modelo: Terra esf�rica, Sistema: Americano
    % C�lculo do raio equivalente a partir do fator K e da varia��o do indice de refra��o
    Kre = 1/(1+((r0/n0)*dndh0));
    r0 = Kre * r0;
    
    % C�lculo do decaimento da terra para cada valor de dx
    y = -(dx.^2)/(2*r0);

    % Altura equivalente do emissor 
    het = he + y(0);
    % Altura equivalente do recetor
    hrt = hr + y(length(dx));
    
    % C�lculo do declive da reta tangente � superficie terrestre no emissor
    dydxe = 0;

    % C�lculo do declive da reta tangente � superficie terrestre no recetor
    dydxr = -d / r0;
end

% Desenho da superf�cie da Terra
plot(dx, y, 'k');

%% Desenho da horizontal � superf�cie da Terra no emissor e no recetor

if model == 2
    % C�lculo da ordenada na origem
    bte = het;

    % C�lculo da equa��o da reta paralela � reta tangente em rela��o � superficie terrestre no emissor
    yte = dydxe .* dx + bte;

    % Desenho da horizontal � superf�cie da Terra no emissor
    plot(dx, yte, [':' 'k']);

    % C�lculo da ordenada na origem
    btr = hrt - dydxr .* d;

    % C�lculo da equa��o da reta paralela � reta tangente em rela��o � superficie terrestre no recetor
    ytr = dydxr .* dx + btr;

    % Desenho da horizontal � superf�cie da Terra no recetor
    plot(dx, ytr, [':' 'k']);
end

%% 1. Desenho do Raio Direto

fprintf('DESENHO DO RAIO DIRETO\n\n');
    
if model == 1 % Modelo: Terra plana
    % C�lculo do declive da equa��o da reta que descreve o raio direto
    md = (hr - he) / d;

    % C�lculo da ordenada na origem
    bd = he;
else % Modelo: Terra esf�rica
    % C�lculo do declive da equa��o da reta que descreve o raio direto
    md = (hrt - het) / d;

    % C�lculo da ordenada na origem
    bd = het;
end

% C�lculo da equa��o da reta que descreve o raio direto (yd = md*x + bd)
rd = md * dx + bd;

% Desenho do raio direto
plot(dx, rd, ['--' 'k']);
    
%% 2. Desenho do primeiro elips�ide de fresnel e indica��o se este se econtra obstru�do

fprintf('DESENHO DO PRIMEIRO ELIPS�IDE DE FRESNEL E INDICA��O SE ESTE SE ECONTRA OBSTRU�DO\n\n');

% Raio de cada sec��o do elips�ide para cada valor de dx (em metros)
r1 = sqrt( (dx .* (d-dx)) * lambda / d ); 

% Equa��es que descrevem o primeiro elips�ide de Fresnel
yup = rd + r1;
ydwn = rd - r1;

% Desenho do primeiro elips�ide de Fresnel
plot(dx, ydwn, 'k');
plot(dx, yup, 'k');

for i = 1 : length(ydwn)
    if ydwn(i) <= y(i)
        fprintf('O Elips�ide de Fresnel encontra-se obstru�do\n\n');
        break;
    end
end

%% 3. C�lculo do ganho das antenas

fprintf('C�LCULO DO GANHO DAS ANTENAS\n\n');

% Ganho da antena emissora (em dBi)
Ge = 20 * log10( (pi * le) / lambda ) + 10 * log10(ne);
fprintf('Ganho da antena emissora: %.5g dBi\n', Ge);

% Ganho da antena recetora (em dBi)
Gr = 20 * log10( (pi * lr) / lambda ) + 10 * log10(nr);
fprintf('Ganho da antena recetora: %.5g dBi\n\n', Gr);

%% 4. C�lculo da atenua��o nos guias de ondas

fprintf('C�LCULO DA ATENUA��O NOS GUIAS DE ONDAS\n\n');

% Atenua��o no guia de onda da antena emissora (em dB)
Ae = (cge/1e3) * Aekm;
fprintf('Atenua��o no guia de onda da antena emissora: %.5g dB\n', -Ae);

% Atenua��o no guia de onda da antena recetora (em dB)
Ar = (cgr/1e3) * Arkm;
fprintf('Atenua��o no guia de onda da antena recetora: %.5g dB\n\n', -Ar);

%% 5. C�lculo da atenua��o em espa�o livre

fprintf('C�LCULO DA ATENUA��O EM ESPA�O LIVRE\n\n');

% Atenua��o em espa�o livre (em dB)
A0 = 32.4 + 20*log10(dKm) + 20*log10(fMhz);
fprintf('Atenua��o em espa�o livre: %.5g dB\n\n', -A0);

%% 6. C�lculo das horizontais � superf�cie da Terra e �ngulos de fogo

fprintf('C�LCULO DOS DECLIVES DAS RETAS TANGENTES � SUPERFICIE DA TERRA\n\n');

if model == 1 % Modelo: Terra plana
    % C�lculo dos declives das retas tangentes � superf�cie da Terra
    dydxe = 0;
    dydxr = 0;
elseif model == 2 && system == 1 % Modelo: Terra esf�rica, Sistema: Europeu
    % C�lculo dos declives das retas tangentes � superf�cie da Terra
    dydxe = d/(2 * r0);
    dydxr = - (d - d/2) / r0;
elseif model == 2 && system == 2 % Modelo: Terra esf�rica, Sistema: Americano
    % C�lculo dos declives das retas tangentes � superf�cie da Terra
    dydxe = 0;
    dydxr = -d / r0;
end
   
% Declive da reta tangente � superf�cie da Terra na antena emissora
fprintf('Declive da reta tangente � superf�cie da Terra na antena emissora: %.5g\n', dydxe);

% Declive da reta tangente � superf�cie da Terra na antena recetora
fprintf('Declive da reta tangente � superf�cie da Terra na antena recetora: %.5g\n\n', dydxe);

fprintf('C�LCULO DOS �NGULOS DE FOGO\n\n');

if model == 1 % Modelo: Terra plana
    if(he < hr)
        % �ngulo de fogo na antena emissora
        afe = abs(atand(md));
        % �ngulo de fogo na antena recetora
        afr = -afe;
    elseif(he > hr)
        % �ngulo de fogo na antena emissora
        afe = atand(md);
        % �ngulo de fogo na antena recetora
        afr = -afe;
    elseif(he == hr)
        % �ngulo de fogo na antena emissora
        afe = 0;
        % �ngulo de fogo na antena recetora
        afr = 0;
    end
elseif model == 2 % Modelo: Terra esf�rica
    if(het < hrt)
        % �ngulo de fogo na antena emissora
        afe = - ( abs(atand(dydxe)) - abs(atand(md)) );
        % �ngulo de fogo na antena recetora
        afr = - ( abs(atand(dydxr)) + abs(atand(md)) );
    elseif(het > hrt)
        % �ngulo de fogo na antena emissora
        afe = - ( abs(atand(dydxe)) + abs(atand(md)) );
        % �ngulo de fogo na antena recetora
        afr = abs(abs(atand(dydxr)) - abs(atand(md)));
    elseif(het == hrt)
        % �ngulo de fogo na antena emissora
        afe = - atand(dydxe);
        % �ngulo de fogo na antena recetora
        afr = atand(dydxr);
    end
end

fprintf('�ngulo de fogo na antena emissora: %.5g�\n', afe);
fprintf('�ngulo de fogo na antena recetora: %.5g�\n\n', afr);

%% 7. C�lculo da localiza��o do ponto de reflex�o e do fator de diverg�ncia
  
fprintf('C�LCULO DA LOCALIZA��O DO PONTO DE REFLEX�O E DO FATOR DE DIVERG�NCIA\n\n');

if model == 1 % Modelo: Terra plana
    % C�lculo da localiza��o do ponto de reflex�o
    xe = (he/(he + hr))*d;
    % C�lculo do �ngulo de incid�ncia
    psi = atand((he+hr)/d);
    % C�lculo do fator de diverg�ncia
    Dv = 1;
elseif model == 2 % Modelo: Terra esf�rica
    % C�lculo da localiza��o do ponto de reflex�o
    xe = roots([1 -1.5*d (0.5*(d^2) - r0*(he+hr)) r0*het*d]);
    for i = 1:length(xe)
        if xe(i) >= 0 && xe(i) <= d
            xe = xe(i);
            break;
        end
    end
    
    % C�lculo do �ngulo de incid�ncia
    psi = atand((he - (xe^2/(2*r0)))/xe);
    
    % C�lculo do fator de diverg�ncia
    Dv = 1/(sqrt(1+((2*xe*(d-xe))/(r0*d*sind(psi)))));    
end

% Ponto de reflex�o (ou especular)
fprintf('Ponto de reflex�o (ou especular): %.5gm\n', xe);

% Fator de Diverg�ncia
fprintf('Fator de Diverg�ncia: %.5g\n\n', Dv);

%% 8. C�lculo da diferen�a de fase entre o campo associado ao raio direto e ao associado ao raio refletido no terreno

fprintf('C�LCULO DA DIFEREN�A DE FASE ENTRE O CAMPO ASSOCIADO AO RAIO DIRETO E AO ASSOCIADO AO RAIO REFLETIDO NO TERRENO\n\n');

% �ndice de reflex�o do solo em rela��o ao ar 
n = sqrt((terrainTable(terrain, 1) * e0 - (terrainTable(terrain, 2) / 2*pi*f) * 1j)/e0);

% Coeficiente de Fresnel
if polarizacao == 1 % Polariza��o: Horizontal
    % C�lculo do do coeficiente de fresnel
    coefFresnel = (sin(psi) - sqrt(n^2 - (cos(psi)^2)) / (sin(psi) + sqrt(n^2 - (cos(psi)^2))));
elseif polarizacao == 2 % Polariza��o: Vertical
    % C�lculo do coeficiente de fresnel
    coefFresnel = (n^2 * sin(psi) - sqrt(n^2 - (cos(psi)^2)) / (n^2 * sin(psi) + sqrt(n^2 - (cos(psi)^2))));
end

% C�lculo do argumento do coeficiente de fresnel (em graus)
argCoefFresnel = angle(coefFresnel);

% C�lculo da constante de propaga��o de fase
k0 = (-2*pi)/ lambda;

if model == 1 % Modelo: Terra plana
    % C�lculo da diferen�a de percurso entre o raio direto e refletido
    deltaR = (2 * he * hr) / d;
elseif model == 2 % Modelo: Terra esf�rica
    
    % C�lculo da localiza��o do ponto de reflex�o
    xe = roots( [1 -1.5*d (0.5*(d^2) - re*(he+hr)) re*he*d] );
    for j = 1 : length(xe)
        if xe(j) >= 0 && xe(j) <= d
            xe = xe(j);
            break;
        end
    end

    % Altura equivalente do emissor 
    hetxe = he - ( (xe^2) / (2 * re) );
    % Altura equivalente do recetor
    hrtxe = hr - ( ((d - xe)^2) / (2 * re) );
    
    % C�lculo da diferen�a de percurso entre o raio direto e refletido
    deltaR = (2 * hetxe * hrtxe) / d;
end 

% Diferen�a de fase entre o campo associado ao raio direto e ao associado ao raio refletido no terreno (em graus)
phi = -k0 * deltaR + argCoefFresnel;
fprintf('Diferen�a de fase entre o campo associado ao raio direto e ao associado ao raio refletido no terreno: %.5g\n\n', phi);

%% 9. Diferen�a de fase entre terra plana e terra esf�rica

fprintf('C�LCULO DA DIFEREN�A DE FASE ENTRE TERRA PLANA E TERRA ESF�RICA\n\n');

if model == 2
    dphi = (2*pi/lambda)*(((2*he*hr)/d)-((2*hetxe*hrtxe)/d));
    fprintf('Diferen�a de fase entre terra plana e terra esf�rica: %.5g\n\n', dphi);
end
    
%% 10. C�lculo da pot�ncia total recebida

fprintf('C�LCULO DA POT�NCIA TOTAL RECEBIDA\n\n');

% C�lculo da pot�ncia associada ao raio direto
Pd = Pe + Ge - Ae + Gr - Ar - A0;

if model == 1 % Modelo: Terra plana
    % C�lculo da pot�ncia total recebida
    Pr = Pd + 20 * log10( abs( 2 * sind( (2*pi*he.*hr) / (lambda*d) ) ) );
elseif model == 2 % Modelo: Terra esf�rica
    % C�lculo da pot�ncia total recebida
    Pr = Pd + 10 * log10( abs( 1 + Dv^2 - 2 * Dv * cosd( (4*pi*hetxe.*hrtxe) / (lambda*d) ) ) );
end
    
fprintf('Pot�ncial total recebida: %.5gdBw\n\n', Pr);

%% 11. C�lculo da atenua��o provocada pelos gases atmosf�ricos

fprintf('C�LCULO DA ATENUA��O PROVOCADA PELOS GASES ATMOSF�RICOS\n\n');

% C�lculo da atenua��o provocada pelo vapor de �gua

% C�lculo da press�o relativa � atmosfera padr�o (em atmosferas)
% 1013: Press�o m�dia da Atmosfera pr�xima da superf�cie terrestre (em hectopascais)
rp = P/1013;

% C�lculo da temperatura relativa � atmosfera padr�o
% 288: Temperatura m�dia da Atmosfera (em kelvins) 
rt = 288/(T+273);

% C�lculo da press�o parcial do vapor de �gua saturado (em hectopascais)
ppvas = 6.1121*exp(17.502*T/(240.97+T));

% C�lculo da press�o parcial do vapor de �gua (em hectopascais)
ppva = H*ppvas;

% C�lculo da densidade do vapor de �gua
rho = 216.7 * (ppva/(T+273.3));

% C�clulo do coeficiente de atenua��o provocado pelo vapor de �gua (em dB/km)
gammaw = ((3.27e-2)*rt + (1.67e-3) * rho*(rt^7/rp) + (7.7e-4) * (fGhz^0.5) + (3.79/(((fGhz-22.235)^2) + 9.81*(rp^2) * rt))+((11.73 * rt)/(((fGhz-183.31)^2) + 11.85*(rp^2)*rt))+((4.01*rt)/(((fGhz-325.153)^2) + 10.44*(rp^2) * rt)))*(fGhz^2)*rho*rp*rt*1e-4; 

% C�lculo da atenua��o provocada pelo vapor de �gua (em dB)
Atw = gammaw * dKm;

fprintf('Atenua��o provocada pelo vapor de �gua: %.5gdB\n', Atw);

% C�lculo da atenua��o provocada pelo oxig�nio

% C�clulo do coeficiente de atenua��o provocado pelo oxig�nio (em dB/km)

% Para freq. menores ou iguais a 57GHz
if fGhz <= 57
    gammao = (((7.27*rt)/((fGhz^2)+0.351*(rp^2)*(rt^2))) + (7.5/(((fGhz - 57)^2)+2.44*(rp^2)*(rt^5))))*(fGhz^2)*(rp^2)*(rt^2)*1e-3;
% Para freq. superiores a 57GHz e inferiores a 63GHz
elseif fGhz > 57 && fGhz < 63
    gammao = (((fGhz-60)*(fGhz-63))/18)*((((7.27*rt)/((57^2)+0.351*(rp^2)*(rt^2))) + (7.5/(((57 - 57)^2)+2.44*(rp.^2).*(rt^5))))*(57^2)*(rp^2)*(rt^2)*1e-3)-1.66*(rp^2)*(rt^8.5)*(fGhz-57)*(fGhz-63) + (((fGhz-57)*(fGhz-60))/18)*((2*1e-4*(rt^1.5)*(1-1.2e-5*(63^1.5))+(4/(((63- 63)^2)+1.5*(rp^2)*(rt^5)))+((0.28*(rt^2))/(((63 - 118.75)^2)+2.84*(rp^2)*(rt^2))))*(63^2)*(rp^2)*(rt^2)*1e-3);
% Para freq. superiores ou iguais a 350 GHz
elseif fGhz >= 350 
    gammao = ((2*1e-4*(rt^1.5)*(1-1.2e-5*(fGhz^1.5))+(4/(((fGhz- 63)^2)+1.5*(rp^2)*(rt^5)))+((0.28*(rt^2))/(((fGhz - 118.75)^2)+2.84*(rp^2)*(rt^2))))*(fGhz^2)*(rp^2)*(rt^2)*1e-3);
end

% C�lculo da atenua��o provocada pelo oxig�nio (em dB)
Ato = gammao * d;
fprintf('Atenua��o provocada pelo oxig�nio: %.5gdB\n', Ato);

% C�lculo da atenua��o suplementar devida � presen�a da atmosfera (em dB)
Atg = Ato + Atw;
fprintf('Atenua��o suplementar devida � presen�a da atmosfera: %.5gdB\n\n', Atg);

%% 12. C�lculo da atenua��o provocada pela precipita��o

fprintf('C�LCULO DA ATENUA��O PROVOCADA PELA PRECIPITA��O\n\n');

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
gammar = k * R.^a;

% C�lculo do comprimento eficaz do percurso 
def = d / (1 + (d / 35 *( exp(-0.015 * R) ) ) );

% C�lculo da atenua��o provocada pela precipita��o n�o excedida em 0.01% do tempo (relativamente � m�dia anual)
Ar1 = gammar * def;

% C�lculo da atenua��o provocada pela precipita��o n�o excedida numa determinada percentagem do tempo (relativamente � m�dia anual)
Arp = Ar1 * 0.12 * p^(-(0.546 + 0.043 * log10(p))); 
fprintf('Atenua��o provocada pela precipita��o n�o excedida numa determinada percentagem do tempo (relativamente � m�dia anual): %.5gdB\n', Arp);

% C�lculo da atenua��o provocada pela precipita��o n�o excedida numa determinada percentagem do tempo (relativamente ao pior m�s)
pm = 0.3 * p^1.15;
Arpm = Ar1 * 0.12 * pm^(-(0.546 + 0.043 * log10(pm))); 
fprintf('Atenua��o provocada pela precipita��o n�o excedida numa determinada percentagem do tempo (relativamente ao pior m�s): %.5gdB\n\n', Arpm);

%% 13. Modela��o dos efeitos refrativos na atmosfera

fprintf('MODELA��O DOS EFEITOS REFRATIVOS NA ATMOSFERA\n\n');

% C�lculo da refratividade
N = (77.6/(T+273))*(P + 4810*(ppva/(T+273)));

% C�lculo do �ndice de refrac��o
n = 1+N*1e-6;

% C�lculo do �ndice de refrac��o modificado M
h = 0 : 100 : 10e3;
M = N + 1e6.*(h./r0);

figure;
plot(M,h);
grid on;
title('M(h)');
xlabel('M');
ylabel('h[km]');

%% 14. Difrac��o em Terra Esf�rica

% Efeito do Solo e Polariza��o
beta = (1 + 1.6*(k^2) + 0.75*(k^4))/(1 + 4.5*(k^2) + 1.35*(k^4));

% Fun��o x
x = beta * d *((pi/(lambda*(r0^2)))^(1/3));
% Fun��o F(x)
Fx = 11 + 10 * log10(x) - 17.6 * x;

% Fun��o y1
y1 = 2 * beta * het *(((pi^2)/((lambda^2)*r0))^(1/3));

% Fun��o y2
y2 = 2 * beta * hrt *(((pi^2)/((lambda^2)*r0))^(1/3));

% C�lculo 
KH = (((2*pi*r0)/lambda)^(-1/3))*(((terrainTable(terrain, 1)-1)^2)+(60*lambda*terrainTable(terrain, 2))^2)^(-1/4);
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

%% 16. Perdas por Difra��o devido ao terreno usando o m�todo de Deygout

if length(perfilTerreno) ~= dKm
    error('N�mero de pontos definidos no ficheiro excel n�o pode ser diferente da dist�ncia (em quil�metros)')
end

distanceArray = 1 : dKm;
plot(distanceArray, perfilTerreno);

AtenuacaoObstaculos = Deygout(d, het, hrt, lambda, dKm, perfilTerreno, re, 0, 0) + A0;
