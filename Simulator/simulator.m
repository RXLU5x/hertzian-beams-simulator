%% Inicialização da aplicação

clearvars
close all
clc

%% Definição das constantes

% Raio da terra (em metros)
r0 = 6370e3;

% Velocidade da luz no vazio (em metros por segundo)
c = 3e8;

% Permitividade elétrica do vácuo (em farad por metro)
e0 = 8.854e-12;

% Indíce de refração a uma altura de 0km
n0 = 1 + 315e-6;

%% Definição das tabelas

% Tabela 1: Tipo de terreno

% 1ª linha: Mar
% 2ª linha: Muito Húmido
% 3ª linha: Normal
% 4ª linha: Seco

% 1ª coluna: Permitividade relativa para cada tipo de terreno
% 2ª coluna: Condutividade elétrica para cada tipo de terreno (em siemens por metro)

terrainTable = [
    81 5 ; 
    25 0.02 ;
    15 0.005 ;
    4 0.001 
];

% Tabela 2: k e alpha

% 1ª linha: 1 gigahertz
% 5ª linha: 5 gigahertz
% 10ª linha: 10 gigahertz
% 20ª linha: 20 gigahertz
% 60ª linha: 60 gigahertz

% 1ª coluna: kH
% 2ª coluna: kV
% 3ª coluna: alphaH
% 4ª coluna: alphaV

kAlphaTable(1, 1:4) = [0.0000387 0.0000352 0.912 0.880];
kAlphaTable(4, 1:4) = [0.000650 0.000591 1.121 1.075];
kAlphaTable(6, 1:4) = [0.00175 0.00155 1.308 1.265];
kAlphaTable(10, 1:4) = [0.0101 0.00887 1.276 1.264];
kAlphaTable(20, 1:4) = [0.0751 0.0691 1.099 1.065];
kAlphaTable(60, 1:4) = [0.707 0.642 0.826 0.824];

% Cálculo por interpolação do valor de k e alpha para a frequência de 5Ghz

alphaH5Ghz = kAlphaTable(4, 3) + ( ( ( kAlphaTable(6, 3) - kAlphaTable(4, 3) ) / ( log10(6) - log10(4) ) ) * ( log10(5) - log10(4) ) );
alphaV5Ghz = kAlphaTable(4, 4) + ( ( ( kAlphaTable(6, 4) - kAlphaTable(4, 4) ) / ( log10(6) - log10(4) ) ) * ( log10(5) - log10(4) ) );

kH5Ghz = 10^( log10(kAlphaTable(4,1)) + ( ( ( log10(kAlphaTable(6,1)) - log10(kAlphaTable(4,1)) ) / ( log10(6) - log10(4) ) ) * ( log10(5) - log10(4) ) ) );
kV5Ghz = 10^( log10(kAlphaTable(4,2)) + ( ( ( log10(kAlphaTable(6,2)) - log10(kAlphaTable(4,2)) ) / ( log10(6) - log10(4) ) ) * ( log10(5) - log10(4) ) ) );

kAlphaTable(5, 1:4) = [kH5Ghz kV5Ghz alphaH5Ghz alphaV5Ghz];

% Tabela 3: 
% Ficheiro excel

perfilTerreno = (xlsread('terrainProfile.xlsx'))';

%% Configuração do modelo (e sistema)

fprintf('CONFIGURAÇÃO DO MODELO\n\n');

% Modelo
fprintf('Selecione o modelo que pretende simular:\n\n1: Terra plana\n2: Terra esférica\n\n');
model = abs(input('Modelo: '));
disp(' ');

switch model
    case 1
        fprintf('Selecionou como modelo: Terra plana\n\n');
    case 2
        fprintf('Selecionou como modelo: Terra esférica\n\n');
        
        fprintf('CONFIGURAÇÃO DO SISTEMA\n\n');
        
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
                error('Sistema selecionado não existe');
        end
    otherwise
        error('Modelo selecionado não existe');
end

%% Configuração das antenas

fprintf('CONFIGURAÇÃO DAS ANTENAS\n\n');

% Polarização
fprintf('Selecione a polarização:\n\n1: Horizontal\n2: Vertical\n\n');
polarizacao = abs(input('Polarização: '));
disp(' ');

switch polarizacao
    case 1
        fprintf('Selecionou como polarização: Horizontal\n\n');
    case 2
        fprintf('Selecionou como polarização: Vertical\n\n');
    otherwise
        error('Polarização selecionada não existe');
end

% Distância entre as antenas (em metros)
d = abs(input('Distância entre as antenas (em metros): '));
disp(' ');

%% Configuração da antena emissora

fprintf('CONFIGURAÇÃO DA ANTENA EMISSORA\n\n');

% Frequência (em hertz)
% 1: 1e9 Hz
% 2: 5e9 Hz
% 3: 10e9 Hz
% 4: 20e9 Hz
% 5: 60e9 Hz
fprintf('Selecione a frequência:\n\n1: 1Ghz\n2: 5Ghz\n3: 10Ghz\n4: 20Ghz\n5: 60Ghz\n\n');
f = abs(input('Frequência: '));
disp(' ');

switch f
    case 1
        fprintf('Selecionou como frequência: 1Ghz\n\n');
        f = 1e9;
    case 2
        fprintf('Selecionou como frequência: 5Ghz\n\n');
        f = 5e9;
    case 3
        fprintf('Selecionou como frequência: 10Ghz\n\n');
        f = 10e9;
    case 4
        fprintf('Selecionou como frequência: 20Ghz\n\n');
        f = 20e9;
    case 5
        fprintf('Selecionou como frequência: 60Ghz\n\n');
        f = 60e9;
    otherwise
        error('Frequência inválida');
end

% Altura (em metros)
he = abs(input('Altura (em metros): '));

% Diâmetro do prato (em metros)
le = abs(input('Diâmetro do prato (em metros): ')); 

% Potência (em decibéis)
Pe = 10*log10(abs(input('Potência (em watts): ')));

% Rendimento (em percentagem)
ne = abs(input('Rendimento (em percentagem): '));

% Comprimento do guia de onda (em metros)
cge = abs(input('Comprimento do guia de onda (em metros): '));

% Atenuação no guia de onda (em decibéis por quilómetro)
Aekm = abs(input('Atenuação no guia de onda (em decibéis por quilómetro): '));

disp(' ');

%% Configuração da antena recetora

fprintf('CONFIGURAÇÃO DA ANTENA RECETORA\n\n');

% Altura (em metros)
hr = abs(input('Altura (em metros): '));

% Diâmetro do prato (em metros)
lr = abs(input('Diâmetro do prato (em metros): ')); 

% Rendimento (em percentagem)
nr = abs(input('Rendimento (em percentagem): '));

% Comprimento do guia de onda (em metros)
cgr = abs(input('Comprimento do guia de onda (em metros): '));

% Atenuação no guia de onda (em decibéis por quilómetro)
Arkm = abs(input('Atenuação no guia de onda (em decibéis por quilómetro): '));

disp(' ');

%% Configuração do meio ambiente

fprintf('CONFIGURAÇÃO DO MEIO AMBIENTE\n\n');

% Terreno:
% 1: Mar
% 2: Muito Húmido
% 3: Normal
% 4: Seco
fprintf('Selecione o tipo de terreno:\n\n1: Mar\n2: Muito Húmido\n3: Normal\n4: Seco\n\n');
terrain = abs(input('Terreno: '));
disp(' ');

switch terrain
    case 1
        fprintf('Selecionou como tipo de terreno: Mar\n\n');
    case 2
        fprintf('Selecionou como tipo de terreno: Muito Húmido\n\n');
    case 3
        fprintf('Selecionou como tipo de terreno: Normal\n\n');
    case 4
        fprintf('Selecionou como tipo de terreno: Seco\n\n');
    otherwise
        error('Terreno selecionado não existe');
end

% Região:
% 1: H
% 2: K
fprintf('Selecione a região:\n\n1: H\n2: K\n\n');
region = abs(input('Região: '));
disp(' ');

switch region
    case 1
        fprintf('Selecionou como região: H\n\n');
    case 2
        fprintf('Selecionou como região: K\n\n');
    otherwise
        error('Região selecionada não existe\n\n');
end

% Pressão atmosférica (em hectopascais)
P = abs(input('Pressão atmosférica (em hectopascais): '));

% Temperatura (em graus centígrados)
T = abs(input('Temperatura (em graus centígrados): '));

% Humidade relativa (em percentagem)
H = abs(input('Humidade relativa (em percentagem): '));

% Variação do índice de refração a uma altura de 0km (dn/dh) (em km^-1)
% Valor para atmosfera padrão: -43e-6
dndh0 = abs(input('Variação do índice de refração (em km^-1): ')) * 1e-3;

% Percentagem de tempo, relativamente à média anual, em que o valor da intensidade da precipitação é excedido (em percentagem)
p = abs(input('Percentagem de tempo, relativamente à média anual, em que o valor da intensidade da precipitação é excedido (em percentagem): '));

disp(' ');

figure;
title('Simulador');
xlabel('Distância (m)');
ylabel('Altura (m)');
grid on;
hold on;

%% Cálculos auxiliares

% Distância (em quilómetros)
dKm = d/1e3;

% Frequência (em gigahertz)
fGhz = f/1e9;

% Frequência (em megahertz)
fMhz = f/1e6;

% Comprimento de onda (em metros)
lambda = c/f;

% Variação do eixo das abcissas (em metros)
dx = 0 : d/1e6 : d;

%% Desenho do referencial

if model == 1 % Modelo: Terra plana
    % Cálculo do decaimento da terra para cada valor de dx
    y = zeros(1, length(dx));
    
    % Cálculo do declive da reta tangente à superficie terrestre no emissor
    dydxe = 0;

    % Cálculo do declive da reta tangente à superficie terrestre no recetor
    dydxr = 0;
elseif model == 2 && system == 1 % Modelo: Terra esférica, Sistema: Europeu
    % Cálculo do raio equivalente a partir do fator K e da variação do indice de refração
    Kre = 1/(1+((r0/n0)*dndh0));
    r0 = Kre * r0;
    
    % Cálculo do decaimento da terra para cada valor de dx
    y = -( (dx - d/2).^2 )/(2*r0);
    
    % Altura equivalente do emissor 
    het = he + y(0);
    % Altura equivalente do recetor
    hrt = hr + y(length(dx));
    
    % Cálculo do declive da reta tangente à superficie terrestre no emissor
    dydxe = d/2 / r0;

    % Cálculo do declive da reta tangente à superficie terrestre no recetor
    dydxr = - (d - d/2) / r0;
    
elseif model == 2 && system == 2 % Modelo: Terra esférica, Sistema: Americano
    % Cálculo do raio equivalente a partir do fator K e da variação do indice de refração
    Kre = 1/(1+((r0/n0)*dndh0));
    r0 = Kre * r0;
    
    % Cálculo do decaimento da terra para cada valor de dx
    y = -(dx.^2)/(2*r0);

    % Altura equivalente do emissor 
    het = he + y(0);
    % Altura equivalente do recetor
    hrt = hr + y(length(dx));
    
    % Cálculo do declive da reta tangente à superficie terrestre no emissor
    dydxe = 0;

    % Cálculo do declive da reta tangente à superficie terrestre no recetor
    dydxr = -d / r0;
end

% Desenho da superfície da Terra
plot(dx, y, 'k');

%% Desenho da horizontal à superfície da Terra no emissor e no recetor

if model == 2
    % Cálculo da ordenada na origem
    bte = het;

    % Cálculo da equação da reta paralela à reta tangente em relação à superficie terrestre no emissor
    yte = dydxe .* dx + bte;

    % Desenho da horizontal à superfície da Terra no emissor
    plot(dx, yte, [':' 'k']);

    % Cálculo da ordenada na origem
    btr = hrt - dydxr .* d;

    % Cálculo da equação da reta paralela à reta tangente em relação à superficie terrestre no recetor
    ytr = dydxr .* dx + btr;

    % Desenho da horizontal à superfície da Terra no recetor
    plot(dx, ytr, [':' 'k']);
end

%% 1. Desenho do Raio Direto

fprintf('DESENHO DO RAIO DIRETO\n\n');
    
if model == 1 % Modelo: Terra plana
    % Cálculo do declive da equação da reta que descreve o raio direto
    md = (hr - he) / d;

    % Cálculo da ordenada na origem
    bd = he;
else % Modelo: Terra esférica
    % Cálculo do declive da equação da reta que descreve o raio direto
    md = (hrt - het) / d;

    % Cálculo da ordenada na origem
    bd = het;
end

% Cálculo da equação da reta que descreve o raio direto (yd = md*x + bd)
rd = md * dx + bd;

% Desenho do raio direto
plot(dx, rd, ['--' 'k']);
    
%% 2. Desenho do primeiro elipsóide de fresnel e indicação se este se econtra obstruído

fprintf('DESENHO DO PRIMEIRO ELIPSÓIDE DE FRESNEL E INDICAÇÃO SE ESTE SE ECONTRA OBSTRUÍDO\n\n');

% Raio de cada secção do elipsóide para cada valor de dx (em metros)
r1 = sqrt( (dx .* (d-dx)) * lambda / d ); 

% Equações que descrevem o primeiro elipsóide de Fresnel
yup = rd + r1;
ydwn = rd - r1;

% Desenho do primeiro elipsóide de Fresnel
plot(dx, ydwn, 'k');
plot(dx, yup, 'k');

for i = 1 : length(ydwn)
    if ydwn(i) <= y(i)
        fprintf('O Elipsóide de Fresnel encontra-se obstruído\n\n');
        break;
    end
end

%% 3. Cálculo do ganho das antenas

fprintf('CÁLCULO DO GANHO DAS ANTENAS\n\n');

% Ganho da antena emissora (em dBi)
Ge = 20 * log10( (pi * le) / lambda ) + 10 * log10(ne);
fprintf('Ganho da antena emissora: %.5g dBi\n', Ge);

% Ganho da antena recetora (em dBi)
Gr = 20 * log10( (pi * lr) / lambda ) + 10 * log10(nr);
fprintf('Ganho da antena recetora: %.5g dBi\n\n', Gr);

%% 4. Cálculo da atenuação nos guias de ondas

fprintf('CÁLCULO DA ATENUAÇÃO NOS GUIAS DE ONDAS\n\n');

% Atenuação no guia de onda da antena emissora (em dB)
Ae = (cge/1e3) * Aekm;
fprintf('Atenuação no guia de onda da antena emissora: %.5g dB\n', -Ae);

% Atenuação no guia de onda da antena recetora (em dB)
Ar = (cgr/1e3) * Arkm;
fprintf('Atenuação no guia de onda da antena recetora: %.5g dB\n\n', -Ar);

%% 5. Cálculo da atenuação em espaço livre

fprintf('CÁLCULO DA ATENUAÇÃO EM ESPAÇO LIVRE\n\n');

% Atenuação em espaço livre (em dB)
A0 = 32.4 + 20*log10(dKm) + 20*log10(fMhz);
fprintf('Atenuação em espaço livre: %.5g dB\n\n', -A0);

%% 6. Cálculo das horizontais à superfície da Terra e ângulos de fogo

fprintf('CÁLCULO DOS DECLIVES DAS RETAS TANGENTES À SUPERFICIE DA TERRA\n\n');

if model == 1 % Modelo: Terra plana
    % Cálculo dos declives das retas tangentes à superfície da Terra
    dydxe = 0;
    dydxr = 0;
elseif model == 2 && system == 1 % Modelo: Terra esférica, Sistema: Europeu
    % Cálculo dos declives das retas tangentes à superfície da Terra
    dydxe = d/(2 * r0);
    dydxr = - (d - d/2) / r0;
elseif model == 2 && system == 2 % Modelo: Terra esférica, Sistema: Americano
    % Cálculo dos declives das retas tangentes à superfície da Terra
    dydxe = 0;
    dydxr = -d / r0;
end
   
% Declive da reta tangente à superfície da Terra na antena emissora
fprintf('Declive da reta tangente à superfície da Terra na antena emissora: %.5g\n', dydxe);

% Declive da reta tangente à superfície da Terra na antena recetora
fprintf('Declive da reta tangente à superfície da Terra na antena recetora: %.5g\n\n', dydxe);

fprintf('CÁLCULO DOS ÂNGULOS DE FOGO\n\n');

if model == 1 % Modelo: Terra plana
    if(he < hr)
        % Ângulo de fogo na antena emissora
        afe = abs(atand(md));
        % Ângulo de fogo na antena recetora
        afr = -afe;
    elseif(he > hr)
        % Ângulo de fogo na antena emissora
        afe = atand(md);
        % Ângulo de fogo na antena recetora
        afr = -afe;
    elseif(he == hr)
        % Ângulo de fogo na antena emissora
        afe = 0;
        % Ângulo de fogo na antena recetora
        afr = 0;
    end
elseif model == 2 % Modelo: Terra esférica
    if(het < hrt)
        % Ângulo de fogo na antena emissora
        afe = - ( abs(atand(dydxe)) - abs(atand(md)) );
        % Ângulo de fogo na antena recetora
        afr = - ( abs(atand(dydxr)) + abs(atand(md)) );
    elseif(het > hrt)
        % Ângulo de fogo na antena emissora
        afe = - ( abs(atand(dydxe)) + abs(atand(md)) );
        % Ângulo de fogo na antena recetora
        afr = abs(abs(atand(dydxr)) - abs(atand(md)));
    elseif(het == hrt)
        % Ângulo de fogo na antena emissora
        afe = - atand(dydxe);
        % Ângulo de fogo na antena recetora
        afr = atand(dydxr);
    end
end

fprintf('Ângulo de fogo na antena emissora: %.5gº\n', afe);
fprintf('Ângulo de fogo na antena recetora: %.5gº\n\n', afr);

%% 7. Cálculo da localização do ponto de reflexão e do fator de divergência
  
fprintf('CÁLCULO DA LOCALIZAÇÃO DO PONTO DE REFLEXÃO E DO FATOR DE DIVERGÊNCIA\n\n');

if model == 1 % Modelo: Terra plana
    % Cálculo da localização do ponto de reflexão
    xe = (he/(he + hr))*d;
    % Cálculo do ângulo de incidência
    psi = atand((he+hr)/d);
    % Cálculo do fator de divergência
    Dv = 1;
elseif model == 2 % Modelo: Terra esférica
    % Cálculo da localização do ponto de reflexão
    xe = roots([1 -1.5*d (0.5*(d^2) - r0*(he+hr)) r0*het*d]);
    for i = 1:length(xe)
        if xe(i) >= 0 && xe(i) <= d
            xe = xe(i);
            break;
        end
    end
    
    % Cálculo do ângulo de incidência
    psi = atand((he - (xe^2/(2*r0)))/xe);
    
    % Cálculo do fator de divergência
    Dv = 1/(sqrt(1+((2*xe*(d-xe))/(r0*d*sind(psi)))));    
end

% Ponto de reflexão (ou especular)
fprintf('Ponto de reflexão (ou especular): %.5gm\n', xe);

% Fator de Divergência
fprintf('Fator de Divergência: %.5g\n\n', Dv);

%% 8. Cálculo da diferença de fase entre o campo associado ao raio direto e ao associado ao raio refletido no terreno

fprintf('CÁLCULO DA DIFERENÇA DE FASE ENTRE O CAMPO ASSOCIADO AO RAIO DIRETO E AO ASSOCIADO AO RAIO REFLETIDO NO TERRENO\n\n');

% Índice de reflexão do solo em relação ao ar 
n = sqrt((terrainTable(terrain, 1) * e0 - (terrainTable(terrain, 2) / 2*pi*f) * 1j)/e0);

% Coeficiente de Fresnel
if polarizacao == 1 % Polarização: Horizontal
    % Cálculo do do coeficiente de fresnel
    coefFresnel = (sin(psi) - sqrt(n^2 - (cos(psi)^2)) / (sin(psi) + sqrt(n^2 - (cos(psi)^2))));
elseif polarizacao == 2 % Polarização: Vertical
    % Cálculo do coeficiente de fresnel
    coefFresnel = (n^2 * sin(psi) - sqrt(n^2 - (cos(psi)^2)) / (n^2 * sin(psi) + sqrt(n^2 - (cos(psi)^2))));
end

% Cálculo do argumento do coeficiente de fresnel (em graus)
argCoefFresnel = angle(coefFresnel);

% Cálculo da constante de propagação de fase
k0 = (-2*pi)/ lambda;

if model == 1 % Modelo: Terra plana
    % Cálculo da diferença de percurso entre o raio direto e refletido
    deltaR = (2 * he * hr) / d;
elseif model == 2 % Modelo: Terra esférica
    
    % Cálculo da localização do ponto de reflexão
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
    
    % Cálculo da diferença de percurso entre o raio direto e refletido
    deltaR = (2 * hetxe * hrtxe) / d;
end 

% Diferença de fase entre o campo associado ao raio direto e ao associado ao raio refletido no terreno (em graus)
phi = -k0 * deltaR + argCoefFresnel;
fprintf('Diferença de fase entre o campo associado ao raio direto e ao associado ao raio refletido no terreno: %.5g\n\n', phi);

%% 9. Diferença de fase entre terra plana e terra esférica

fprintf('CÁLCULO DA DIFERENÇA DE FASE ENTRE TERRA PLANA E TERRA ESFÉRICA\n\n');

if model == 2
    dphi = (2*pi/lambda)*(((2*he*hr)/d)-((2*hetxe*hrtxe)/d));
    fprintf('Diferença de fase entre terra plana e terra esférica: %.5g\n\n', dphi);
end
    
%% 10. Cálculo da potência total recebida

fprintf('CÁLCULO DA POTÊNCIA TOTAL RECEBIDA\n\n');

% Cálculo da potência associada ao raio direto
Pd = Pe + Ge - Ae + Gr - Ar - A0;

if model == 1 % Modelo: Terra plana
    % Cálculo da potência total recebida
    Pr = Pd + 20 * log10( abs( 2 * sind( (2*pi*he.*hr) / (lambda*d) ) ) );
elseif model == 2 % Modelo: Terra esférica
    % Cálculo da potência total recebida
    Pr = Pd + 10 * log10( abs( 1 + Dv^2 - 2 * Dv * cosd( (4*pi*hetxe.*hrtxe) / (lambda*d) ) ) );
end
    
fprintf('Potêncial total recebida: %.5gdBw\n\n', Pr);

%% 11. Cálculo da atenuação provocada pelos gases atmosféricos

fprintf('CÁLCULO DA ATENUAÇÃO PROVOCADA PELOS GASES ATMOSFÉRICOS\n\n');

% Cálculo da atenuação provocada pelo vapor de água

% Cálculo da pressão relativa à atmosfera padrão (em atmosferas)
% 1013: Pressão média da Atmosfera próxima da superfície terrestre (em hectopascais)
rp = P/1013;

% Cálculo da temperatura relativa à atmosfera padrão
% 288: Temperatura média da Atmosfera (em kelvins) 
rt = 288/(T+273);

% Cálculo da pressão parcial do vapor de água saturado (em hectopascais)
ppvas = 6.1121*exp(17.502*T/(240.97+T));

% Cálculo da pressão parcial do vapor de água (em hectopascais)
ppva = H*ppvas;

% Cálculo da densidade do vapor de água
rho = 216.7 * (ppva/(T+273.3));

% Cáclulo do coeficiente de atenuação provocado pelo vapor de àgua (em dB/km)
gammaw = ((3.27e-2)*rt + (1.67e-3) * rho*(rt^7/rp) + (7.7e-4) * (fGhz^0.5) + (3.79/(((fGhz-22.235)^2) + 9.81*(rp^2) * rt))+((11.73 * rt)/(((fGhz-183.31)^2) + 11.85*(rp^2)*rt))+((4.01*rt)/(((fGhz-325.153)^2) + 10.44*(rp^2) * rt)))*(fGhz^2)*rho*rp*rt*1e-4; 

% Cálculo da atenuação provocada pelo vapor de àgua (em dB)
Atw = gammaw * dKm;

fprintf('Atenuação provocada pelo vapor de àgua: %.5gdB\n', Atw);

% Cálculo da atenuação provocada pelo oxigénio

% Cáclulo do coeficiente de atenuação provocado pelo oxigénio (em dB/km)

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

% Cálculo da atenuação provocada pelo oxigénio (em dB)
Ato = gammao * d;
fprintf('Atenuação provocada pelo oxigénio: %.5gdB\n', Ato);

% Cálculo da atenuação suplementar devida à presença da atmosfera (em dB)
Atg = Ato + Atw;
fprintf('Atenuação suplementar devida à presença da atmosfera: %.5gdB\n\n', Atg);

%% 12. Cálculo da atenuação provocada pela precipitação

fprintf('CÁLCULO DA ATENUAÇÃO PROVOCADA PELA PRECIPITAÇÃO\n\n');

% Cálculo da intensidade de precipitação para uma percentagem de tempo de 0.01 (em milímetros por hora)
if region == 1 % Região: H
    R = 32;
elseif region == 2 % Região: K
    R = 42;
end

% Cálculo do k e do alpha
if polarizacao == 1 % Polarização: Horizontal
    k = kAlphaTable(fGhz, 1); 
    a = kAlphaTable(fGhz, 3);
elseif polarizacao == 2 % Polarização: Vertical    
    k = kAlphaTable(fGhz, 2); 
    a = kAlphaTable(fGhz, 4);
end

% Cáculo do coeficiente de atenuação provocada pela precipitação (em dB/km)
gammar = k * R.^a;

% Cálculo do comprimento eficaz do percurso 
def = d / (1 + (d / 35 *( exp(-0.015 * R) ) ) );

% Cálculo da atenuação provocada pela precipitação não excedida em 0.01% do tempo (relativamente à média anual)
Ar1 = gammar * def;

% Cálculo da atenuação provocada pela precipitação não excedida numa determinada percentagem do tempo (relativamente à média anual)
Arp = Ar1 * 0.12 * p^(-(0.546 + 0.043 * log10(p))); 
fprintf('Atenuação provocada pela precipitação não excedida numa determinada percentagem do tempo (relativamente à média anual): %.5gdB\n', Arp);

% Cálculo da atenuação provocada pela precipitação não excedida numa determinada percentagem do tempo (relativamente ao pior mês)
pm = 0.3 * p^1.15;
Arpm = Ar1 * 0.12 * pm^(-(0.546 + 0.043 * log10(pm))); 
fprintf('Atenuação provocada pela precipitação não excedida numa determinada percentagem do tempo (relativamente ao pior mês): %.5gdB\n\n', Arpm);

%% 13. Modelação dos efeitos refrativos na atmosfera

fprintf('MODELAÇÃO DOS EFEITOS REFRATIVOS NA ATMOSFERA\n\n');

% Cálculo da refratividade
N = (77.6/(T+273))*(P + 4810*(ppva/(T+273)));

% Cálculo do índice de refracção
n = 1+N*1e-6;

% Cálculo do índice de refracção modificado M
h = 0 : 100 : 10e3;
M = N + 1e6.*(h./r0);

figure;
plot(M,h);
grid on;
title('M(h)');
xlabel('M');
ylabel('h[km]');

%% 14. Difracção em Terra Esférica

% Efeito do Solo e Polarização
beta = (1 + 1.6*(k^2) + 0.75*(k^4))/(1 + 4.5*(k^2) + 1.35*(k^4));

% Função x
x = beta * d *((pi/(lambda*(r0^2)))^(1/3));
% Função F(x)
Fx = 11 + 10 * log10(x) - 17.6 * x;

% Função y1
y1 = 2 * beta * het *(((pi^2)/((lambda^2)*r0))^(1/3));

% Função y2
y2 = 2 * beta * hrt *(((pi^2)/((lambda^2)*r0))^(1/3));

% Cálculo 
KH = (((2*pi*r0)/lambda)^(-1/3))*(((terrainTable(terrain, 1)-1)^2)+(60*lambda*terrainTable(terrain, 2))^2)^(-1/4);
if polarizacao == 0
    K = KH;
else
    K = KH*sqrt((terrainTable(terrain, 1)^2)+(60*lambda*terrainTable(terrain, 2))^2);
end

% Função G(y2)
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

% Função G(y1)
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

% Difracção
Ad = -(Fx + Gy1 + Gy2);

%% 16. Perdas por Difração devido ao terreno usando o método de Deygout

if length(perfilTerreno) ~= dKm
    error('Número de pontos definidos no ficheiro excel não pode ser diferente da distância (em quilómetros)')
end

distanceArray = 1 : dKm;
plot(distanceArray, perfilTerreno);

AtenuacaoObstaculos = Deygout(d, het, hrt, lambda, dKm, perfilTerreno, re, 0, 0) + A0;
