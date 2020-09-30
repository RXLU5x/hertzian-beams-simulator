function [AtenuacaoObstaculos] = Deygout(d, het, hrt, lambda, dKm, perfilTerreno, re, c, flag)

idx = 1;

if(c == 0)% Caso seja a primeira chamada
    for i = 1 : dKm
        if perfilTerreno(i) ~= 0
            obstaclesDistances(idx) = i;
        	obstaclesHeights(idx) = perfilTerreno(i);
            idx = idx + 1;
        end
    end
elseif(c >= 1)% Caso não seja a primeira chamada
    if(flag == 1)% Relativa à primeira ligação
        for i = 1 : c-1
            if perfilTerreno(i) ~= 0
                obstaclesDistances(idx) = i;
                obstaclesHeights(idx) = perfilTerreno(i);
                idx = idx + 1;
            end
        end
    elseif(flag==2)% Relativa à segunda ligação
        for i = c+1 : dKm
            if perfilTerreno(i) ~= 0
                obstaclesDistances(idx) = i-c;
                obstaclesHeights(idx) = perfilTerreno(i);
                idx = idx + 1;
            end
        end
    end
end

dx = 0:d;
y = ( (dx - d/2).^2 )./(2*re);

for j = 1 : length(obstaclesHeights)
    % Alturas equivaletntes dos obstáculos
    obstaclesHeightsEq(j) = obstaclesHeights(j)-y((obstaclesDistances(j)*1e3));
    % Cálculo do v para cada obstáculo
    obstaclesV(j) = (obstaclesHeightsEq(j)-het)*sqrt((2*d)/(lambda*obstaclesDistances(j)*1e3*(d-obstaclesDistances(j)*1e3)));
end


% Obstáculo Principal
imax = find(obstaclesV == max(obstaclesV));

At1 = 0;
At2 = 0;
count1 = 0;
count2 = 0;

% 1ª ligação

% Cálculo do declive da equação da reta que descreve o raio direto
md1 = (obstaclesHeightsEq(imax) - het) / obstaclesDistances(imax)*1e3;

% Cálculo da ordenada na origem
bd1 = het;

AtenuacaoTotalLigacao1 = 0;
for i = 1 : imax
   rd1 = md1 * obstaclesDistances(i)*1e3 + bd1;
   raioRD1 = sqrt ( obstaclesDistances(i)*1e3 * ((obstaclesDistances(imax)*1e3 - obstaclesDistances(i)*1e3) / obstaclesDistances(imax)*1e3) * lambda);
   if obstaclesHeightsEq(i) >= rd1 - raioRD1
        if(obstaclesDistances(i)*1e3 ~= 0)
           count1 = count1 + 1;
        end
        if(count1 > 1)
          At1 = At1 + Deygout(obstaclesDistances(imax)*1e3, obstaclesHeightsEq(i), hrt, lambda, dKm, perfilTerreno, re, obstaclesDistances(imax), 1);
        end
        if(count1 == 1)
            v = (obstaclesHeightsEq(i)-y1((obstaclesDistances(i))*1e3))*sqrt((2*(obstaclesDistances(imax)*1e3))/(lambda*(obstaclesDistances(i))*1e3*(obstaclesDistances(imax)-obstaclesDistances(i)*1e3)));
            AtenuacaoTotalLigacao1 = AtenuacaoTotalLigacao1 + 6.9 + 20*log10(sqrt((v - 0.1)^2 +1) + v - 0.1);
        end
   end
end

% 2ª ligação

% Cálculo do declive da equação da reta que descreve o raio direto
md2 = (hrt - obstaclesHeightsEq(imax)) / (d-obstaclesDistances(imax)*1e3);

% Cálculo da ordenada na origem
bd2 = obstaclesHeightsEq(imax);

dx2 = obstaclesDistances(imax)*1e3:d;
y2 = md2*(d-dx2)+bd2;

AtenuacaoTotalLigacao2 = 0;
for i = imax+1: length(obstaclesHeights)
   rd2 = md2 * obstaclesDistances(i)*1e3 + bd2;
   raioRD2 = sqrt ( (obstaclesDistances(i)*1e3 - obstaclesDistances(imax)*1e3)* ((d - ((obstaclesDistances(i)*1e3 - obstaclesDistances(imax)*1e3)) / d) * lambda));
   if obstaclesHeightsEq(i) >= rd2 - raioRD2
        if(obstaclesDistances(i)*1e3 ~= 0)
           count2 = count2 + 1;
        end
        if(count2 > 1)
          At2 = At2 + Deygout(d-obstaclesDistances(imax)*1e3, obstaclesHeightsEq(imax), hrt, lambda, dKm, perfilTerreno, re, obstaclesDistances(imax), 2);
        end
        if(count2 == 1)
            v = (obstaclesHeightsEq(i)-y2((obstaclesDistances(i)*1e3-obstaclesDistances(imax)*1e3)))*sqrt((2*(d-obstaclesDistances(imax)*1e3))/(lambda*((obstaclesDistances(i)-obstaclesDistances(imax)))*1e3*(d-obstaclesDistances(i)*1e3)));
            AtenuacaoTotalLigacao2 = AtenuacaoTotalLigacao2 + 6.9 + 20*log10(sqrt((v - 0.1)^2 +1) + v - 0.1);
        end
   end
end

AtenuacaoObstaculoPrincipal = 6.9 + 20*log10(sqrt((obstaclesV(imax) - 0.1)^2 +1) + obstaclesV(imax) - 0.1);

AtenuacaoObstaculos = AtenuacaoObstaculoPrincipal + AtenuacaoTotalLigacao2 + AtenuacaoTotalLigacao1 + At1 + At2;

end

