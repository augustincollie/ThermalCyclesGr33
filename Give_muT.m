function [ mu_T ] = Give_muT( T )

% Cette fonction a pour but de donner le coefficient µ_T defini comme
% (dH/dp) à T constant et qui est necessaire aux calculs relatifs aux pompes
% comme la pompe alimentaire.
% Valeurs extraites du tableau present dans les exercices sur les
% installations motrices.
% INPUT : Temperature (valable de 0.02 a 230) [°C]
% OUTPUT : µ_T(T) [kJ/(kg*bar)]
temperature = [0.02 10:10:230]; % [°C]
muT = [0.1014 0.0972 0.09385 0.09101 0.08847 0.08614 0.08392 0.08175 ...
        0.0796 0.07743 0.07520 0.07288 0.07044 0.06784 0.06504 0.062 ...
        0.05865 0.05494 0.05079 0.0461 0.04073 0.03454 0.0273 0.01871]; % [kJ/kg bar]

p = polyfit(temperature,muT,3);
T_vec = [T^3 T^2 T 1];

mu_T = T_vec*p';

%x = linspace(0,230,500);
%y = polyval(p,x);


end

