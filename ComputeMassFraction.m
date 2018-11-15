% Calcule les fractions massiques des fumees, en [kg_espece / kg_fumee]
% Retourne 4 arguments minimum, qui sont les fractions massiques. A la
% demande, un 5eme argument (R*_g) peut etre delivre.
% INPUTS :  - L : Coefficient d'exces d'air (lambda) [/]
%           - x : rapport molaire O/C du combustible [mol_c/mol_c]
%           - y : rapport molaire H/C du combustible [mol_c/mol_c]
% OUTPUTS:  - [MO2, MCO2, MN2, MH2O] : fractions massiques [kg_espece/kg_fumee]
%           - Rg : Constante des gazs des fumees [kJ/kg/K]
function [MO2, MCO2, MN2, MH2O, Rg] = ComputeMassFraction(L,x,y)
    
    MmO2 = 0.032;
    MmN2 = 0.028;
    MmH2O = 0.018;
    MmCO2 = 0.044;
    
    % Fraction molaire
    A = (y-2*x)/4; % Simplicite d'ecriture
    denom = y/2 + (L-1)*(1+A)+3.76*L*(1+A);
    O2 = (L-1)*(1+A)/denom;
    CO2 = 1/denom;
    N2 = 3.76*L*(1+A)/denom;
    H2O = (y/2)/denom;
    
    % Masse molaire des fumees
    M_g = MmO2*O2+MmCO2*CO2+MmN2*N2+MmH2O*H2O; % kg_g/mol_g
    if nargout > 4
        R = 8.314462*1e-3; % [kJ/mol/K]
        Rg = R/M_g; % [kJ/kg/K]
    end
    % fraction massique
    MO2 = MmO2*O2/M_g;
    MCO2 = MmCO2*CO2/M_g;
    MN2 = MmN2*N2/M_g;
    MH2O = MmH2O*H2O/M_g;
    
end