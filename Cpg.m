% Renvoie la chaleur specifique des fumees a T2 ou moyennee entre T2 et T1.
% INPUTS :  - M : vecteur (taille 4) des fractions massiques
%               M(1) : M_O2
%               M(2) : C_CO2
%               M(3) : M_N2
%               M(4) : M_H2O
%           - T2 : Temperature initiale/finale
%           - T1 : Temperature finale/initiale
% L'ordre des temperatures ne change rien.
% N'accepte pas les temperatures <300K
function [Cpg_m] = Cpg(M,T2,T1)
    
    if nargin == 3
        Tvec = linspace(T1,T2,100); % Travail avec un vecteur
    elseif nargin == 2
        Tvec = T2; % Chaleur spec. a un point donne
    end
    CpO2 = janaf('c','O2',Tvec);
    CpCO2 = janaf('c','CO2',Tvec);
    CpN2 = janaf('c','N2',Tvec);
    CpH2O = janaf('c','H2O',Tvec);
    Cpg_vec = M(1)*CpO2 + M(2)*CpCO2 + M(3)*CpN2 + M(4)*CpH2O; % Vecteur ou scalaire selon T
    if nargin == 3
        Cpg_m = sum(Cpg_vec)/length(Tvec);
    elseif nargin == 2
        Cpg_m = Cpg_vec;
    end    
end
