% Fonction qui renvoie une chaleur specifique de l'air a la temp. T
% ou la chaleur spec. moyenne entre T2 et T1. En ajoutant un argument de
% sortie, renvoie la constange des gazs de l'air (R*_a)
%
% Si nargout == 2, retourne aussi la constante R* de l'air [kJ/kg/K]
% Si nargin == 1, retourne la chaleur specifique a la T° demandee
% Si nargin == 2, retourne la chaleur specifique moyennee entre T1 et T2
%
% Si T1 ou T2 est <300K , Calcule le Cpa a 300K (sans erreur notable
% commise, le Cpa ne varie presque pas entre 273.15 et 300K)
function [Cpa_m, Ra] = Cpa(T2,T1)
    
    if T2 < 300
        T2 = 300;
    end
    if nargin > 1
        if T1 < 300
            T1 = 300;
        end
    end
    MmN2 = 0.028;
    MmO2 = 0.032;
    
    Ma = 0.79*MmN2 + 0.21*MmO2; % Masse molaire de l'air 

    %Fraction massique de l'air [-]
    MO2 = MmO2*0.21/Ma;
    MN2 = MmN2*0.79/Ma;
    if nargin ==2
        Tvec = linspace(T1,T2,100);
    elseif nargin ==1
        Tvec = T2;
    end
    Cpa_vec = janaf('c','O2',Tvec)*MO2+janaf('c','N2',Tvec)*MN2;
    if nargin == 2
        Cpa_m = sum(Cpa_vec)/length(Tvec);
    else
        Cpa_m = Cpa_vec;
    end
    
    if nargout > 1
        R = 8.314462*1e-3; % [kJ/mol/K]
        Ra = R/Ma;
    end
end