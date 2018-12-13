function [ CpgM ] = CpgMoyen( MassFr, T2, T1 )
% Retourne la chaleur massique moyenne entre deux temperatures comme etant
% integral(Cpg/T)dT / (log(T2/T1)). Cette fonction est approximative, mais
% reflete mieux la chaleur massique que la fonction Cpg.
% La definition de cette chaleur massique se trouve p.30, eq (1.63).
% Utilise la methode d'integrale numerique trapezoidale (trapz in Matlab)

% CONDITIONS : T2 > T1 && fumees resultantes de la combustion 

T = linspace(T1,T2,300);
CPG = Cpg(MassFr,T); % VECTEUR !

CpgM = trapz(T,CPG./T) / log(T2/T1);
end

