% Fais le calcul de la detente dans une turbine, en comprenant la
% correction du a la loi de Baumann.
% Ne fais que des detentes de vapeur d'eau. Gazs issues de fumees
% non-autorisees
% INPUTS :  Etats d'entree 
%           Pression de sortie 
%           Rendement isentropique de la turbine
% OUTPUTS : Etats de sortie
%           Vecteurs IO (input output), vecteur contenant tous les etats
%           intermediaires de la turbine.
function [ t_out, h_out, s_out, tIO, hIO, sIO] = detenteTurb( t_in, p_in, h_in, s_in, p_out, eta)
    sIO = zeros(1,300); % Initialisation des vecteurs
    hIO = sIO;
    tIO = sIO;
    sIO(1) = s_in; % Conditions initiales
    hIO(1) = h_in;
    tIO(1) = t_in;
    pIO = linspace(p_in,p_out,length(sIO)); % On travaille par petits dp
    for i=1:length(sIO)-1
        eta_prime = eta*(1-(1-XSteam('x_ph',pIO(i),hIO(i)))); % Baumann
        s_is = sIO(i);
        h_is = XSteam('h_ps',pIO(i+1),s_is);
        hIO(i+1) = hIO(i) - eta_prime*(hIO(i) - h_is); % calcul du vrai h
        tIO(i+1) = XSteam('T_ph',pIO(i+1),hIO(i+1));
        sIO(i+1) = XSteam('s_ph',pIO(i+1),hIO(i+1));
    end
    
    t_out = tIO(length(tIO));
    s_out = sIO(length(tIO));
    h_out = hIO(length(tIO));


end

