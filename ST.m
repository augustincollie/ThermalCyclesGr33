function [ETA XMASSFLOW DATEN DATEX DAT MASSFLOW COMBUSTION FIG] = ST(P_e,options,display)
% ST Steam power plants modelisation
% ST(P_e,options,display) compute the thermodynamics states for a Steam
% power plant (combustion, exchanger, cycle) turbine based on several 
% inputs (given in OPTION) and based on a given electricity production P_e.
% It returns the main results. It can as well plots graphs if input 
% argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_E = electrical power output target [kW]
% OPTIONS is a structure containing :
%   -options.nsout     [-] : Number of feed-heating 
%   -options.reheat    [-] : Number of reheating
%   -options.T_max     [°C] : Maximum steam temperature
%   -options.T_cond_out[°C] : Condenseur cold outlet temperature
%   -options.p3_hp     [bar] : Maximum pressure
%   -options.drumFlag  [-] : if =1 then drum if =0 => no drum. 
%   -options.eta_mec   [-] : mecanic efficiency of shafts bearings
%   -options.comb is a structure containing combustion data : 
%       -comb.Tmax     [°C] : maximum combustion temperature
%       -comb.lambda   [-] : air excess
%       -comb.x        [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
%       -comb.y        [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
%   -options.T_exhaust [°C] : Temperature of exhaust gas out of the chimney
%   -options.p_3       [-] : High pressure after last reheating
%   -options.x4        [-] : Vapor ratio [gaseous/liquid] (in french : titre)
%   -options.T_0       [°C] : Reference temperature
%   -options.TpinchSub [°C] : Temperature pinch at the subcooler
%   -options.TpinchEx  [°C] : Temperature pinch at a heat exchanger
%   -options.TpinchCond[°C] : Temperature pinch at condenser 
%   -options.Tdrum     [°C] : minimal drum temperature
%   -option.eta_SiC    [-] : Isotrenpic efficiency for compression
%   -option.eta_SiT    [-] : Isotrenpic efficiency for Turbine. It can be a vector of 2 values :
%             	             eta_SiT(1)=eta_SiT_HP,eta_SiT(2)=eta_SiT_others
% DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then 
%          do not plot.
%
%OUPUTS : 
% ETA is a vector with :
%   -eta(1) : eta_cyclen, cycle energy efficiency
%   -eta(2) : eta_toten, overall energy efficiency
%   -eta(3) : eta_cyclex, cycle exegy efficiency
%   -eta(4) : eta_totex, overall exergie efficiency
%   -eta(5) : eta_gen, Steam generator energy efficiency
%   -eta(6) : eta_gex, Steam generator exergy efficiency
%   -eta(7) : eta_combex, Combustion exergy efficiency
%   -eta(8) : eta_chemex, Chimney exergy efficiency (losses)
%   -eta(9) : eta_transex, Heat exchanger overall exergy efficiency
%   FYI : eta(i) \in [0;1] [-]
% Xmassflow is a vector with each feedheating massflow [kg/s] (respect to figure 
%           2.33, page 91 "Thermal Power Plants" English version).
%           Xmassflow(1) = mass flow at 6_1 etc...
% DATEN is a vector with : 
%   -daten(1) : perte_gen [kW]
%   -daten(2) : perte_mec [kW]
%   -daten(3) : perte_cond [kW]
% DATEX is a vector with :
%   -datex(1) : perte_mec    [kW]
%   -datex(2) : perte_totex  [kW]
%   -datex(3) : perte_rotex  [kW]
%   -datex(4) : perte_combex [kW]
%   -datex(5) : perte_condex [kW]
%   -datex(6) : perte_chemex [kW]
%   -datex(7) : perte_transex[kW]
% DAT is a matrix containing :
% dat = {T_1       , T_2       , ...       , T_6_I,     T_6_II, ... ;  [°C]
%        p_1       , p_2       , ...       , p_6_I,     p_6_II, ... ;  [bar]
%        h_1       , h_2       , ...       , h_6_I,     h_6_II, ... ;  [kJ/kg]
%        s_1       , s_2       , ...       , s_6_I,     s_6_II, ... ;  [kJ/kg/K]
%        e_1       , e_2       , ...       , e_6_I,     e_6_II, ... ;  [kJ/kg]
%        x_1       , x_2       , ...       , x_6_I,     x_6_II, ... ;   };[-]
% MASSFLOW is a vector containing : 
%   -massflow(1) = m_a, air massflow [kg/s]
%   -massflow(2) = m_v, water massflow at 2 [kg/s]
%   -massflow(3) = m_c, combustible massflow [kg/s] 
%   -massflow(4) = m_f, exhaust gas massflow [kg/s]
% 
% COMBUSTION is a structure with :
%   -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
%   -combustion.e_c    : the combustible exergie         [kJ/kg]
%   -combustion.lambda : the air excess                   [-]
%   -combustion.Cp_g   : heat capacity of exhaust gas     [kJ/kg/K]
%   -combustion.fum  : is a vector of the exhaust gas composition :
%       -fum(1) = m_O2f  : massflow of O2 in exhaust gas [kg/s]
%       -fum(2) = m_N2f  : massflow of N2 in exhaust gas [kg/s]
%       -fum(3) = m_CO2f : massflow of CO2 in exhaust gas [kg/s]
%       -fum(4) = m_H2Of : massflow of H2O in exhaust gas [kg/s] 
%
% FIG is a vector of all the figure you plot. Before each figure, define a
% figure environment such as:  
%  "FIG(1) = figure;
%  plot(x,y1);
%  [...]
%   FIG(2) = figure;
%  plot(x,y2);
%  [...]"
%  Your vector FIG will contain all the figure plot during the run of this
%  code (whatever the size of FIG).

%% Conditions initiales


if nargin<3
    display = 1;
    if nargin<2
        options = struct();
        if nargin<1
            P_e = 500e3; % [kW] Puissance énergétique de l'installation
        end
    end
end


% Importation des options
for d = 1 %Useless loop to fold up initial conditions
    if isfield(options,'nsout')
        nsout = options.nsout;
    else
        nsout = 1;  % [-]
    end
    if isfield(options,'reheat')
        reheat = options.reheat;
    else
        reheat = 1;  % [-]
    end
    if isfield(options,'T_max')
        T_max = options.T_max;
    else
        T_max = 525;  % [C]
    end
    % TCOND_OUT ?
    if isfield(options,'p3_hp')
        p3_hp = options.p3_hp;
    else
        p3_hp = 200;  % [bar]
    end
    if isfield(options,'drumFlag')
        drumFlag = options.drumFlag;
    else
        drumFlag = 1;  % [-]
    end
    if isfield(options,'eta_mec')
        eta_mec = options.eta_mec;
    else
        eta_mec = 0.98;  % [-]
    end
    if isfield(options,'comb.Tmax')
        comb.Tmax = options.comb.Tmax;
    else
        comb_Tmax = 0;  % A FAIRE
    end
    if isfield(options,'comb.lambda')
        comb_lambda = options.comb.lambda;
    else
        comb_lambda = 1.1;  % A FAIRE
    end
    if isfield(options,'comb.x')
        comb_x = options.comb.x;
    else
        comb_x = 0;  % [O/C]
    end
    if isfield(options,'comb.y')
        comb_y = options.comb.y;
    else
        comb_y = 4;  % [H/C]
    end
    if isfield(options,'T_exhaust')
        T_exhaust = options.T_exhaust;
    else
        T_exhaust = 0;  % A FAIRE
    end
    if isfield(options,'p3')
        p3 = options.p3;
    else
        switch reheat
            case 1
                p3 = 0.14*p3_hp; % Livre p85 : rapport optimal
            case 0
                p3 = 0; % Pas de resurchauffe
            case 2
                p3 = 0; % A DEFINIR
            otherwise
                error('A maximum of 2 reheats are allowed.');
        end
    end
    if isfield(options,'x4')
        x4 = options.x4;
    else
        x4 = 0.89;  % [-]
    end
    if isfield(options,'T_0')
        T_0 = options.T_0;
    else
        T_0 = 15.0;  % [C]
    end

    if isfield(options,'TPinchSub')
        TPinchSub = options.TPinchSub;
    else
        TPinchSub = 15.0;  % [C] % PAS FINI
    end
    if isfield(options,'TPinchEx')
        TPinchEx = options.TPinchEx;
    else
        TPinchEx = 15.0;  % [C] % PAS FINI
    end
    if isfield(options,'TPinchSub')
        TPinchCond = options.TPinchCond;
    else
        TPinchCond = 15.0;  % [C] % PAS FINI
    end
    if isfield(options,'TDrum')
        TDrum = options.TDrum;
    else
        TDrum = 15.0;  % [C] % PAS FINI
    end


    if isfield(options,'eta_SiC')
        eta_SiC = options.eta_SiC;
    else
        eta_SiC = 0.9;  % [-]
    end
    if isfield(options,'eta_SiT')
        eta_SiT(1) = options.eta_SiT(1);
    else
        eta_SiT = [0.9 0.9]; % [/]
    end
% A REFAIRE
%     if isfield(options,'T_cond_out')
%         T_cond_out = options.T_cond_out;
%     else
%         T_cond_out = TPinchCond + T_0;  % [C]
%     end
end

%% Partie 1 : ch. comb. -> sortie turb. BP -> sortie condenseur

%Vecteurs des etats pour le passage dans les turbines
% L1 & L2 : non definis dans cette partie
% L3 : entree HP
% L4 : sortie HP
% L5 : entree MP apres resurchauffe
% L6 : sortie BP
% L7 : sortie condenseur
p = zeros(1,7); % [p,t,x,s,h,e]
t = p;
x = p; % x = NaN si vapeur surchauffe ou liquide sous-refroidi 
s = p;
h = p;
e = p;

% Initialisation de la matrice avec les donnees initiales
p(3) = p3_hp;
p(5) = p3;
p(4) = p(5);
t(3) = T_max;
t(5) = T_max;
x(7) = 0;
t(7) = T_0 + TPinchCond;
t(6) = T_0 + TPinchCond;
p(7) = XSteam('psat_T',t(7));
p(6) = p(7);

    % Calcul Etat 3
    h(3) = XSteam('h_pT',p(3),t(3));
    s(3) = XSteam('s_pT',p(3),t(3));
    x(3) = XSteam('x_ph',p(3),h(3));
    e(3) = Exergie(h(3),s(3));


    % RESURCHAUFFE
    if reheat == 1
        % Calcul Etat 4
        s4s = s(3);
        h4s = XSteam('h_ps',p(4),s4s);
        h(4) = h(3) - eta_SiT(1)*(h(3) - h4s);
        t(4) = XSteam('T_ph',p(4),h(4));
        s(4) = XSteam('s_ph',p(4),h(4));
        x(4) = XSteam('x_ph',p(4),h(4));
        e(4) = Exergie(h(4),s(4));

        % Calcul Etat 5 : t et p connus
        h(5) = XSteam('h_pT',p(5),t(5));
        s(5) = XSteam('s_ph',p(5),h(5));
        x(5) = XSteam('x_ph',p(5),h(5));        
        e(5) = Exergie(h(5),s(5));

    elseif reheat == 0
        % On zappe les etats 4 et 5
        p([4 5]) = [NaN p(3)];
        t([4 5]) = [NaN t(3)];
        x([4 5]) = [NaN x(3)];
        s([4 5]) = [NaN s(3)];
        h([4 5]) = [NaN h(3)];
        e([4 5]) = [NaN e(3)];
    end

    % SORTIE DE TURBINE BP
    % Calcul Etat 6
    s6s = s(5);
    x6s = XSteam('x_ps',p(6),s6s);
    h6s = XSteam('h_px',p(6),x6s);
    h(6) = h(5) - eta_SiT(2)*(h(5) - h6s);
    x(6) = XSteam('x_ph',p(6),h(6));
    t(6) = XSteam('T_ph',p(6),h(6));
    s(6) = XSteam('s_ph',p(6),h(6));
    e(6) = Exergie(h(6),s(6));
    % Verification du titre
    if ISNAN(x(6))
        error('ERREUR : Le titre en sortie de turbine BP est > 1');
    end
    if x(6) < 0.88
        error('ERREUR : Le titre en sortie de turbine BP est < 0.88');
    end

    % SORTIE DE CONDENSEUR
    % Calcul Etat 7 - p,t,x connus - etat liquide
    s(7) = XSteam('sL_T',t(7));
    h(7) = XSteam('hL_T',t(7));

%% Soutirages : Sortie de cond. -> Sortie de pompe alim.

    switch nsout
        case 0 % Pas de soutirage - court-circuitage comme cycle Rankin-Hirn
        p(1) = p(7);
        t(1) = t(7);
        x(1) = x(7);
        s(1) = s(7);
        h(1) = h(7);
        e(1) = e(7);
        case 1 % Soutirage obligatoire en sortie de HP
            % TO DO
        case 2 % Soutirage en sortie de HP et MP ?? (reheat == 2)
            % TO DO - sot sure !
        otherwise
            % TO DO
    end
    
    % Pompe alimentaire
    p(2) = p(3);
    % TO DO

end


% Retourne l'exergie a un etat donne comme une difference avec 
% INPUT = - etat : vecteur contenant les variables associees a cet etat
% (p,t,x,s,h)
% OUTPUT : - e : exergie de l'etat [kJ/kg], comparee a l'exergie à 15°C
function e = Exergie(h , s)
    h0= XSteam('hL_T',T_0);
    s0= XSteam('sL_T',T_0);
    e = (h-h0) - T_0*(s-s0); 
end
