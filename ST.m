function [ETA, XMASSFLOW, DATEN, DATEX, DAT, MASSFLOW, COMBUSTION, FIG] = ST(P_e,options,display)
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
            P_e = 500e3; % [kW] Puissance energetique de l'installation
        end
    end
end


% Importation des options
for d = 1 %Useless loop to fold up initial conditions
    if isfield(options,'nsout')
        nsout = options.nsout;
    else
        nsout = 8;  % [-]
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
        eta_mec = 0.99;  % [-]
    end
    if isfield(options,'comb.Tmax')
        comb_Tmax = options.comb.Tmax;
    else
        comb_Tmax = 1400;  % p.111
    end
    if isfield(options,'comb.lambda')
        comb_lambda = options.comb.lambda;
    else
        comb_lambda = 1.1;  %
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
        TPinchSub = 4.0;  % [C]
    end
    if isfield(options,'TPinchEx')
        TPinchEx = options.TPinchEx;
    else
        TPinchEx = 10.0;  % [C]
    end
    if isfield(options,'TPinchSub')
        TPinchCond = options.TPinchCond;
    else
        TPinchCond = 15.0;  % [C]
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
        if isvec(options.eta_SiT)
            eta_SiT(1) = options.eta_SiT(1);
            eta_SiT(2) = options.eta_SiT(2);
        else
            eta_SiT(1) = options.eta_SiT;
            eta_SiT(2) = options.eta_SiT;
        end
    else
        eta_SiT = [0.9 0.9]; % [/]
    end
    if isfield(options,'T_cond_out')
        T_cond_out = options.T_cond_out;
    else
        T_cond_out = TPinchCond + T_0;  % [C]
    end
end

%% Partie 1 : ch. comb. -> sortie turb. BP -> sortie condenseur

% Vecteurs des etats pour le passage dans les turbines
% L1 : avant Pompe alimentaire
% L2 : apres Pa
% L3 : entree HP
% L4 : sortie HP
% L5 : entree MP apres resurchauffe
% L6 : sortie BP
% L7 : sortie condenseur
% L8 : sortie de pompe Pe
% L9 : de post-Pe a la bache
% L10 : sortie de bache
% L11 : post pompe Pb
% L12 : de post Pb a entrée Pa
p = zeros(1,11); % [p,t,x,s,h,e]
t = zeros(1,11);
x = zeros(1,11); % x = NaN si vapeur surchauffe ou liquide sous-refroidi 
s = zeros(1,11);
h = zeros(1,11);
e = zeros(1,11);

% Initialisation des vecteurs avec les donnees initiales
p(3) = p3_hp;
p(5) = p3;
p(4) = p(5);
t(3) = T_max;
t(5) = T_max;
x(7) = 0;
t(7) = T_cond_out;
t(6) = t(7);
p(7) = XSteam('psat_T',t(7));
p(6) = p(7);

% Calcul Etat 3
h(3) = XSteam('h_pT',p(3),t(3));
s(3) = XSteam('s_pT',p(3),t(3));
x(3) = XSteam('x_ph',p(3),h(3));
e(3) = Exergie(h(3),s(3),T_0);


% RESURCHAUFFE
if reheat == 1
    % Calcul Etat 4
    s4s = s(3);
    h4s = XSteam('h_ps',p(4),s4s);
    h(4) = h(3) - eta_SiT(1)*(h(3) - h4s);
    t(4) = XSteam('T_ph',p(4),h(4));
    s(4) = XSteam('s_ph',p(4),h(4));
    x(4) = XSteam('x_ph',p(4),h(4));
    e(4) = Exergie(h(4),s(4),T_0);

    % Calcul Etat 5 : t et p connus
    h(5) = XSteam('h_pT',p(5),t(5));
    s(5) = XSteam('s_ph',p(5),h(5));
    x(5) = XSteam('x_ph',p(5),h(5));        
    e(5) = Exergie(h(5),s(5),T_0);

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
e(6) = Exergie(h(6),s(6),T_0);
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
e(7) = Exergie(h(7),s(7),T_0);
    
%% Soutirages 
% On considere quelques hypotheses :
% 1) Il y a toujours un soutirage en sortie de HP (si nsout > 0).
% 2) Pour les autres soutirages, ils sont repartis de manière
% "equidistants" au niveau enthalpique.
% 3) La détente dans les turbines est isentropique.
% 4) on peut considerer que p_6is = p_6i en sortie des bleeders.
% 5) a la sortie des éventuelles desurchauffes, etat de vapeur saturee.
% 6) On met la bache de degazage lorsque le soutirage i a une température
% superieure a 120°C a l'etat de liquide -> bache à l'indice i.
% 7) le titre à l'entrée de chaque pompe est x = 0.
% 8) La sortie d'un échangeur est à l'état de liquide saturé.
    
if nsout == 0 % Pas de soutirage - court-circuitage comme cycle Rankin-Hirn
    p(1) = p(7);
    t(1) = t(7);
    x(1) = x(7);
    s(1) = s(7);
    h(1) = h(7);
    e(1) = e(7);
else

    % Allocation et init. pour les differentes etapes des soutirages
    % 4 etapes : extraction vap, vapeur sat., liquide sat., post detente vanne
    % 5eme etape pour le subcooler uniquement
    tbleed = zeros(5,nsout);
    pbleed = zeros(5,nsout);
    hbleed = zeros(5,nsout);
    sbleed = zeros(5,nsout);
    ebleed = zeros(5,nsout);
    xbleed = zeros(5,nsout);

    % Premier soutirage en sortie HP
    tbleed(1,nsout) = t(3);
    pbleed(1,nsout) = p(3);
    hbleed(1,nsout) = h(3);
    sbleed(1,nsout) = s(3);
    ebleed(1,nsout) = e(3);
    xbleed(1,nsout) = x(3);

    % Autres soutirages
    hdiv = (h(5) - h(6))/((nsout - 1) + 1); % "Pas" de difference d'enthalpie

    for i=1:(nsout-1) % soutirage i=nsout est en sortie de HP
        hbleed(1,i) = h(6) + hdiv*(i);
        s_is = s(5);
        h_is = (hbleed(1,i) - h(6))/eta_SiT + h(6);
        p_is = XSteam('p_hs',h_is,s_is);
        pbleed(1,i) = p_is; % p = p_is, voir hypotheses

        % Reste des variables a l'etat 1
        tbleed(1,i) = XSteam('T_ph',pbleed(1,i),hbleed(1,i));
        sbleed(1,i) = XSteam('s_ph',pbleed(1,i),hbleed(1,i));
        ebleed(1,i) = Exergie(hbleed(1,i),sbleed(1,i),T_0);
        xbleed(1,i) = XSteam('x_ph',pbleed(1,i),hbleed(1,i));
    end
    % On definit la pression partout
    pbleed([2 3],:) = ones(2,1)*pbleed(1,:); % isobare
    for i=2:nsout
        pbleed(4,i) = pbleed(3,i-1); % detente isenthalpique ds la vanne
    end
    pbleed(4,1) = p(7); % temporaire
    
    %On definit les autres etats
    bache = 0;
    for i=1:nsout
        if hbleed(1,i) > XSteam('hV_p',pbleed(1,i))
            hbleed(2,i) = XSteam('hV_p',pbleed(1,i));
        else
            hbleed(2,i) = hbleed(1,i); %soutirages non-surchauffes
        end
        hbleed(3,i) = XSteam('hL_p',pbleed(1,i));
        hbleed(4,i) = hbleed(3,i); % detente isenthalpique
        for j=2:4
            tbleed(j,i) = XSteam('T_ph',pbleed(j,i),hbleed(j,i));
            sbleed(j,i) = XSteam('s_ph',pbleed(j,i),hbleed(j,i));
            ebleed(j,i) = Exergie(hbleed(j,i) , sbleed(j,i),T_0);
            xbleed(j,i) = XSteam('x_ph',pbleed(j,i),hbleed(j,i));
        end

        % Recherche de l'emplacement de la bache : le flot principal arrive
        % un peu sous-refroidi p/ a la T de départ (130-150 degC, slide S6)
        if (tbleed(3,i) > 120) && (bache == 0)
            bache = i; % On garde l'indice en memoire.
        elseif nsout == 1
            bache = 1;
        end
    end
    
    %La bache n'a pas de vanne de detente
    tbleed(4,bache) = tbleed(3,bache);
    sbleed(4,bache) = sbleed(3,bache);
    pbleed(4,bache) = pbleed(3,bache);
    xbleed(4,bache) = xbleed(3,bache);
    ebleed(4,bache) = ebleed(3,bache);
    
    % On s'occupe du 1er soutirage, qui a le subcooler
    tbleed(4,1) = t(7) + TPinchSub; % Livre p.69
    pbleed(4,1) = pbleed(3,1);
    hbleed(4,1) = XSteam('h_pT',pbleed(3,1),tbleed(3,1));
    pbleed(5,1) = p(7);
    tbleed(5,1) = t(7);
    hbleed(5,1) = hbleed(4,1);
    for j=4:5
        sbleed(j,1) = XSteam('s_ph',pbleed(j,1),hbleed(j,1));
        ebleed(j,1) = Exergie(hbleed(j,1) , sbleed(j,1),T_0);
        xbleed(j,1) = XSteam('x_ph',pbleed(j,1),hbleed(j,1));
    end    
end

%% De condenseur a pompe alimentaire
% Point 8 : sortie de pompe Pe
% vecteur 9 : de post-Pe a la bache
% Point 10 : sortie de bache
% Point 11 : post pompe Pb
% Vecteur 12 : de post Pb a post-desurchauffeurs
% Point 1 : devant Pa
% Point 2 : apres Pa
if nsout > 0
    % 8 Pompe en sortie de condenseur
    p(8) = pbleed(3,bache);
    s8s = s(7);
    h8s = XSteam('h_ps',p(8),s8s);
    h(8) = h(7) + (h8s-h(7))/eta_SiC;
    s(8) = XSteam('s_ph',p(8),h(8));
    t(8) = XSteam('T_ph',p(8),h(8));
    x(8) = XSteam('x_ph',p(8),h(8));
    e(8) = Exergie(h(8),s(8),T_0);


    % 10 : sortie de bache
    t(10) = tbleed(3,bache);
    p(10) = pbleed(3,bache);
    x(10) = 0;
    h(10) = XSteam('h_px',p(10),x(10));
    s(10) = XSteam('s_ph',p(10),h(10));
    e(10) = Exergie(h(10),s(10),T_0);

    % 1 : devant Pa
    tref = sort(tbleed(2,:),'descend');
    t(1) = tref(1) + TpinchEx;
    p(1) = XSteam('psat_T',t(1)); % On suppose liq. sature a l'entree
    x(1) = 0;
    h(1) = XSteam('h_pT',p(1),t(1));
    s(1) = XSteam('s_pT',p(1),t(1));
    e(1) = Exergie(h(1),s(1),T_0);
end
% 2 : Pompe alimentaire
p(2) = p(3);
s2s = s(1);
h2s = XSteam('h_ps',p(2),s2s);
h(2) = h(1) + (h2s-h(1))/eta_SiC;
s(2) = XSteam('s_ph',p(2),h(2));
t(2) = XSteam('T,ph',p(2),h(2));
x(2) = XSteam('x,ph',p(2),h(2));
e(2) = exergy(h(2),s(2));

if nsout > 0
    % 11 : Sortie Pb
    p(11) = p(1);
    s11s = s(10);
    h11s = XSteam('h_ps',p(11),s11s);
    h(11) = h(10) + (h11s-h(10))/eta_SiC;
    s(11) = XSteam('s_ph',p(11),h(11));
    t(11) = XSteam('T,ph',p(11),h(11));
    x(11) = XSteam('x,ph',p(11),h(11));
    e(11) = Exergie(h(11),s(11),T_0);
end

%% FLUX MASSIQUES
%Calcul des Xflow
if nsout > 1
    DHliq = zeros(nsout,1); % Var. enthalpie du flux principal
    DHbleed = zeros(nsout,1); % Var. enthalpie soutirages
    DHres = zeros(nsout,1); % Residus d'enthalpie "post" soutirage
    for i=1:nsout
        DHbleed(i) = hbleed(1,i) - hbleed(3,i);
        DHliq(i) = XSteam('h,pT',pbleed(1,i),tbleed(3,i)-TPinchEx) - XSteam('h_pT',tbleed(1,i)-TPinchEx);
        if i == (bache-1) || i == nsout
            DHres(i) = 0;
        else
            DHres(i) = hbleed(4,i-1) - hbleed(3,i);
        end
    end
    %Creation de la matrice necessaire au systeme lineaire
    DHRES = [DHres(1:nsout-bache)*ones(1,nsout-bache).*triu(nsout-bache,1) zeros(nsout-bache,bache) ; ...
            zeros(bache,nsout-bache) DHres(bache:nsout)*ones(1,bache).*triu(bache,1)];
    DHLIQ = [DHliq(1:nsout-bache)*ones(1,nsout-bache) zeros(nsout-bache,bache) ; DHliq(bache:nsout)*ones(1,nsout)];
    A = diag(DHbleed) + DHRES - DHLIQ;
    b = DHliq;
    Xflow = A\b;
elseif nsout ==1
    Xflow = (h(2)-h(7))/(hbleed(1,nsout)-h(2));
else 
    %Rankin-Hirn
    Xflow = 0;
end

%% Debit de l'installation
% P_E = ((HP*flow + MP et BP * flow variables) - Pe*flow - Pa*flow -
% Pb*flow)*eta_mec
HP = (h(3)-h(4))*(1+sum(Xflow));
eHP = (e(3)-e(4))*(1+sum(Xflow));
MBP = hbleed(1,1) - h(7); % Dernier etape BP
eMBP = ebleed(1,1)-e(7);
MBP = MBP + h(5)-h(6); % Flux principl
eMBP = eMBP + (e(5)-e(6));
for i=2:nsout-1
    MBP = MBP + (hbleed(i)-hbleed(i-1))*(sum(Xflow(1:i))); % Soutirages
    eMBP = eMBP + (ebleed(i)-ebleed(i-1))*(sum(Xflow(1:i)));
end
P_Pa = (h(2)-h(1))*(1+sum(Xflow));
P_Pe = (h(8)-h(7))*(1+sum(Xflow(1:bache-1)));
P_Pb = (h(11)-h(10))*(1+sum(Xflow));
e_Pa = (e(2)-e(1))*(1+sum(Xflow));
e_Pe = (e(8)-e(7))*(1+sum(Xflow(1:bache-1)));
e_Pb = (e(11)-e(10))*(1+sum(Xflow));


Pm = HP+MBP-P_Pa-P_Pe-P_Pb; % Puissance motrice
ePm = eHP+eMBP-e_Pa-e_Pe-e_Pb; % "Exergie motrice"
MASSFLOW(2) = P_e/(Pm*eta_mec);
m_tot = MASSFLOW(2);
XMASSFLOW = MASSFLOW(2)*Xflow;


%% CHAMBRE DE COMBUSTION

    % Sur base de l'etat 2 et de l'etat 3, evaluation du transfert de
    % chaleur entre le combustible (CH4) et l'eau du circuit.
    
    %Données utiles
    Tsat_p2 = XSteam('Tsat_p',p(2));
    h_2prime = XSteam('hL_p',p(2));
    h_2secnd = XSteam('hV_p',p(2));
    hLV_p2 = (h_2secnd - h_2prime);
    s_2prime = XSteam('sL_p',p(2));
    s_2scnd = XSteam('sV_p',p(2));
    combustion.LHV = 50.15; % [MJ/kg] Pouvoir Calorifique Inférieur du Méthane
    
    
    % Masse Molaire des composants [kg/mol]

    M_molO2 = 0.033;
    M_molCO2 =  0.044;
    M_molN2 = 0.028;
    M_molH2O = 0.018;
    
    M_molAir = 0.21*M_molO2 + 0.79*M_molN2;
    
    %Air -- Comburant 
    
    fracMol_O2Air = 0.21*M_mol_O2 / M_molAir;
    fracMol_N2Air = 0.79*M_mol_N2 / M_molAir;
    R_air = 8.3145/M_molAir/1000;
    
    m_a1 = ((32 + 3.76*28)*(1+((comb_y-2*comb_x)/4)))/(12+comb_y+comb_x*16); % Pouvoir Comburivore
    
    %Combustible
    
    M_molComb = 0.012 + 0.001*comb_y + 0.016*comb_x;
        % PCI = ?????????
    
    %Vecteurs des états du flux principal
    t_exch = [t(2),Tsat_p2,Tsat_p2,t(3)];
    p_exch = p(2);
    h_exch = [h(2), h_2prime, h_2secnd, h(3)];
    s_exch = [s(2), s_2prime, s_2secnd, s(3)];
    e_exch = Exergie(h_exch,s_exch,T_0);
    x_exch = [x(2), 0, 1, x(3)];
    
    %Calcul de l'échange de chaleur à l'échangeur (delta_h)
    dh_exch = [h_exch(2)-h_exch(1), hLV_p2, h_exch(3)-h_exch(4)];
    
    %Evaluation de la composition des fumées
    %comb_y , comb_x, comb_lambda, comb_Tmax
        
        %Fractions molaires
    denom = 1+(comb_y/2)+4.76*comb_lambda*(1+((comb_y-2*comb_x)/4))-(1+((comb_y-2*comb_x)/4));
    fracMol_CO2f = (1)/denom;
    fracMol_H2Of = (comb_y/2)/denom;
    fracMol_O2f = ((comb_lambda-1)*(1+((comb_y-2*comb_x)/4)))/denom;
    fracMol_N2f = (3.76*comb_lambda*(1+((comb_y-2*comb_x)/4)))/denom;
        
        %Fractions massiques des fumées
    massFum = fracMol_CO2f*M_molCO2 + fracMol_H2Of*M_molH2O + fracMol_O2f*M_molO2 + fracMol_N2f*M_molN2;
    mass_CO2f = fracMol_CO2f/massFum;
    mass_H2Of = fracMol_H2Of/massFum;
    mass_O2f = fracMol_O2f/massFum;
    mass_N2f = fracMol_N2f/massFum;
    R_fum = 803145/massFum/1000;
    
    
    %Capacité calorifique des fumées? --> Evaluer le transfert de chaleur
    %des fumées vers le liquide chauffé.
    
        %Calcul du Cp moyen pour chaque gaz d'echappement entre 2 et 2'
    
    T_lin1 = linspace(t_exch(1),t_exch(2))
    c_CO2f = janaf('c','CO2', T_lin1 );
    c_H2Of = janaf('c','H2O', T_lin1);
    c_O2f = janaf('c','O2', T_lin1);
    c_N2f = janaf('c','N2', T_lin1);
    
    c_moy_2 = (1/(t_exch(1)-t_exch(2)))*[sum(c_CO2f), sum(c_H2Of), sum(c_O2f), sum(c_N2f)];

        %Estimation de l'enthalpie des fumées entre 2 et 2'
    
        % CONTINUE HERE (ref = lignes 356+) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
    %Calcul du Cp moyen pour chaque gaz d'echappement entre 2' et 2"
    
    c_CO2f = janaf('c','CO2', t_exch(2));
    c_H2Of = janaf('c','H2O', t_exch(2));
    c_O2f = janaf('c','O2', t_exch(2));
    c_N2f = janaf('c','N2', t_exch(2));
    
    cp_LV = [c_CO2f, c_H2Of, c_O2f, c_N2f];
    
    %Calcul du Cp moyen pour chaque gaz d'echappement entre 2" et 3
    
    c_CO2f = janaf('c','CO2', linspace(t_exch(3),t_exch(4)));
    c_H2Of = janaf('c','H2O', linspace(t_exch(3),t_exch(4)));
    c_O2f = janaf('c','O2', linspace(t_exch(3),t_exch(4)));
    c_N2f = janaf('c','N2', linspace(t_exch(3),t_exch(4)));
    
    c_moy_2 = (1/t_exch(3)-t_exch(4))*[sum(c_CO2f), sum(c_H2Of), sum(c_O2f), sum(c_N2f)];
    
    %Enthalpie des gaz en [kg/mol]
    

    

    %       -fum(1) = m_O2f  : massflow of O2 in exhaust gas [kg/s]
%       -fum(2) = m_N2f  : massflow of N2 in exhaust gas [kg/s]
%       -fum(3) = m_CO2f : massflow of CO2 in exhaust gas [kg/s]
%       -fum(4) = m_H2Of : massflow of H2O in exhaust gas [kg/s] 
    
    
    %Evaluation du débit de combustible nécessaire à la combustion
    
    fracMol_CHyOx = fracMol_CO2f; %Stoechiométriquement identiques. 
    
    
%% RENDEMENTS
% Energetique
eta_gen = m_tot*(h(3)-h(2))/(m_comb*PCI);
eta_cyclen = (HP+MBP)/(MBP-(h(5)-h(6))-(h(3)-h(4)) + (h(3)-h(2)));
eta_toten = eta_mec*eta_gen*eta_cyclen;

% Exergetique
eta_combex = ef*(lambda*ma1 + 1)/ec;
eta_chemnex = (ef - ech)/(ef-er);
eta_transex = m_tot*((e(3)-e(2))+(e(6)-e(5)))/(eta_combex*eta_chemnex*m_comb*ec);
eta_gex = eta_transex*eta_chemnex*eta_combex;
eta_rotex = Pm/ePm;
eta_cyclex = eta_rotex*ePm/(m_tot*(e(3)-e(2)));

eta_totex = eta_mec*eta_gex*eta_cyclex;

%% Pertes
Ltot = P_e/eta_totex;
Lgen = (1-eta_gen)*Ltot;
Lcyclen = (1-eta_cyclen)*Ltot;
Lmec = (1-eta_mec)*Ltot;
Ltotex = (1-eta_totex)*Ltot;
Lrotex = (1-eta_rotex)*Ltot;
Lcombex = (1-eta_combex)*Ltot;
Lcyclex = (1-eta_cyclex)*Ltot;
Lchemnex = (1-eta_chemnex)*Ltot;
Ltransex = (1-eta_transex)*Ltot;

%% OUTPUT
ETA = [eta_cyclen eta_toten eta_cyclex eta_totex eta_gen eta_gex eta_combex eta_chemnex eta_transex];
DATEN = [Lgen Lmec Lcyclen];
DATEX = [Lmec Ltotex Lrotex Lcombex Lcyclex Lchemnex Ltransex];

FIG(1) = figure;
pie(eta_toten,eta_mec,eta_gen,eta_cyclen);
legend('Puissance effective','Pertes mécaniques','Pertes à la cheminée','Pertes au condenseur');

FIG(2) = figure;
pie(eta_totex,eta_mec,eta_combex,eta_chemnex,eta_transex,eta_rotex,eta_cyclex);
legend('Puissance effective','Pertes mécaniques','Pertes à la combustion', ...
    'Pertes à la cheminée','Pertes au générateur de vapeur', ...
    'Pertes à la turbine et aux pompes','Pertes au condenseur');


end

% Retourne l'exergie a un etat donne comme une difference avec 
% INPUT =  - h : enthalpie de l'etat [kJ/kg]
%          - s : entropie de l'etat [kJ/kg]
%          - T0 : temperature de reference 
% OUTPUT = - e : exergie de l'etat [kJ/kg], comparee a l'exergie à T_0 °C
function e = Exergie(h , s , T0)
    h0= XSteam('hL_T',T0);
    s0= XSteam('sL_T',T0);
    e = (h-h0) - (273.15+T0)*(s-s0); 
end

%Retourne le PCI d'un combustible du type CH_yO_x selon les données
%disponibles dans des tables issues du cours LMECA2160 - Combustion and
%fuels.
% OUTPUT : - Lhv : PCI du carburant exprimé en [kJ/kmol]
function Lhv = LHV(y,x)
    if (y == 0 && x == 0) %Carbon Graphite
        Lhv = 393400;
    else if (y == 4 && x == 0) % Methane
        Lhv = 802400;
    else if (y == 8/3 && x == 0) % Octane
        Lhv = 2044400;
    else if (y == 0 && x == 1) % Carbon monoxyde
        Lhv = 282400;
end
