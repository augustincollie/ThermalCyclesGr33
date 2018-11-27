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
close all;

if nargin<3
    display = 1;
    if nargin<2
        options = struct('nsout',8,'reheat',1,'T_max',525,'T_cond_out',30,'p3_hp',200, ...
            'DrumFlag',1,'eta_mec',0.98,'comb',struct('Tmax',1400,'lambda',1.2,'x',0,'y',4) ...
            ,'T_exhaust',100,'x4',0.89,'T_0',15,'TpinchSub',4,'TpinchEx',10, ...
            'TpinchCond',15,'Tdrum',120,'eta_SiC',0.9,'eta_SiT',[0.9 0.9]);
        if nargin<1
            P_e = 220e3; % [kW] Puissance energetique de l'installation
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
        eta_mec = 0.98;  % [-]
    end
    if isfield(options.comb,'Tmax')
        comb_Tmax = options.comb.Tmax;
    else
        comb_Tmax = 1400;  % p.111
    end
    if isfield(options.comb,'lambda')
        comb_lambda = options.comb.lambda;
    else
        comb_lambda = 1.2;  %
    end
    if isfield(options.comb,'x')
        comb_x = options.comb.x;
    else
        comb_x = 0;  % [O/C]
    end
    if isfield(options.comb,'y')
        comb_y = options.comb.y;
    else
        comb_y = 4;  % [H/C]
    end
    if isfield(options,'T_exhaust')
        T_exhaust = options.T_exhaust;
    else
        T_exhaust = 400;  % A FAIRE
    end
    if isfield(options,'p_3')
        p3 = options.p_3;
    else
        switch reheat
            case 1
                p3 = 0.14*p3_hp; % Livre p85 : rapport optimal
            case 0
                p3 = 0; % Pas de resurchauffe
            case 2
                p3 = sqrt(0.14)*p3_hp; 
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

    if isfield(options,'TpinchSub')
        TpinchSub = options.TpinchSub;
    else
        TpinchSub = 4.0;  % [C]
    end
    if isfield(options,'TpinchEx')
        TpinchEx = options.TpinchEx;
    else
        TpinchEx = 10.0;  % [C]
    end
    if isfield(options,'TpinchSub')
        TpinchCond = options.TpinchCond;
    else
        TpinchCond = 15.0;  % [C]
    end
    if isfield(options,'TDrum')
        TDrum = options.TDrum;
    else
        TDrum = 120.0;  % [C]
    end
    if isfield(options,'eta_SiC')
        eta_SiC = options.eta_SiC;
    else
        eta_SiC = 0.9;  % [-]
    end
    if isfield(options,'eta_SiT')
        if isvector(options.eta_SiT)
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
        T_cond_out = TpinchCond + T_0;  % [C]
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
% L9 : sortie de bache
% L10 : post pompe Pb
p = zeros(1,10); % [p,t,x,s,h,e]
t = zeros(1,10);
x = zeros(1,10); % x = NaN si vapeur surchauffe ou liquide sous-refroidi 
s = zeros(1,10);
h = zeros(1,10);
e = zeros(1,10);

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
if abs(x(6) - 1) <= 1e-4  || x(6) < 0.88
    sprintf('ATTENTION : solution actuelle non-valable car x_sortie_turbine = %f avec les configurations actuelles. \n Nouveau x(6) = x4 = %.2f \n',x(6),x4)
    x(6) = x4;
    h(6) = XSteam('h_px',p(6),x(6));
    s(6) = XSteam('s_ph',p(6),h(6));
    t(6) = XSteam('T_ph',p(6),h(6));
    e(6) = Exergie(h(6),s(6));
    eta_SiT(2) = (h(5)-h(6))/(h(5)-h6s);
    sprintf('Nouveau rendement isentr. de la turbine MP et BP : %.3f',eta_SiT(2))    
end


% SORTIE DE CONDENSEUR
% Calcul Etat 7 - p,t,x connus - etat liquide
s(7) = XSteam('sL_T',t(7));
h(7) = XSteam('hL_T',t(7));
e(7) = Exergie(h(7),s(7));
    
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
    DHres = zeros(1,nsout);

    % Premier soutirage en sortie HP
    tbleed(1,1) = t(4);
    pbleed(1,1) = p(4);
    hbleed(1,1) = h(4);
    sbleed(1,1) = s(4);
    ebleed(1,1) = e(4);
    xbleed(1,1) = x(4);

    % Autres soutirages
    hdiv = (h(5) - h(6))/((nsout - 1) + 1); % "Pas" de difference d'enthalpie

    for i=2:nsout % soutirage i=nsout est en sortie de HP
        hbleed(1,i) = h(5)-(i-1)*hdiv;
        s_is = s(5);
        h_is = h(5)-(h(5)-hbleed(1,i))/eta_SiT(2);
        p_is = XSteam('p_hs',h_is,s_is);
        pbleed(1,i) = p_is; % p = p_is, voir hypotheses
        
        % Reste des variables a l'etat 1
        tbleed(1,i) = XSteam('T_ph',pbleed(1,i),hbleed(1,i));
        sbleed(1,i) = XSteam('s_ph',pbleed(1,i),hbleed(1,i));
        ebleed(1,i) = Exergie(hbleed(1,i),sbleed(1,i));
        xbleed(1,i) = XSteam('x_ph',pbleed(1,i),hbleed(1,i));
    end

    % On definit la pression partout
    pbleed([2 3],:) = ones(2,1)*pbleed(1,:); % isobare
    for i=1:nsout-1
        pbleed(4,i) = pbleed(3,i+1); % detente isenthalpique ds la vanne
    end
    pbleed(4,nsout) = pbleed(3,nsout);
    
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
            ebleed(j,i) = Exergie(hbleed(j,i) , sbleed(j,i));
            xbleed(j,i) = XSteam('x_ph',pbleed(j,i),hbleed(j,i));
        end
        tbleed(5,i) = tbleed(4,i);
        hbleed(5,i) = hbleed(4,i);
        sbleed(5,i) = sbleed(4,i);
        ebleed(5,i) = ebleed(4,i);
        xbleed(5,i) = xbleed(4,i);
        % Recherche de l'emplacement de la bache : le flot principal arrive
        % un peu sous-refroidi p/ a la T de départ (130-150 degC, slide S6)
        if (tbleed(3,i) > 120)
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
    tbleed(4,nsout) = t(7) + TpinchSub; % Livre p.69
    pbleed(4,nsout) = pbleed(3,nsout);
    hbleed(4,nsout) = XSteam('h_pT',pbleed(3,nsout),tbleed(4,nsout));
    pbleed(5,nsout) = p(7);
    tbleed(5,nsout) = t(7);
    hbleed(5,nsout) = hbleed(4,nsout);
    for j=4:5
        sbleed(j,nsout) = XSteam('s_ph',pbleed(j,nsout),hbleed(j,nsout));
        ebleed(j,nsout) = Exergie(hbleed(j,nsout) , sbleed(j,nsout));
        xbleed(j,nsout) = XSteam('x_ph',pbleed(j,nsout),hbleed(j,nsout));
    end    
end
for i=1:nsout
    if i == 1
        DHres(i) = 0;
    else
        DHres(i) = xbleed(4,i-1)*(XSteam('hV_p',pbleed(4,i))-XSteam('hL_p',pbleed(4,i)));
    end
end
% tbleed
% pbleed
% xbleed
% hbleed
% sbleed
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
    e(8) = Exergie(h(8),s(8));


    % 9 : sortie de bache
    t(9) = tbleed(3,bache);
    p(9) = pbleed(3,bache);
    x(9) = 0;
    h(9) = XSteam('h_px',p(9),x(9));
    s(9) = XSteam('s_ph',p(9),h(9));
    e(9) = Exergie(h(9),s(9));
    
    % 10 : Sortie Pb
    p(10) = p(3)/3;
    s10s = s(9);
    h10s = XSteam('h_ps',p(10),s10s);
    h(10) = h(9) + (h10s-h(9))/eta_SiC;
    s(10) = XSteam('s_ph',p(10),h(10));
    t(10) = XSteam('T_ph',p(10),h(10));
    x(10) = XSteam('x_ph',p(10),h(10));
    e(10) = Exergie(h(10),s(10));
    

    % 1 : devant Pa
    t(1) = tbleed(3,1) - TpinchEx;
    p(1) = p(10);
    x(1) = 0;
    h(1) = XSteam('h_pT',p(1),t(1));
    s(1) = XSteam('s_pT',p(1),t(1));
    e(1) = Exergie(h(1),s(1));       
end
% 2 : Pompe alimentaire
p(2) = p(3);
s2s = s(1);
h2s = XSteam('h_ps',p(2),s2s);
h(2) = h(1) + (h2s-h(1))/eta_SiC;
s(2) = XSteam('s_ph',p(2),h(2));
t(2) = XSteam('T_ph',p(2),h(2));
x(2) = XSteam('x_ph',p(2),h(2));
e(2) = Exergie(h(2),s(2));


%% FLUX MASSIQUES
%Calcul des Xflow
if nsout > 1
    n = nsout;
    b = bache;
    DHliq = zeros(n,1); % Var. enthalpie du flux principal

    DHbleed = zeros(n,1); % Var. enthalpie soutirages
    DHres = DHres'; % Residus d'enthalpie "post" soutirage
    for i=1:nsout
        % Pour DHbleed
        DHbleed(i) = hbleed(1,i) - hbleed(3,i);
        
        % Pour DHliq
        if i == n %Subcooler
            DHbleed(i) = hbleed(1,i) - hbleed(3,i);
            DHliq(i) = XSteam('h_pT',p(8),tbleed(3,i) - TpinchEx) - XSteam('h_pT',p(8),tbleed(4,i) - TpinchSub);
            DHliq(i) = DHliq(i) + XSteam('h_pT',p(8),tbleed(4,i) - TpinchSub) - h(8);
        elseif i == b % bache
            DHliq(i) = XSteam('h_pT',p(8),tbleed(3,i)-TpinchEx) - XSteam('h_pT',p(8),tbleed(3,i+1) - TpinchEx);
        elseif i == b-1 % Echangeur apres pompe Pb
            DHliq(i) = XSteam('h_pT',p(10),tbleed(3,i) - TpinchEx) - h(10); 
        elseif i < b-1 % apres bache
            DHliq(i) = XSteam('h_pT',p(10),tbleed(3,i) - TpinchEx) - XSteam('h_pT',p(10),tbleed(3,i+1) - TpinchEx);
        else % avant bache
            DHliq(i) = XSteam('h_pT',p(8),tbleed(3,i) - TpinchEx) - XSteam('h_pT',p(8),tbleed(3,i+1) - TpinchEx); 
        end
    end

    %Creation de la matrice necessaire au systeme lineaire
    DHRES = DHres*ones(1,n);
    DHRES = tril(DHRES,-1);
    DHRES(b+1:n,1:b) = 0;
    DHRES(:,b) = 0;
    DHLIQ = DHliq*ones(1,n);
    DHLIQ(b+1:n,1:b) = 0;
    DHLIQ(b+1:n,b) = 0;
    DHBLEED = diag(DHbleed); 
    A = DHBLEED + DHRES - DHLIQ;
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

if nsout > 0
    
    HP = (h(3)-h(4))*(1+sum(Xflow));
    eHP = (e(3)-e(4))*(1+sum(Xflow));
    MBP = h(5)-h(6) + Xflow(nsout)*(hbleed(1,nsout)-h(6)) + (h(5)-hbleed(1,2))*sum(Xflow(2:nsout));
    eMBP = e(5)-e(6) + Xflow(nsout)*(ebleed(1,nsout)-e(6)) + (e(5)-ebleed(1,2))*sum(Xflow(2:nsout));
    for i=2:nsout-1
        MBP = MBP + (hbleed(1,i)-hbleed(1,i+1))*sum(Xflow(i:nsout));
        eMBP = eMBP + (ebleed(1,i)-ebleed(1,i+1))*sum(Xflow(i:nsout));
    end
    P_Pa = (h(2)-h(1))*(1+sum(Xflow));
    P_Pe = (h(8)-h(7))*(1+sum(Xflow(bache+1:nsout)));
    P_Pb = (h(10)-h(9))*(1+sum(Xflow));
    e_Pa = (e(2)-e(1))*(1+sum(Xflow));
    e_Pe = (e(8)-e(7))*(1+sum(Xflow(bache+1:nsout)));
    e_Pb = (e(10)-e(9))*(1+sum(Xflow));
else
    HP = h(3)-h(4);
    eHP = e(3)-e(4);
    MBP = h(5)-h(6);
    eMBP = e(5)-e(6);
    P_Pa = (h(2)-h(1));
    P_Pe = 0;
    P_Pb = 0;
    e_Pa = (e(2)-e(1));
    e_Pe = 0;
    e_Pb = 0;
end
    
Wm = HP+MBP-P_Pa-P_Pe-P_Pb; % Travail moteur
eWm = eHP+eMBP-e_Pa-e_Pe-e_Pb; % "Exergie motrice"
m_cond = P_e/(Wm*eta_mec);
m_tot = m_cond*(1+sum(Xflow));
MASSFLOW(2) = m_tot;
XMASSFLOW = m_cond*Xflow;


%% Verif

% tbleed
% pbleed
% xbleed
% hbleed
% sbleed
% t
% p
% h
% s
% x
% e
% % 
%  1+sum(Xflow)
%% CHAMBRE DE COMBUSTION
    
    comb_Tmax = comb_Tmax + 273.15;
    T_exhaust = T_exhaust + 273.15;
    T0 = T_0 + 273.15;
    
    % Pouvoir Comburivore
    m_a1 = ((32 + 3.76*28)*(1+(comb_y-2*comb_x)/4)) / (12 + comb_y + 16*comb_x);    
    %Combustible
    [e_c, PCI]= LHV(comb_y,comb_x);
    
    % BONUS : On trouve un lambda convenable a la temperature max de comb.
    function f = lambdaAppropriate(L)
       [mass_O2f, mass_CO2f, mass_N2f, mass_H2Of] =  ComputeMassFraction(L,comb_x,comb_y);
       MassFr = [mass_O2f mass_CO2f mass_N2f mass_H2Of];
       Cpgas = Cpg(MassFr,comb_Tmax,T0);
       f = Cpgas*(comb_Tmax-T0)+(1/(L*m_a1))*(Cpgas*(comb_Tmax-T0)-PCI);
    end
    lambda2 = fsolve(@lambdaAppropriate,1.2,optimset('display','off'));
    comb_lambda = lambda2;
    
    %Evaluation de la composition des fumées
    [mass_O2f, mass_CO2f, mass_N2f, mass_H2Of] = ComputeMassFraction(comb_lambda,comb_x,comb_y);
    MassFr = [mass_O2f mass_CO2f mass_N2f mass_H2Of];
    
    Tin = T_exhaust - TpinchEx; % Temperature d'entree de l'air dans la chambre
    h_comb = Cpg(MassFr,comb_Tmax,T0)*(comb_Tmax-273.15);
    s_comb = CpBiz(comb_Tmax,T0)*log(comb_Tmax/T0);
    h_exhaust = Cpg(MassFr,T_exhaust,T0)*(T_exhaust-273.15);
    s_exhaust = CpBiz(T_exhaust,T0)*log(T_exhaust/T0);
    % MATRICE pour resolution des debits massiques
    % Var = [m_a m_c m_g]'
    VECTEUR = zeros(3,1);
    switch reheat
        case 0
            VECTEUR(1) = m_tot*(h(3)-h(2));
        case 1
            VECTEUR(1) = m_tot*(h(3)-h(2)) + (m_tot-XMASSFLOW(1))*(h(5)-h(4));
    end
    % (Cpa(Tin,T0)*(Tin-T0))
    MATRICE = [0     0                  (h_comb-h_exhaust); ...
              -1    -1                  1 ; ...
               1    -comb_lambda*m_a1   0];
    VECTEUR(2:3) = [0;0];
    m_comb = MATRICE\VECTEUR;
    m_a = m_comb(1); m_c = m_comb(2); m_g = m_comb(3);
    e_f = h_comb - T0*s_comb;
    e_exh = h_exhaust - T0*s_exhaust;

    % Etat de reference de l'air de combustion
    
    Cpa_0 = Cpa(T0,273.15);
    h0 = Cpa_0*(T0-273.15); % Par rapport a 0°C (p.118)
    s0 = Cpa_0*log(T0/273.15);
    h_a = Cpa(Tin,273.15)*(Tin-273.15);
    s_a = Cpa(Tin,273.15)*log(Tin/273.15); %Entropie de l'air à la combustion
    
    %e_r = (h_a - h0)- T0*(s_a - s0); %exergie des reactifs (ou l'exergie du combustible a ete negligee)    
    e_r = 0; % Negligeable
    %e_f = PCI/(comb_lambda*m_a1) - Cpg(MassFr,T_exhaust,T0)*T0*log(1+PCI/((comb_lambda*m_a1 + 1)*Cpg(MassFr,comb_Tmax,T0)*T0));
    %e_exh = Cpg(MassFr,T_exhaust,T0)*(T_exhaust - T0) - Cpg(MassFr,T_exhaust,T0)*T0*log(T_exhaust/T0);
    
%     % PLOT CPG
%     Tvec = linspace(300,1500,1200);
%     plot(Tvec,Cpg(MassFr,Tvec))
    
    %Vecteur MASSFLOW
    MASSFLOW(1) = m_a;
    MASSFLOW(3) = m_c;
    MASSFLOW(4) = m_g;
    
    %Vecteur combustion.fum
    COMBUSTION.fum(1) = mass_O2f*m_g;
    COMBUSTION.fum(2) = mass_N2f*m_g;
    COMBUSTION.fum(3) = mass_CO2f*m_g;
    COMBUSTION.fum(4) = mass_H2Of*m_g;
        
    %Structure Combustion
    COMBUSTION.LHV = PCI;
    COMBUSTION.e_c = e_c;
    COMBUSTION.lambda = comb_lambda;
    COMBUSTION.Cp_g = Cpg(MassFr,T_exhaust);
    
    
%% RENDEMENTS
% Energetique
eta_gen = (m_tot*(h(3)-h(2)) + (m_tot-XMASSFLOW(1))*(h(5)-h(4)))/(m_c*PCI);
eta_cyclen = Wm/((h(3)-h(2))*(1+sum(Xflow)) + (1+sum(Xflow(2:nsout)))*(h(5)-h(4)));
eta_toten = eta_mec*eta_gen*eta_cyclen;    

eta_combex = (e_f-e_r)*m_g/(m_c*e_c);
eta_chemnex = (e_f - e_exh)/(e_f-e_r);
eta_transex = (m_tot*(e(3)-e(2)) + (m_tot-XMASSFLOW(1))*(e(5)-e(4)))/(m_g*(e_f-e_exh));
eta_gex = eta_transex*eta_chemnex*eta_combex;
eta_rotex = Wm/eWm;
eta_cyclex = m_cond*Wm/(m_tot*(e(3)-e(2)) + (m_tot-XMASSFLOW(1))*(e(5)-e(4))); % NON

eta_totex = eta_mec*eta_gex*eta_cyclex;
%% Pertes
% Flux primaire
Pprim_en = m_c*PCI;
Pprim_ex = m_c*e_c;

% Energie
kmec = 1/eta_mec - 1;
%Pgen = (m_c*PCI) - (m_tot*(h(3)-h(2)) + (m_tot-XMASSFLOW(1))*(h(5)-h(4)));
Pgen = h_exhaust*m_g;
Pcyclen = m_cond*(h(6)-h(7));

% Exergie
eta_pex = (e_Pe + e_Pa + e_Pb) / (P_Pe + P_Pa + P_Pb);
eta_tex = (HP + MBP) / (eHP + eMBP);
Protex = (1-eta_pex)*(P_Pe + P_Pa + P_Pb)*m_cond + (1-eta_tex)*(HP+MBP)*m_cond;
Pcyclex = (e(6)-e(7))*m_cond;
Pmec = (P_Pe + P_Pa + P_Pb + HP + MBP)*m_cond*kmec;
Pcombex = m_c*e_c - m_g*e_f;
Pchemnex = e_exh*m_g;
Ptransex = Pprim_ex -P_e-Pchemnex-Pcombex-Pmec-Pcyclex-Protex;
Ptotex = Pprim_ex - P_e;

%% OUTPUT
ETA = [eta_cyclen eta_toten eta_cyclex eta_totex eta_gen eta_gex eta_combex eta_chemnex eta_transex];
DATEN = [Pgen Pmec Pcyclen];
DATEX = [Pmec Ptotex Protex Pcombex Pcyclex Pchemnex Ptransex];

if display ==1
    %%%%%%%%%%%%%%%%%%% linearisation
    % 1 - 2 : Pa
    s12 = linspace(s(1),s(2),10);
    p12 = linspace(p(1),p(2),10);
    ts12 = arrayfun( @(p,s) XSteam('T_ps',p,s),p12,s12);
    hs12 = arrayfun( @(p,s) XSteam('h_ps',p,s),p12,s12);
    
    % 2 - 3 : Chambre combustion
    s23 = linspace(s(2),s(3),400);
    p23 = p(2);
    ts23 = arrayfun( @(s) XSteam('T_ps',p23,s),s23);
    hs23 = arrayfun( @(s) XSteam('h_ps',p23,s),s23);    
    
    % 3 - 4 : Premiere detente
    s34 = zeros(1,100);
    hs34 = s34;
    ts34 = s34;
    s34(1) = s(3);
    hs34(1) = h(3);
    ts34(1) = t(3);
    p34 = linspace(p(3),p(4),length(s34));
    for i=1:length(s34)-1
        s4s = s34(i);
        h4s = XSteam('h_ps',p34(i+1),s4s);
        hs34(i+1) = hs34(i) - eta_SiT(1)*(hs34(i) - h4s);
        ts34(i+1) = XSteam('T_ph',p34(i+1),hs34(i+1));
        s34(i+1) = XSteam('s_ph',p34(i+1),hs34(i+1));
    end
    
    % 4 - 5 : Resurchauffe
    p45 = p(4);
    ts45 = linspace(t(4),t(5),100);
    s45 = arrayfun( @(t) XSteam('s_pT',p45,t),ts45);
    hs45 = arrayfun( @(t) XSteam('h_pT',p45,t),ts45);
    
    % 5 - 6 : Deuxieme detente
    hs56 = hbleed(1,2:nsout);
    s56 = sbleed(1,2:nsout);
    ts56 = tbleed(1,2:nsout);
    
    % 6 - 7 : Condenseur
    s67 = linspace(s(6),s(7),30);
    p67 = p(6);
    hs67 = arrayfun( @(s) XSteam('h_ps',p67,s),s67);
    ts67 = arrayfun( @(s) XSteam('T_ps',p67,s),s67);
    
    % 7 - 8 : Pompe Pe --> Trop petit ! juste 2 points et une ligne
    s78 = [s(7) s(8)];
    hs78 = [h(7) h(8)];
    ts78 = [t(7) t(8)];
    
    % 8 - 9 : Parcours avant la bache
    p89 = p(8);
    ts89 = [t(8) fliplr(tbleed(3,bache+1:nsout) - TpinchEx)];
    s89 = arrayfun( @(t) XSteam('s_pT',p89,t),ts89);
    hs89 = arrayfun( @(t) XSteam('h_pT',p89,t),ts89);
    ts89 = [ts89 t(9)];
    s89 = [s89 s(9)];
    hs89 = [hs89 h(9)];
    
    % 9 - 10 : Pompe Pb
    s90 = linspace(s(9),s(10),10);
    p90 = linspace(p(9),p(10),10);
    ts90 = arrayfun( @(p,s) XSteam('T_ps',p,s),p90,s90);
    hs90 = arrayfun( @(p,s) XSteam('h_ps',p,s),p90,s90);    
    
    % 10 - 1 : Parcours post-bache
    p01 = p(10);
    ts01 = [t(10) fliplr(tbleed(3,1:bache-1) - TpinchEx)];
    s01 = arrayfun( @(t) XSteam('s_pT',p01,t),ts01);
    hs01 = arrayfun( @(t) XSteam('h_pT',p01,t),ts01);
    
    %T-s et h-s graphes
    T = linspace(0,400,400);
    sliq = arrayfun( @(t) XSteam('sL_T',t),T);
    svap = arrayfun( @(t) XSteam('sV_T',t),T);
    hliq = arrayfun( @(t) XSteam('hL_T',t),T);
    hvap = arrayfun( @(t) XSteam('hV_T',t),T);    
    FIG(1) = figure; % diagramme T-s
    plot([sliq svap],[T T],'-k');
    hold on
    %plot(svap,T,'-k');
    plot([s12 s23 s34 s45 s56 s67 s78 s89 s90 s01],[ts12 ts23 ts34 ts45 ts56 ts67 ts78 ts89 ts90 ts01],'r');    
    plot(sbleed,tbleed,'--')
    legend('Courbe de saturation','Fluide principal','Soutirages','Location','NorthWest')
    title('Diagramme T-s de l''installation')
    xlabel('Entropie [kJ/kg/K]')
    ylabel('Température [°C]')
    
    FIG(2) = figure; % diagramme h-s
    plot([sliq svap],[hliq hvap],'-k');
    hold on
    %plot(svap,hvap,'-k');
    plot([s12 s23 s34 s45 s56 s67 s78 s89 s90 s01],[hs12 hs23 hs34 hs45 hs56 hs67 hs78 hs89 hs90 hs01],'r');
    plot(sbleed,hbleed,'--')
    legend('Courbe de saturation','Fluide principal','Soutirages','Location','SouthEast')    
    title('Diagramme h-s de l''installation')
    xlabel('Entropie [kJ/kg/K]')
    ylabel('Enthalpie [kJ/kg]')
    
    
    
    
    % Pie graphes
    FIG(3) = figure; % Energie
    legends = {sprintf('Pertes au générateur \n de vapeur %.1f MW',Pgen/1000), ...
        sprintf('Pertes au condenseur \n %.1f MW',Pcyclen/1000), ...
        sprintf('Pertes mécaniques \n %.1f MW',Pmec/1000), ...
        sprintf('Puissance effective \n %.1f MW',P_e/1000)};
    pie([Pgen Pcyclen Pmec P_e],legends);
    title('Flux d''énergie de l''installation vapeur');
    FIG(4) = figure; % Exergie
    legends = {sprintf('Pertes au condenseur \n %.1f MW',Pcyclex/1000), ...
        sprintf('Irréversibilité de la combustion \n %.1f MW',Pcombex/1000),...
        sprintf('Pertes mécaniques \n %.1f MW',Pmec/1000), ...
        sprintf('Irréversibilités turb. et compr. \n %.1f MW',Protex/1000)...
        sprintf('Irréversibilités transf. de chaleur \n %.1f MW',Ptransex/1000), ...
        sprintf('Pertes à la cheminée \n %.1f MW',Pchemnex/1000), ...
        sprintf('Puissance effective \n %.1f MW',P_e/1000)};
    pie([Pcyclex Pcombex Pmec Protex Ptransex Pchemnex P_e],legends);  
    title('Flux d''exergie de l''installation vapeur');
end



%Retourne le PCI d'un combustible du type CH_yO_x selon les données
%disponibles dans des tables issues du cours LMECA2160 - Combustion and
%fuels.
% OUTPUT : - LHV : PCI du carburant exprimé en [kJ/kg]
%          - e_c : exergie du combustible en [kJ/kg]
function [e_c, LHV] = LHV(y,x)
    
    % Valeurs LHV en kJ/kmol
    LHV_CO = 282400;
    LHV_CH4 = 802400;
    LHV_C = 393400;
    LHV_H2Ovap = 241800;
    
    % LMECA2160 eq (4.87)
    a = x/(1+y/2); % ecriture
    LHV = a*LHV_CO + a*(y/2)*LHV_H2Ovap + (1-a)*((1-y/4)*LHV_C + (y/4)*LHV_CH4);
    Mmol_combustible = 12+y*1+x*16; % [kg/kmol]
    LHV = LHV/Mmol_combustible;
    
    e_c = 1.04*LHV; % "Element de thermodynamique technique", P.Wauters.
end

% Retourne le Cp moyen pour l'entropie. Demande que T2 > T1 (en K)
function Cp = CpBiz(T2,T1)
    if T2 < T1
        error('Erreur dans CpBiz : T2 < T1');
    end
    Tvec = linspace(T2,T1,500);
    Cpgg = Cpg(MassFr,Tvec);
    integral = - trapz(Tvec,Cpgg./Tvec);
    Cp = integral/log(T2/T1);
    
end


% Retourne l'exergie a un etat donne comme une difference avec 
% INPUT =  - h : enthalpie de l'etat [kJ/kg]
%          - s : entropie de l'etat [kJ/kg]
%          - T0 : temperature de reference 
% OUTPUT = - e : exergie de l'etat [kJ/kg], comparee a l'exergie à T_0 °C
function e = Exergie(h , s)
    h_0= XSteam('hL_T',T_0);
    s_0= XSteam('sL_T',T_0);
    e = (h-h_0) - (273.15+T_0)*(s-s_0); 
end

end
