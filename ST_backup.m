function [ETA, XMASSFLOW, DATEN, DATEX, DAT, MASSFLOW, COMBUSTION, FIG] = ST4(P_e,options,display)
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
            'drumFlag',1,'eta_mec',0.98,'comb',struct('Tmax',1200,'x',0.2,'y',0.4) ...
            ,'T_exhaust',120,'x4',0.89,'T_0',15,'TpinchSub',4,'TpinchEx',15, ...
            'TpinchCond',15,'Tdrum',120,'eta_SiC',0.9,'eta_SiT',[0.9 0.9]);
        if nargin<1
            P_e = 250e3; % [kW] Puissance energetique de l'installation
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
    if isfield(options,'Tdrum')
        Tdrum = options.Tdrum;
    else
        Tdrum = 120.0;  % [C]
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


t(3) = T_max;
x(7) = 0;
t(7) = T_cond_out;
%t(6) = t(7);
p(7) = XSteam('psat_T',t(7));
p(6) = p(7);

% Calcul Etat 3
p(3) = p3_hp;
h(3) = XSteam('h_pT',p(3),t(3));
s(3) = XSteam('s_pT',p(3),t(3));
x(3) = XSteam('x_ph',p(3),h(3));
e(3) = Exergie(h(3),s(3));


% RESURCHAUFFE
if reheat >= 1
    % Calcul Etat 4    
    p(4) = p3;
    [t(4), h(4), s(4),~,~,~] = detenteTurb(t(3),p(3),h(3),s(3),p(4),eta_SiT(1));
    s4s = s(3);
    h4s = XSteam('h_ps',p(4),s4s);
    h(4) = h(3) - eta_SiT(1)*(h(3) - h4s);
    t(4) = XSteam('T_ph',p(4),h(4));
    s(4) = XSteam('s_ph',p(4),h(4));
    x(4) = XSteam('x_ph',p(4),h(4));
    e(4) = Exergie(h(4),s(4));
    p(5) = p(4);
    
    if reheat == 2
        % premiere resurchauffe sur 2
        t3_2 = T_max;
        p3_2 = p(4);
        h3_2 = XSteam('h_pT',p3_2,t3_2);
        s3_2 = XSteam('s_pT',p3_2,t3_2);
        x3_2 = XSteam('x_ph',p3_2,h3_2);
        e3_2 = Exergie(h3_2,s3_2);
        
        % deuxieme detente        
        p4_2 = sqrt(0.14)*p3;
        [t4_2, h4_2, s4_2,~,~,~] = detenteTurb(t3_2,p3_2,h3_2,s3_2,p4_2,eta_SiT(2));
        x4_2 = XSteam('x_ph',p4_2,h4_2);
        e4_2 = Exergie(h4_2,s4_2);
        
        p(5) = p4_2;
    end
    
    % Calcul Etat 5 : t et p connus
    % Resurchauffe classique ou 2eme resurchauffe
    t(5) = T_max;
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
[t(6), h(6), s(6),~,~,~] = detenteTurb(t(5),p(5),h(5),s(5),p(6),eta_SiT(2));
x(6) = XSteam('x_ph',p(6),h(6));

% Verification du titre
if (abs(x(6) - 1) <= 1e-4)  || (x(6) < 0.88) || (x(6) < x4)
    disp(['ATTENTION : solution actuelle non-valable car x_sortie_turbine = ', num2str(x(6),3) ,' avec les configurations actuelles.']) 
end
while (abs(x(6) - 1) <= 1e-4)  || (x(6) < 0.88) || (x(6) < x4)
    eta_SiT(2) = eta_SiT(2)-0.005;
    [t(6), h(6), s(6),~,~,~] = detenteTurb(t(5),p(5),h(5),s(5),p(6),eta_SiT(2));
    x(6) = XSteam('x_ph',p(6),h(6));
end
disp(['Nouveau rendement isentr. de la turbine MP et BP : ',num2str(eta_SiT(2))])
e(6) = Exergie(h(6),s(6));


% SORTIE DE CONDENSEUR
% Calcul Etat 7 - p,t,x connus - etat liquide
s(7) = XSteam('sL_T',t(7));
h(7) = XSteam('hL_T',t(7));
e(7) = Exergie(h(7),s(7));
    
%% Soutirages 
% On considere quelques hypotheses :
% 1) Il y a toujours un soutirage en sortie de HP (si nsout > 0, reheat > 0).
% 2) Pour les autres soutirages, ils sont repartis de manière
% "equidistants" au niveau enthalpique.
% 3) La détente dans les turbines est isentropique.
% 4) on peut considerer que p_6is = p_6i en sortie des bleeders.
% 5) La sortie d'un échangeur est à l'état de liquide saturé (avant la vanne).
    
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
    hdiv = (h(5) - h(6))/((nsout - reheat) + 1); % "Pas"
    if reheat >= 1
        tbleed(1,1) = t(4);
        pbleed(1,1) = p(4);
        hbleed(1,1) = h(4);
        sbleed(1,1) = s(4);
        ebleed(1,1) = e(4);
        xbleed(1,1) = x(4);
        if reheat == 2 && nsout > 1
            tbleed(1,2) = t4_2;
            pbleed(1,2) = p4_2;
            hbleed(1,2) = h4_2;
            sbleed(1,2) = s4_2;
            ebleed(1,2) = e4_2;
            xbleed(1,2) = x4_2;
        end
    end
    % Autres soutirages    
    for i=(1+reheat):nsout
         hbleed(1,i) = h(5)-(i-reheat)*hdiv;
         pbleed(1,i) = fsolve(@pressionBleeders,p(5),optimset('display','off','TolX',1e-3));
        
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
        if (tbleed(3,i) > Tdrum) && drumFlag == 1
            bache = i; % On garde l'indice en memoire.
        elseif nsout == 1
            bache = 1;
        end
    end
    if drumFlag == 1
        %La bache n'a pas de vanne de detente
        tbleed(4,bache) = tbleed(3,bache);
        sbleed(4,bache) = sbleed(3,bache);
        pbleed(4,bache) = pbleed(3,bache);
        xbleed(4,bache) = xbleed(3,bache);
        ebleed(4,bache) = ebleed(3,bache);
    end
    
    % ----------------------------------- PARTIE POMPES ET BACHE
    % 1 : devant Pa
    t(1) = tbleed(3,1) - TpinchEx;
    % On doit prendre en compte le NPSH ! On suppose un NPSH de 40 [m]
    % arbitraire plus un facteur de securite de 0.5 [m].
    % On rajoute par la suite encore 1 bar de securite, purement
    % gratuit...
    Pin = (1000*9.81*(40 + 0.5) + XSteam('psat_T',t(1))*1e5)/1e5;
    p(1) = Pin + 1;
    h(1) = XSteam('h_pT',p(1),t(1));
    x(1) = XSteam('x_ph',p(1),h(1));
    s(1) = XSteam('s_pT',p(1),t(1));
    e(1) = Exergie(h(1),s(1));
    
    % 8 Pompe en sortie de condenseur
    if drumFlag == 1
        p(8) = pbleed(3,bache);
    else
        p(8) = p(1);
    end
    s8s = s(7);
    h8s = XSteam('h_ps',p(8),s8s);
    h(8) = h(7) + (h8s-h(7))/eta_SiC;
    s(8) = XSteam('s_ph',p(8),h(8));
    t(8) = XSteam('T_ph',p(8),h(8));
    x(8) = XSteam('x_ph',p(8),h(8));
    e(8) = Exergie(h(8),s(8));
    
    if drumFlag == 1
        % 9 : sortie de bache
        t(9) = tbleed(3,bache);
        p(9) = pbleed(3,bache);
        x(9) = 0;
        h(9) = XSteam('h_px',p(9),x(9));
        s(9) = XSteam('s_ph',p(9),h(9));
        e(9) = Exergie(h(9),s(9));

        % 10 : Sortie Pb
        p(10) = p(1);
        s10s = s(9);
        h10s = XSteam('h_ps',p(10),s10s);
        h(10) = h(9) + (h10s-h(9))/eta_SiC;
        s(10) = XSteam('s_ph',p(10),h(10));
        t(10) = XSteam('T_ph',p(10),h(10));
        x(10) = XSteam('x_ph',p(10),h(10));
        e(10) = Exergie(h(10),s(10));       
    end
    
    % -------------------------------------- SUBCOOLER
    tbleed(4,nsout) = t(8) + TpinchSub; % Livre p.69
    pbleed(4,nsout) = pbleed(3,nsout);
    hbleed(4,nsout) = XSteam('h_pT',pbleed(4,nsout),tbleed(4,nsout));
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

%% Pompe alimentaire

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
    h9 = zeros(n+1,1); % Enthalpie du flux principal
    t9 = tbleed(3,:) - TpinchEx;
    t9 = [t9 0];
    for i=1:nsout
        if i == n %Subcooler
            h9(n) = XSteam('h_pT',p(8),t9(n));
        elseif i == b-1 % Echangeur apres pompe Pb
            h9(i) = h(10);  
            t9(i) = t(10);
        elseif i < b-1 % apres bache
            h9(i) = XSteam('h_pT',p(10),t9(i));
        else % avant bache
            h9(i) = XSteam('h_pT',p(8),t9(i));
        end
    end    
    H9 = [];
    T9 = [];
    for i=1:nsout
        % Pour DHbleed
        DHbleed(i) = hbleed(1,i) - hbleed(3,i);
        
        % Pour DHliq
        if i == n %Subcooler
            DHliq(i) = XSteam('h_pT',p(8),tbleed(3,i) - TpinchEx) - h(8);
            H9 = [H9 XSteam('h_pT',p(8),tbleed(3,i) - TpinchEx) h(8)];
            T9 = [T9 tbleed(3,i)-TpinchEx t(8)];
        elseif i == b % bache
            DHliq(i) = XSteam('h_pT',p(8),tbleed(3,i)-TpinchEx) - XSteam('h_pT',p(8),tbleed(3,i+1) - TpinchEx);
            H9 = [H9 XSteam('h_pT',p(8),tbleed(3,i)-TpinchEx) XSteam('h_pT',p(8),tbleed(3,i+1) - TpinchEx)];
            T9 = [T9 tbleed(3,i)-TpinchEx tbleed(3,i+1)-TpinchEx];
        elseif i == b-1 % Echangeur apres pompe Pb
            DHliq(i) = XSteam('h_pT',p(10),tbleed(3,i) - TpinchEx) - h(10); 
            H9 = [H9 XSteam('h_pT',p(10),tbleed(3,i) - TpinchEx) h(10)];
            T9 = [T9 tbleed(3,i)-TpinchEx t(10)];
        elseif i < b-1 % apres bache
            DHliq(i) = XSteam('h_pT',p(10),tbleed(3,i) - TpinchEx) - XSteam('h_pT',p(10),tbleed(3,i+1) - TpinchEx);
            H9 = [H9 XSteam('h_pT',p(10),tbleed(3,i) - TpinchEx) XSteam('h_pT',p(10),tbleed(3,i+1) - TpinchEx)];
            T9 = [T9 tbleed(3,i)-TpinchEx tbleed(3,i+1)-TpinchEx];
        else % avant bache
            DHliq(i) = XSteam('h_pT',p(8),tbleed(3,i) - TpinchEx) - XSteam('h_pT',p(8),tbleed(3,i+1) - TpinchEx); 
            H9 = [H9 XSteam('h_pT',p(8),tbleed(3,i) - TpinchEx) XSteam('h_pT',p(8),tbleed(3,i+1) - TpinchEx)];
            T9 = [T9 tbleed(3,i)-TpinchEx tbleed(3,i+1)-TpinchEx];
        end
    end

    %Creation de la matrice necessaire au systeme lineaire
    DHRES = DHres*ones(1,n);
    DHRES = tril(DHRES,-1);
    DHLIQ = DHliq*ones(1,n);
    DHLIQ(n,:) = DHLIQ(n,:) + (hbleed(3,nsout)-hbleed(4,nsout))*ones(1,n);
    if drumFlag == 1
        DHRES(b+1:n,1:b) = 0;
        DHRES(:,b) = 0;
        DHLIQ(b+1:n,1:b) = 0;
        DHLIQ(b+1:n,b) = 0;
    end
    DHBLEED = diag(DHbleed); 
    A = DHBLEED + DHRES - DHLIQ;
    B = DHliq;
    Xflow = A\B;
elseif nsout ==1
    Xflow = (h(2)-h(7))/(hbleed(1,nsout)-h(2));
else 
    %Rankin-Hirn
    Xflow = 0;
end

%% Debit de l'installation
if nsout > 0
    switch reheat
        case 0
        Turb = h(5)-h(6) + Xflow(nsout)*(hbleed(1,nsout)-h(6)) + ...
                (h(5)-hbleed(1,2))*sum(Xflow);
        eTurb = e(5)-e(6) + Xflow(nsout)*(ebleed(1,nsout)-e(6)) + ...
                (e(5)-ebleed(1,2))*sum(Xflow);
        case 1
        Turb = h(5)-h(6) + Xflow(nsout)*(hbleed(1,nsout)-h(6)) + ...
                (h(3)-h(4))*(1+sum(Xflow)) + (h(5)-hbleed(1,2))*sum(Xflow(2:nsout));
        eTurb = e(5)-e(6) + Xflow(nsout)*(ebleed(1,nsout)-e(6)) + ...
                (e(3)-e(4))*(1+sum(Xflow)) + (e(5)-ebleed(1,2))*sum(Xflow(2:nsout));
        case 2
        Turb = (h(3)-h(4))*(1+sum(Xflow)) + (h3_2-h4_2)*(1+sum(Xflow(2:nsout))) + ...
                h(5)-h(6) + Xflow(nsout)*(hbleed(1,nsout)-h(6)) + (h(5)-hbleed(1,2))*sum(Xflow(3:nsout));
        eTurb = (e(3)-e(4))*(1+sum(Xflow)) + (e3_2-e4_2)*(1+sum(Xflow(2:nsout))) + ...
                e(5)-e(6) + Xflow(nsout)*(ebleed(1,nsout)-e(6)) + (e(5)-ebleed(1,2))*sum(Xflow(3:nsout));
    end

    for i=(1+reheat):nsout-1
        Turb = Turb + (hbleed(1,i)-hbleed(1,i+1))*sum(Xflow(i:nsout));
        eTurb = eTurb + (ebleed(1,i)-ebleed(1,i+1))*sum(Xflow(i:nsout));
    end
    P_Pa = (h(2)-h(1))*(1+sum(Xflow));
    P_Pe = (h(8)-h(7))*(1+sum(Xflow(bache+1:nsout)));
    P_Pb = (h(10)-h(9))*(1+sum(Xflow));
    e_Pa = (e(2)-e(1))*(1+sum(Xflow));
    e_Pe = (e(8)-e(7))*(1+sum(Xflow(bache+1:nsout)));
    e_Pb = (e(10)-e(9))*(1+sum(Xflow));
else
    switch reheat
        case 0
        Turb = (h(5)-h(6));
        eTurb = (e(5)-e(6));
        case 1
        Turb = (h(3)-h(4)) + (h(5)-h(6));
        eTurb = (e(3)-e(4)) + (e(5)-e(6));
        case 2
        Turb = (h(3)-h(4)) + (h3_2-h4_2) + (h(5)-h(6));
        eTurb = (e(3)-e(4)) + (e3_2-e4_2) + (e(5)-e(6));
    end

    P_Pa = (h(2)-h(1));
    P_Pe = 0;
    P_Pb = 0;
    e_Pa = (e(2)-e(1));
    e_Pe = 0;
    e_Pb = 0;
end
Pumps = P_Pa + P_Pe + P_Pb;
ePumps = e_Pa + e_Pe + e_Pb;
Wm = Turb - Pumps;
eWm = eTurb - ePumps;
m_cond = P_e/(Wm*eta_mec);
m_tot = m_cond*(1+sum(Xflow));
MASSFLOW(2) = m_tot;
XMASSFLOW = m_cond*Xflow;

%% CHAMBRE DE COMBUSTION
    
    comb_Tmax = comb_Tmax + 273.15;
    T_exhaust = T_exhaust + 273.15;
    T0 = T_0 + 273.15;
    
    % Pouvoir Comburivore
    m_a1 = ((32 + 3.76*28)*(1+(comb_y-2*comb_x)/4)) / (12 + comb_y + 16*comb_x);    
    %Combustible
    [e_c, PCI]= LHV(comb_y,comb_x);
    
    % On importe lambda OU on le trouve (plus optimal !) si lambda n'est pas
    % donne
    if isfield(options.comb,'lambda')
        comb_lambda = options.comb.lambda;
    else
        comb_lambda = fsolve(@lambdaAppropriate,1.2,optimset('display','off'));
    end
    %Evaluation de la composition des fumées
    [mass_O2f, mass_CO2f, mass_N2f, mass_H2Of] = ComputeMassFraction(comb_lambda,comb_x,comb_y);
    MassFr = [mass_O2f mass_CO2f mass_N2f mass_H2Of];
    
    Tin = T_exhaust - TpinchEx; % Temperature d'entree de l'air dans la chambre

    
    h_comb = Cpg(MassFr,T0,comb_Tmax)*(comb_Tmax-T0);    
    s_comb = Cpg(MassFr,T0,comb_Tmax)*log(comb_Tmax/T0);
    h_exhaust = Cpg(MassFr,T0,T_exhaust)*(T_exhaust-T0);
    s_exhaust = Cpg(MassFr,T0,T_exhaust)*log(T_exhaust/T0);
    
    h_a = Cpa(T0,Tin)*(Tin-T0);
    s_a = Cpa(T0,Tin)*log(Tin/T0);
    e_r = h_a-T0*s_a;
    
    % MATRICE pour resolution des debits massiques
    % Var = [m_a m_c m_g]'
    VECTEUR = zeros(3,1);
    switch reheat
        case 0
            VECTEUR(1) = m_tot*(h(3)-h(2));
        case 1
            VECTEUR(1) = m_tot*(h(3)-h(2)) + (m_tot-XMASSFLOW(1))*(h(5)-h(4));
        case 2
            VECTEUR(1) = m_tot*(h(3)-h(2)) + (m_tot-XMASSFLOW(1))*(h3_2-h4_2) + ...
                    (m_tot-sum(XMASSFLOW(1:2)))*(h(5)-h(4));
    end

    MATRICE = [h_a     0                  (h_comb-h_exhaust); ...
              -1    -1                  1 ; ...
               1    -comb_lambda*m_a1   0];
    VECTEUR(2:3) = [0;0];
    m_comb = MATRICE\VECTEUR;
    m_a = m_comb(1); m_c = m_comb(2); m_g = m_comb(3);
    e_f = h_comb - T0*s_comb;
    e_exh = h_exhaust - T0*s_exhaust;
    
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
if reheat == 0
    Qcomb = m_tot*(h(3)-h(2));
    eQcomb = m_tot*(e(3)-e(2));
elseif reheat == 1
    Qcomb = m_tot*(h(3)-h(2)) + (m_tot-XMASSFLOW(1))*(h(5)-h(4));
    eQcomb = m_tot*(e(3)-e(2)) + (m_tot-XMASSFLOW(1))*(e(5)-e(4));
else
    Qcomb = m_tot*(h(3)-h(2)) + (m_tot-XMASSFLOW(1))*(h3_2-h4_2) + (m_tot-sum(XMASSFLOW(1:2)))*(h(5)-h(4));
    eQcomb = m_tot*(e(3)-e(2)) + (m_tot-XMASSFLOW(1))*(e3_2-e4_2) + (m_tot-sum(XMASSFLOW(1:2)))*(e(5)-e(4));
end
eta_gen = Qcomb/(m_c*PCI);
eta_cyclen = (Wm*m_cond)/Qcomb;
eta_toten = eta_mec*eta_gen*eta_cyclen;    

eta_combex = (e_f-e_r)*m_g/(m_c*e_c);
eta_chemnex = (e_f - e_exh)/(e_f-e_r);
eta_transex = eQcomb/(m_g*(e_f-e_exh));
eta_gex = eta_transex*eta_chemnex*eta_combex;
eta_rotex = Wm/eWm;
eta_cyclex = m_cond*Wm/eQcomb;

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
eta_pex = (ePumps) / (Pumps);
eta_tex = (Turb) / (eTurb);
Protex = (1-eta_pex)*(Pumps)*m_cond + (1-eta_tex)*(Turb)*m_cond;
Pcyclex = (e(6)-e(7))*m_cond;
Pmec = (Pumps + Turb)*m_cond*kmec;
Pcombex = m_c*e_c - m_g*e_f;
Pchemnex = e_exh*m_g;
Ptransex = Pprim_ex -P_e-Pchemnex-Pcombex-Pmec-Pcyclex-Protex;
Ptotex = Pprim_ex - P_e;

%% OUTPUT
ETA = [eta_cyclen eta_toten eta_cyclex eta_totex eta_gen eta_gex eta_combex eta_chemnex eta_transex];
DATEN = [Pgen Pmec Pcyclen];
DATEX = [Pmec Ptotex Protex Pcombex Pcyclex Pchemnex Ptransex];
XMASSFLOW = flipud(XMASSFLOW);

if display ==1
    %%%%%%%%%%%%%%%%%%% linearisation
    % 1 - 2 : Pa
    s12 = linspace(s(1),s(2),10);
    p12 = linspace(p(1),p(2),10);
    ts12 = arrayfun( @(p,s) XSteam('T_ps',p,s),p12,s12);
    hs12 = arrayfun( @(p,s) XSteam('h_ps',p,s),p12,s12);
    
    % 2 - 3 : Chambre combustion
    if reheat == 0
        s23 = linspace(s(2),s(5),400);
        p23 = p(2);
        ts23 = arrayfun( @(s) XSteam('T_ps',p23,s),s23);
        hs23 = arrayfun( @(s) XSteam('h_ps',p23,s),s23);
    else
        s23 = linspace(s(2),s(3),400);
        p23 = p(2);
        ts23 = arrayfun( @(s) XSteam('T_ps',p23,s),s23);
        hs23 = arrayfun( @(s) XSteam('h_ps',p23,s),s23);    
    end
    
    if reheat >= 1
        % 3 - 4 : Premiere detente
        [~,~,~,ts34,hs34,s34] = detenteTurb(t(3),p(3),h(3),s(3),p(4),eta_SiT(1));

        if reheat == 2
            % 4 - 3' : Premiere resurchauffe
            p43p = p3_2;
            ts43p = linspace(t(4),t3_2,100);
            s43p = arrayfun( @(t) XSteam('s_pT',p43p,t),ts43p);
            hs43p = arrayfun( @(t) XSteam('h_pT',p43p,t),ts43p);
            
            % 3' - 4' = Deuxieme detente
            [~,~,~,ts3p4p,hs3p4p,s3p4p] = detenteTurb(t3_2,p3_2,h3_2,s3_2,p4_2,eta_SiT(2));
           
            % 4' - 5 : Deuxieme resurchauffe
            p4p5 = p4_2;
            ts4p5 = linspace(t4_2,t(5),100);
            s4p5 = arrayfun( @(t) XSteam('s_pT',p4p5,t),ts4p5);
            hs4p5 = arrayfun( @(t) XSteam('h_pT',p4p5,t),ts4p5);
            
            s45 = [];
            ts45 = [];
            hs45 = [];
        else
            % 4 - 5 : Resurchauffe
            p45 = p(4);
            ts45 = linspace(t(4),t(5),100);
            s45 = arrayfun( @(t) XSteam('s_pT',p45,t),ts45);
            hs45 = arrayfun( @(t) XSteam('h_pT',p45,t),ts45);
            
            s43p = []; ts43p = []; hs43p = [];
            s3p4p = []; ts3p4p = []; hs3p4p = [];
            s4p5 = []; ts4p5 = []; hs4p5 = [];
        end
    else % PAS DE RESURCHAUFFE
        s34 = []; s43p = []; s3p4p = []; s4p5 = []; s45 = [];
        ts34 = []; ts43p = []; ts3p4p = []; ts4p5 = []; ts45 = [];
        hs34 = []; hs43p = []; hs3p4p = []; hs4p5 = []; hs45 = [];
    end
    
    % 5 - 6 : Deuxieme/troisieme detente
    [~,~,~,ts56,hs56,s56] = detenteTurb(t(5),p(5),h(5),s(5),p(6),eta_SiT(2));
    
    % 6 - 7 : Condenseur
    s67 = linspace(s(6),s(7),30);
    p67 = p(6);
    hs67 = arrayfun( @(s) XSteam('h_ps',p67,s),s67);
    ts67 = arrayfun( @(s) XSteam('T_ps',p67,s),s67);
    
    if nsout > 0
        % 7 - 8 : Pompe Pe --> Trop petit ! juste 2 points et une ligne
        s78 = [s(7) s(8)];
        hs78 = [h(7) h(8)];
        ts78 = [t(7) t(8)];

        % 8 - 9 : Parcours avant la bache
        p89 = p(8);
        ts89 = [t(8) fliplr(tbleed(3,bache+1:nsout) - TpinchEx)];
        s89 = arrayfun( @(t) XSteam('s_pT',p89,t),ts89);
        hs89 = arrayfun( @(t) XSteam('h_pT',p89,t),ts89);
        if drumFlag == 1
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
        else
            s90 = []; s01 = [];
            ts90 = []; ts01 = [];
            hs90 = []; hs01 = [];
        end
    else
        s78 = [];
        s89 = [];
        s90 = [];
        s01 = [];
        ts78 = [];
        ts89 = [];
        ts90 = [];
        ts01 = [];
        hs78 = [];
        hs89 = [];
        hs90 = [];
        hs01 = [];
    end
    
    %T-s et h-s graphes
    T = linspace(0,400,400);
    sliq = arrayfun( @(t) XSteam('sL_T',t),T);
    svap = arrayfun( @(t) XSteam('sV_T',t),T);
    hliq = arrayfun( @(t) XSteam('hL_T',t),T);
    hvap = arrayfun( @(t) XSteam('hV_T',t),T);    
    FIG(1) = figure; % diagramme T-s
    plot([sliq svap],[T T],'-k'); % Cloche de saturation
    hold on
    plot([s12 s23 s34 s43p s3p4p s4p5 s45 s56 s67 s78 s89 s90 s01] , [ts12 ts23 ts34 ts43p ts3p4p ts4p5 ts45 ts56 ts67 ts78 ts89 ts90 ts01],'r');    
    if nsout > 0
        plot(sbleed,tbleed,'--')
    end
    legend('Courbe de saturation','Fluide principal','Soutirages','Location','NorthWest')
    title('Diagramme T-s de l''installation')
    xlabel('Entropie [kJ/kg/K]')
    ylabel('Température [°C]')
    
    FIG(2) = figure; % diagramme h-s
    plot([sliq svap],[hliq hvap],'-k'); % Cloche de saturation
    hold on
    plot([s12 s23 s34 s43p s3p4p s4p5 s45 s56 s67 s78 s89 s90 s01] , [hs12 hs23 hs34 hs43p hs3p4p hs4p5 hs45 hs56 hs67 hs78 hs89 hs90 hs01],'r');
    if nsout > 0
        plot(sbleed,hbleed,'--')
    end
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

% Trouve la pression en debut de soutirage
function f = pressionBleeders(ptest)
    [~,hbl,~,~,~,~] = detenteTurb(t(5),p(5),h(5),s(5),ptest,eta_SiT(2));
    f = hbl - hbleed(1,i);
end

% Trouve le coefficient d'exces d'air adequat a la temperature maximale de
% la chambre de combustion et au combustible.
% renvoie une fonction a utiliser avec fsolve.
function f = lambdaAppropriate(L)
   [mass_O2f, mass_CO2f, mass_N2f, mass_H2Of] =  ComputeMassFraction(L,comb_x,comb_y);
   MassFr = [mass_O2f mass_CO2f mass_N2f mass_H2Of];
   Cpgas = Cpg(MassFr,comb_Tmax,T0);
   f = Cpgas*(comb_Tmax-T0)+(1/(L*m_a1))*(Cpgas*(comb_Tmax-T0)-PCI);
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

% Retourne l'exergie POUR L'EAU a un etat donne comme une difference avec 
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


