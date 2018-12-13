function [ETA, MASSFLOW, FIG] = CCGT3P(P_eg,options,display)
% CCGT3P is a Combine cycle Gas Turbine with 2 pressure level
% CCGT3P(P_e,options,display) compute the thermodynamics states for a CCGT
% with 3 pressure level (cfr p166 english reference book) including
% combustion, exchanger and cycles. This is done based on several inputs
% (given in OPTION) and based on a given electricity production P_e.
% It returns the main results. It can as well plots graphs if input 
% argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_EG = electrical power output target for gas turbine [kW]
% OPTIONS is a structure containing :
%   -options.T0       [°C] : Reference temperature
%   -options.T_ext    [°C] : External temperature
%   -options.T_STmax  [°C] : maximum temperature on ST cycle
%   -options.eta_mec  [-] : mecanic efficiency of shafts bearings
%   -options.pdrum   [bar]: Drum pressure
%   -options.pmid    [bar]: Intermediary pressure level
%   -options.x7       [-] : Vapor ratio [gaseous/liquid] (titre)
%   -option.eta_SiC   [-] : Isotrenpic efficiency for compression
%   -option.eta_SiT   [-] : Isotrenpic efficiency for compression
%   -options.GT    [struct] : options for Gas turbine (see GT function) 
% DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then 
%          do not plot.
%
%OUPUTS : 
% ETA is a vector with :
%   -eta(1)  : eta_STcyclen, cycle energy efficiency
%   -eta(2)  : eta_GTcyclen, cycle energy efficiency
%   -eta(3)  : eta_toten, overall energy efficiency
%   -eta(4)  : eta_STcyclex, cycle exegy efficiency
%   -eta(5)  : eta_GTcyclex, cycle exegy efficiency
%   -eta(6)  : eta_totex, overall exergie efficiency
%   -eta(7)  : eta_gen, Steam generator energy efficiency
%   -eta(8)  : eta_gex, Steam generator exergy efficiency
%   -eta(9)  : eta_combex, Combustion exergy efficiency
%   -eta(10) : eta_chemex, Chimney exergy efficiency (losses)
%   -eta(11) : eta_transex, Heat exchanger overall exergy efficiency
%   FYI : eta(i) \in [0;1] [-]
% MASSFLOW is a vector containing : 
%   -massflow(1) [kg/s]: water massflow at high pressure turbine inlet
%   -massflow(2) [kg/s]: water massflow at medium pressure turbine inlet
%   -massflow(3) [kg/s]: water massflow at low pressure turbine inlet
%   -massflow(4) [kg/s]: air massflow at gas turbine inlet 
%   -massflow(5) [kg/s]: combustible massflow
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

close all;
%% Recuperation des donnees

if nargin<3
    display=1;
   if nargin<2
       GasTurb=struct('k_mec',0.015,'T_0',15,'T_ext',15,'r',15,'k_cc',0.95, ...
           'T_3',1250,'eta_PiC',0.9,'eta_PiT',0.9);
       options=struct('T_0',15,'T_ext',15,'T_STmax',565,'pdrum',122.8,'pmid',27.3, ...
           'x7',0.95,'eta_SiC',0.9,'eta_SiT',0.9,'GT',GasTurb);
       if nargin<1
           P_eg=283.7*1e3;
       end
   end
end

if isfield(options,'T_ext')
    T_ext = options.T_ext;
else
    T_ext = 15; % [°C]
end
if isfield(options,'T_STmax')
    T_STmax = options.T_STmax;
else
    T_STmax = 565; % [°C]
end
if isfield(options,'pdrum')
    pdrum = options.pdrum;
else
    pdrum = 122.8; % [bar]
end
if isfield(options,'pmid')
    pmid = options.pmid;
else
    pmid = 27.3; % [bar]
end
if isfield(options,'eta_mec')
    eta_mec = options.eta_mec;
else
    eta_mec = 0.98; % [-]
end
if isfield(options,'x7')
    x7 = options.x7;
else
    x7 = 0.95; % [-]
end
if isfield(options,'T_0')
    T_0 = options.T_0;
else
    T_0 = 15; % [°C]
end
if isfield(options,'eta_SiC')
    eta_SiC = options.eta_SiC;
else
    eta_SiC = 0.9; % [-]
end
if isfield(options,'eta_SiT')
    eta_SiT = options.eta_SiT;
else
    eta_SiT = 0.9; % [-]
end
if isfield(options,'GT')
    GasTurb = options.GT;
else
    GasTurb=struct('k_mec',1/eta_mec-1,'T_0',T_0,'T_ext',T_ext,'r',15,'k_cc',0.95, ...
           'T_3',1250,'eta_PiC',0.9,'eta_PiT',0.9);
end


%% Partie Gas Turbine
[ETA_GT,~, DATEX_GT, DAT_GT, MASSFLOW_GT, COMBUSTION_GT,~] = GT(P_eg,GasTurb,0);

% Tcheminee = DAT_GT(1,4) - 273.15;
% if Tcheminee > T_STmax     
%     warning('Attention, température Tcheminee = %.1f [°C] trop haute pour le cycle vapeur!', Tcheminee);
% end

%Recuperation des donnees de la GT
m_a = MASSFLOW_GT(1);   m_c = MASSFLOW_GT(2);   m_g = MASSFLOW_GT(3);
LHV = COMBUSTION_GT.LHV;    e_c = COMBUSTION_GT.e_c;
Tgas_in = DAT_GT(1,4); % EN KELVIN !!!
MassFr = COMBUSTION_GT.fum/m_g; % [O2 CO2 N2 H2O]

%% Steam Turbine : Allocation et determination des etats

% Allocation des vecteurs d'etats
nbre_states = 7; % Nombre d'etats dans le cycle ST
% 1 : Sortie condenseur
% 2 : Sortie pompe alimentaire
% 3 : Entree turbine HP
% 4 : Sortie turbine HP
% 5 : Entree turbine MP
% 6 : Sortie turbine MP - Entree turbine BP
% 7 : Sortie turbine BP
t = zeros(1,nbre_states);
p = zeros(1,nbre_states);
s = zeros(1,nbre_states);
h = zeros(1,nbre_states);
e = zeros(1,nbre_states);

plow = 4; % [bar], donnees de l'enonce.

% --- 3 : ENTREE TURBINE HP ---
t(3) = T_STmax;
p(3) = pdrum;
h(3) = XSteam('h_pt',p(3),t(3));
s(3) = XSteam('s_ph',p(3),h(3));
e(3) = Exergie(h(3),s(3));

% --- 4 : SORTIE TURBINE HP ---
p(4) = pmid;
[t(4),h(4),s(4),ts34,hs34,s34] = detenteTurb(t(3),p(3),h(3),s(3),p(4),eta_SiT);
e(4) = Exergie(h(4),s(4));

% --- 5 : ENTREE TURBINE MP ---
t(5) = T_STmax;
p(5) = pmid;
h(5) = XSteam('h_pt',p(5),t(5));
s(5) = XSteam('s_ph',p(5),h(5));
e(5) = Exergie(h(5),s(5));

% --- 6 : SORTIE TURBINE MP // ENTREE TURBINE BP ---
p(6) = plow;
[t(6),h(6),s(6),ts56,hs56,s56] = detenteTurb(t(5),p(5),h(5),s(5),p(6),eta_SiT);
e(6) = Exergie(h(6),s(6));

% --- 7 : SORTIE TURBINE BP ---
p(7) = fsolve(@compute_p7,0.04,optimset('display','off'));
[t(7),h(7),s(7),ts67,hs67,s67] = detenteTurb(t(6),p(6),h(6),s(6),p(7),eta_SiT);
e(7) = Exergie(h(7),s(7));

% --- 1 : SORTIE CONDENSEUR ---
p(1) = p(7);
h(1) = XSteam('h_px',p(1),0);
t(1) = XSteam('t_ph',p(1),h(1));
s(1) = XSteam('s_ph',p(1),h(1));
e(1) = Exergie(h(1),s(1));

% --- 2 : SORTIE POMPE ALIMENTAIRE ---
p(2) = plow;
s2s = s(1);
h2s = XSteam('h_ps',p(2),s2s);
h(2) = h(1) + (h2s-h(1))/eta_SiC;
s(2) = XSteam('s_ph',p(2),h(2));
t(2) = XSteam('T_ph',p(2),h(2));
e(2) = Exergie(h(2),s(2));


%% Steam Turbine : Etats supplementaires
TpinchEx = 10; % entre 8 et 15°C (p155-158)
T0K = T_0 + 273.15;
% Etats liquide et vapeur sature
h8p = XSteam('hL_p',p(6));
s8p = XSteam('sL_p',p(6));
h8pp = XSteam('hV_p',p(6));
s8pp = XSteam('sV_p',p(6));
T8sat = XSteam('Tsat_p',p(6));
h9p = XSteam('hL_p',p(4));
s9p = XSteam('sL_p',p(4));
h9pp = XSteam('hV_p',p(4));
s9pp = XSteam('sV_p',p(4));
T9sat = XSteam('Tsat_p',p(4));
h10p = XSteam('hL_p',p(3));
s10p = XSteam('sL_p',p(3));
h10pp = XSteam('hV_p',p(3));
s10pp = XSteam('sV_p',p(3));
T10sat = XSteam('Tsat_p',p(3));
h9 = XSteam('h_pt',p(4),t(6));
s9 = XSteam('s_pt',p(4),t(6));
T9 = t(6);
T8 = XSteam('Tsat_p',p(4));
h8 = XSteam('h_pt',p(6),T8);
s8 = XSteam('s_pt',p(6),T8);
e8p = Exergie(h8p,s8p);
e9p = Exergie(h9p,s9p);
e9 = Exergie(h9,s9);

%Etats pour le gaz
% vecteur chemney = [4g HPg MPg BPg 5g]
h_4g = Cpg(MassFr,Tgas_in,T0K)*(Tgas_in - T0K);
s_4g = CpgMoyen(MassFr,Tgas_in,T0K)*log(Tgas_in/T0K);
e_4g = h_4g - T0K*s_4g;
h_HPg = Cpg(MassFr,T10sat + TpinchEx + 273.15,T0K)*(T10sat + TpinchEx - T_0);
h_MPg = Cpg(MassFr,T9sat + TpinchEx + 273.15,T0K)*(T9sat + TpinchEx - T_0);
h_LPg = Cpg(MassFr,T8sat + TpinchEx + 273.15,T0K)*(T8sat + TpinchEx - T_0);

%% Steam Turbine : bilan d'enthalpie des pompes intermediaires

% Pompe LP to MP 
s_P11 = XSteam('sL_p',p(6));
s_is = s_P11;
h_is = XSteam('h_ps',p(4),s_is);
h_P11 = h8p;
h_P12 = h_P11 + (h_is-h_P11)/eta_SiC;
s_P12 = XSteam('s_ph',p(4),h_P12);
T_P12 = XSteam('T_ph',p(4),h_P12);
e_P11 = Exergie(h_P11,s_P11);
e_P12 = Exergie(h_P12,s_P12);

% Pompe MP to HP
s_P21 = XSteam('sL_p',p(4));
s_is = s_P21;
h_is = XSteam('h_ps',p(3),s_is);
h_P21 = h9p;
e_P21 = Exergie(h9p,s_P21);
h_P22 = h_P21 + (h_is-h_P21)/eta_SiC;
s_P22 = XSteam('s_ph',p(3),h_P22);
e_P22 = Exergie(h_P22,s_P22);
T_P22 = XSteam('T_ph',p(3),h_P22);

%% Steam Turbine : Calcul des Massflows par bilan aux echangeurs

% X = [m_HP ; m_MP ; m_LP]
B = m_g*[(h_HPg - h_4g) ; (h_MPg - h_HPg) ; (h_LPg - h_MPg)];
A = [(h(3)-h10p+h(5)-h(4))          (h(5)-h9)               0 ; ...
     (h10p-h9p-(h_P22-h_P21))       (h9-h9p)               (h(6)-h8) ; ...
     (h9p-h8p-(h_P12-h_P11))  (h9p-h8p-(h_P12-h_P11))      (h8-h8p)];
Massflows = A\-B;
m_HP = Massflows(1);    m_MP = Massflows(2);   m_LP = Massflows(3);
m_tot = sum(Massflows);

%% Steam Turbine : Etat des fumees a la sortie de la cheminee

h_5g = -m_tot*(h8p-h(2))/m_g + h_LPg;
f = @(T5g) h_5g - Cpg(MassFr,T5g,T0K)*(T5g-T0K);
T_5g = fsolve(f,400,optimset('display','off'));

s_5g = CpgMoyen(MassFr,T_5g,T0K)*log(T_5g/T0K);
e_5g = h_5g-T0K*s_5g;

%% RENDEMENTS ET PERTES

% Donnees diverses
Pprim_en = m_c*LHV;
Pprim_ex = m_c*e_c;
P_ev = m_HP*(h(3)-h(4)) + (m_HP+m_MP)*(h(5)-h(6)) + m_tot*(h(6)-h(7)) - ...
        m_tot*(h(2)-h(1)) - (h_P12-h_P11)*(m_tot-m_LP) - (h_P22-h_P21)*m_HP; 
P_ev = P_ev*eta_mec;
P_TGV = P_ev + P_eg;
Qcomb_ST = m_g*(h_4g - h_5g);
eQcomb_ST = m_g*(e_4g - e_5g);
Qgen_ST = m_tot*(h8p-h(2)) + m_LP*(h(6)-h8p) + (m_tot-m_LP)*(h9p-h8p - (h_P12-h_P11)) + ...
        m_HP*(h(3)-h9p-(h_P22-h_P21)) + m_MP*(h9 - h9p);
Qgex_ST = m_tot*(e8p-e(2)) + m_LP*(e(6)-e8p) + (m_tot-m_LP)*(e9p-e8p - (e_P12-e_P11)) + ...
        m_HP*(e(3)-e9p-(e_P22-e_P21)) + m_MP*(e9 - e9p); 

% Calculs pour les pertes rotoriques
eta_HP = (h(3)-h(4))/(e(3)-e(4));
eta_MP = (h(5)-h(6))/(e(5)-e(6));
eta_LP = (h(6)-h(7))/(e(6)-e(7));
eta_P1 = (e_P12-e_P11)/(h_P12-h_P11);
eta_P2 = (e_P22-e_P21)/(h_P22-h_P21);
eta_Pa = (e(2)-e(1))/(h(2)-h(1));
P_rotorST = (1-eta_HP)*m_HP*(h(3)-h(4)) + (1-eta_MP)*(m_HP+m_MP)*(h(5)-h(6)) + ...
    (1-eta_LP)*m_tot*(h(6)-h(7)) + (1-eta_P1)*(m_HP+m_MP)*(h_P12-h_P11) + ...
    (1-eta_P2)*m_HP*(h_P22-h_P21) + (1-eta_Pa)*m_tot*(h(2)-h(1));
P_rotorGT = DATEX_GT(2);
P_rotor = P_rotorST+P_rotorGT;

% Pertes    
P_cond = m_tot*(h(7)-h(1));
eP_cond = m_tot*(e(7)-e(1));
P_mecaST = (1-eta_mec)*P_ev;
P_cheminee = h_5g*m_g;
eP_cheminee = e_5g*m_g;
P_combex = DATEX_GT(3);
P_mecaGT = DATEX_GT(1);
P_meca = P_mecaST + P_mecaGT;
P_transf = eQcomb_ST-Qgex_ST;

% Rendements
eta_GTcyclen = ETA_GT(1);
eta_STcyclen = P_ev/Qcomb_ST;
eta_GTcyclex = ETA_GT(3);
eta_STcyclex = P_ev/eQcomb_ST;
eta_gen = Qgen_ST/Qcomb_ST;
eta_transex = Qgex_ST/eQcomb_ST;
eta_toten = P_TGV/(Pprim_en);
eta_totex = P_TGV/(Pprim_ex);
eta_combex = ETA_GT(6);
eta_chemnex = 1-e_5g/e_4g;
eta_gex = eta_combex*eta_chemnex*eta_transex;

%% OUTPUTS

ETA = [eta_STcyclen eta_GTcyclen eta_toten eta_STcyclex eta_GTcyclex ...
    eta_totex eta_gen eta_gex eta_combex eta_chemnex eta_transex];
MASSFLOW = [Massflows' m_a m_c];

%% PLOTS
if display == 1  
    
    % -------------------------- PIE CHARTS
    FIG(4) = figure; % Energie
    legende = {sprintf('Pertes à la cheminée \n %.1f MW',P_cheminee/1000), ...
        sprintf('Pertes au \n condenseur \n %.1f MW',P_cond/1000), ...
        sprintf('Pertes \n mécaniques \n %.1f MW',P_meca/1000), ...
        sprintf('Puissance \n effective ST \n %.1f MW',P_ev/1000), ...
        sprintf('Puissance \n effective GT \n %.1f MW',P_eg/1000)};
    pie([P_cheminee P_cond P_meca P_ev P_eg],legende);
    title('Flux d''énergie de l''installation TGV');
    FIG(5) = figure; % Exergie
    legende = {sprintf('Irréversibilité de \n la combustion \n %.1f MW',P_combex/1000),...
        sprintf('Irréversibilités \n transfert thermique \n %.1f MW',P_transf/1000), ...
        sprintf('Pertes à la cheminée \n %.1f MW',eP_cheminee/1000), ...
        sprintf('Irréversibilités au \n complexe rotorique \n %.1f MW',P_rotor/1000)...
        sprintf('Pertes au \n condenseur \n %.1f MW',eP_cond/1000), ...
        sprintf('Pertes \n mécaniques \n %.1f MW',P_meca/1000), ...
        sprintf('Puissance \n effective ST \n %.1f MW',P_ev/1000)...
        sprintf('Puissance \n effective GT \n %.1f MW',P_eg/1000)};
    pie([P_combex P_transf eP_cheminee P_rotor eP_cond P_meca P_ev P_eg],legende);  
    title('Flux d''exergie de l''installation TGV');
    
    
    % ----------------------  Plot des echangeurs de chaleur
    gasX = [0 (h_4g-h_5g)*m_g/m_HP];
    gasY = [Tgas_in-273.15 T_5g-273.15];
    steamY = [t(3) ; T10sat ; T10sat ; T9sat ; T9sat ; T8sat ; T8sat ; t(2)];
    steamX = [0 ; (h8p-h(2))*m_tot/m_HP ; (h8pp-h8p)*m_LP/m_HP ; ...
              ((h8-h8pp)*m_LP/m_HP + (h9p-h8p)*(m_tot-m_LP)/m_HP - (h_P12-h_P11)*(m_tot-m_LP)/m_HP); ...
              (h9pp - h9p)*m_MP/m_HP ; ...
              ((h(6)-h8)*m_LP/m_HP + (h9-h9pp)*m_MP/m_HP + (h10p-h9p) - (h_P22-h_P21)); ...
              (h10pp-h10p) ; ...
              ((h(3)-h10pp) + (m_MP/m_HP)*(h(5)-h9) + h(5)-h(4))];
    CumulativeSteamX = flipud(cumsum(steamX));
    steamX = sum(steamX) - CumulativeSteamX;
    FIG(1) = figure;
    plot(gasX,gasY);
    hold on;
    plot(steamX,steamY);
    text(0,t(3),'3')
    text(CumulativeSteamX(1),t(2),'2')
    text(0,Tgas_in-273.15,'4g')
    text((h_4g-h_5g)*m_g/m_HP,T_5g-273.15,'5g')
    title('Diagramme des échangeurs de chaleur dans la chaudière')
    legend('gaz','vapeur')
    xlabel('Quantité d''énergie échangée [kJ/kg_{vHP}]')
    ylabel('température [°C]')
    hold off;
    
    % ------------------------- PLOT DES DIAGRAMMES HS // TS
    % 1 - 2
    s12 = [s(1) s(2)];
    hs12 = [h(1) h(2)];
    ts12 = [t(1) t(2)];
    % 2 - 8'
    p28p = p(2);
    s28p = linspace(s(2),s8p,50);
    hs28p = arrayfun( @(s) XSteam('h_ps',p28p,s),s28p);
    ts28p = arrayfun( @(s) XSteam('T_ps',p28p,s),s28p);
    % 8' - post pompe P1
    s8pP1 = [s8p s_P12];
    hs8pP1 = [h8p h_P12];
    ts8pP1 = [T8sat T_P12];
    % Post P1 - 9'
    pP19p = p(4);
    sP19p = linspace(s_P12,s9p,50);
    hsP19p = arrayfun( @(s) XSteam('h_ps',pP19p,s),sP19p);
    tsP19p = arrayfun( @(s) XSteam('T_ps',pP19p,s),sP19p);
    % 9' - 3
    p9p3 = p(3);
    s9p3 = linspace(s9p,s(3),200);
    hs9p3 = arrayfun( @(s) XSteam('h_ps',p9p3,s),s9p3);
    ts9p3 = arrayfun( @(s) XSteam('T_ps',p9p3,s),s9p3);

    % 4 - 5 : resurchauffe
    p45 = p(4);
    s45 = linspace(s(4),s(5),100);
    hs45 = arrayfun( @(s) XSteam('h_ps',p45,s),s45); 
    ts45 = arrayfun( @(s) XSteam('T_ps',p45,s),s45); 

    % 7 - 1 : Condenseur
    s71 = [s(7) s(1)];
    hs71 = [h(7) h(1)];
    ts71 = [t(7) t(1)];
    
    % 9' - 4
    p9p4 = p(4);
    s9p4 = linspace(s9p,s(4));
    hs9p4 = arrayfun( @(s) XSteam('h_ps',p9p4,s),s9p4);
    ts9p4 = arrayfun( @(s) XSteam('T_ps',p9p4,s),s9p4);
    
    % 8' - 6
    p8p6 = p(6);
    s8p6 = linspace(s8p,s(6));
    hs8p6 = arrayfun( @(s) XSteam('h_ps',p8p6,s),s8p6);
    ts8p6 = arrayfun( @(s) XSteam('T_ps',p8p6,s),s8p6);
    
    FIG(2) = figure; % diagramme T-s
    T = linspace(0,400,400);
    sliq = arrayfun( @(t) XSteam('sL_T',t),T);
    svap = arrayfun( @(t) XSteam('sV_T',t),T);
    hliq = arrayfun( @(t) XSteam('hL_T',t),T);
    hvap = arrayfun( @(t) XSteam('hV_T',t),T);    
    plot([sliq svap],[T T],'-k'); % Cloche de saturation
    hold on
    sdraw = [s(1) s(2) s8p s_P12 s9p s_P22 s10p s10pp s(3) s(4) s(5) s(6) s(7) s(1) s9 s8];
    tdraw = [t(1) t(2) T8sat T_P12 T9sat T_P22 T10sat T10sat t(3) t(4) t(5) t(6) t(7) t(1) T9 T8];    
    plot([s12 s28p s8pP1 sP19p s9p3 s34 s45 s56 s67 s71],[ts12 ts28p ts8pP1 tsP19p ts9p3 ts34 ts45 ts56 ts67 ts71],'r')
    plot(sdraw,tdraw,'*')
    plot(s9p4,ts9p4,'r')
    plot(s8p6,ts8p6,'r')
    text([s(1) s(3) s(4) s(5) s(6) s(7) s8 s9],[t(1) t(3) t(4) t(5) t(6) t(7) T8 T9], ...
        ['  1';'  3';'  4';'  5';'  6';'  7';'  8';'  9'])
    hold off
    title('Diagramme T-s de l''installation');
    legend('Courbe de saturation de l''eau','Cycle vapeur','Location','NorthWest');
    xlabel('Entropie [kJ/kg/K]')
    ylabel('Température [°C]')
    
    
    FIG(3) = figure; % diagramme h-s
    plot([sliq svap],[hliq hvap],'-k'); % Cloche de saturation
    hold on
    sdraw = [s(1) s(2) s8p s_P12 s9p s_P22 s10p s10pp s(3) s(4) s(5) s(6) s(7) s(1) s9 s8];
    tdraw = [h(1) h(2) h8p h_P12 h9p h_P22 h10p h10pp h(3) h(4) h(5) h(6) h(7) h(1) h9 h8];
    plot([s12 s28p s8pP1 sP19p s9p3 s34 s45 s56 s67 s71],[hs12 hs28p hs8pP1 hsP19p hs9p3 hs34 hs45 hs56 hs67 hs71],'r')
    plot(sdraw,tdraw,'*')    
    plot(s9p4,hs9p4,'r')
    plot(s8p6,hs8p6,'r')
    text([s(1) s(3) s(4) s(5) s(6) s(7) s8 s9],[h(1) h(3) h(4) h(5) h(6) h(7) h8 h9], ...
        ['  1';'  3';'  4';'  5';'  6';'  7';'  8';'  9'])
    hold off    
    title('Diagramme h-s de l''installation')
    legend('Courbe de saturation de l''eau','Cycle vapeur','Location','NorthWest')
    xlabel('Entropie [kJ/kg/K]')
    ylabel('Enthalpie [kJ/kg]')    
end
%% Fonctions annexes

% Trouve la pression en sortie de turbine BP en fonction du titre voulu a
% la sortie. A combiner avec un fsolve
function f = compute_p7(p)
    s7s = s(6);
    h7s = XSteam('h_ps',p,s7s);
    h(7) = h(6) - eta_SiT*(h(6) - h7s);
    h_test = XSteam('h_px',p,x7);
    f = h(7) - h_test;
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