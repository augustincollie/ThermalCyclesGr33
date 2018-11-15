function [ETA, DATEN, DATEX, DAT, MASSFLOW, COMBUSTION, FIG] = GT(P_e,options,display)
% GT Gas turbine modelisation GT(P_e,options,display) compute the
% thermodynamics states for a Gas turbine based on several inputs (given in
% OPTION) and based on a given electricity production P_e. It returns the
% main results. It can as well plots graphs if input argument DISPLAY =
% true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_E = electrical power output target [kW] OPTIONS is a structure
% containing :
%   -options.k_mec [-] : Shaft losses -options.T_0   [°C] : Reference
%   temperature -options.T_ext [°C] : External temperature -options.r
%   [-] : Comperssion ratio -options.k_cc  [-] : Coefficient of pressure
%   losses due to combustion
%                        chamber
%   -options.T_3   [°C] : Temperature after combustion (before turbine)
%   -option.eta_PiC[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for compression
%   -option.eta_PiT[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for expansion
%DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then the
%          do not plot.
%
%OUPUTS :
% ETA is a vector with :
%   -eta(1) : eta_cyclen, cycle energy efficiency -eta(2) : eta_toten,
%   overall energy efficiency -eta(3) : eta_cyclex, cycle exegy efficiency
%   -eta(4) : eta_totex, overall exergie efficiency -eta(5) : eta_rotex,
%   compressor-turbine exergy efficiency -eta(6) : eta_combex, Combustion
%   exergy efficiency FYI : eta(i) \in [0;1] [-]
% DATEN is a vector with :
%   -daten(1) : perte_mec [kW] -daten(2) : perte_ech [kW]
% DATEX is a vector with :
%   -datex(1) : perte_mec [kW] -datex(2) : perte_rotex [kW] -datex(3) :
%   perte_combex [kW] -datex(4) : perte_echex  [kW]
% DAT is a matrix containing : dat = {T_1       , T_2       , T_3       ,
% T_4; [°C]
%        p_1       , p_2       , p_3       , p_4; [bar] h_1       , h_2
%        , h_3       , h_4; [kJ/kg] s_1       , s_2       , s_3       ,
%        s_4; [kJ/kg/K] e_1       , e_2       , e_3       , e_4;};[kJ/kg]
% MASSFLOW is a vector containing :
%   -massflow(1) = m_a, air massflow [kg/s] -massflow(2) = m_c, combustible
%   massflow [kg/s] -massflow(3) = m_f, exhaust gas massflow [kg/s]
% 
% COMBUSTION is a structure with :
%   -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
%   -combustion.e_c    : the combuistible exergie         [kJ/kg]
%   -combustion.lambda : the air excess                   [-]
%   -combustion.Cp_g   : heat capacity of exhaust gas at 400 K [kJ/kg/K]
%   -combustion.fum  : is a vector of the exhaust gas composition :
%       -fum(1) = m_O2f  : massflow of O2 in exhaust gas [kg/s] -fum(2) =
%       m_N2f  : massflow of N2 in exhaust gas [kg/s] -fum(3) = m_CO2f :
%       massflow of CO2 in exhaust gas [kg/s] -fum(4) = m_H2Of : massflow
%       of H2O in exhaust gas [kg/s]
%
% FIG is a vector of all the figure you plot. Before each figure, define a
% figure environment such as:
%  "FIG(1) = figure; plot(x,y1); [...]
%   FIG(2) = figure;
%  plot(x,y2); [...]" Your vector FIG will contain all the figure plot
%  during the run of this code (whatever the size of FIG).
%

%% Your Work
% https://che.utah.edu/~geoff/air_poll/lectures_handouts/combustion_basics.pdf
% Slide 13.
close all;
opt=optimset('display','off');
% Exemple of how to use 'nargin' to check your number of inputs
if nargin<3
    display=1;
   if nargin<2
       % Par defaut, prendre l'exemple p.126
       options=struct('k_mec',0.015,'T_0',0,'T_ext',15,'r',18,'k_cc',0.95, ...
           'T_3',1400,'eta_PiC',0.9,'eta_PiT',0.9);
       if nargin<1
           P_e=230e3;%230MW
       end
   end
end


% Prise des donnees dans la structure
for d=1
    if isfield(options,'T_0')
        T0 = options.T_0;
    else
        T0 = 0; % [°C]  
    end
    if isfield(options,'k_mec')
        k_mec = options.k_mec;
    else
        k_mec = 0.015; % [-] 1 a 2% (p.123)  
    end
    if isfield(options,'T_ext')
        T_ext = options.T_ext;
    else
        T_ext = 15; % [°C]  
    end
    if isfield(options,'k_cc')
        k_cc = options.k_cc;
    else
        k_cc = 0.95; % [-] toujours < 1
    end
    if isfield(options,'r')
        r = options.r;
    else
        r = 10; % [-] ratio. compression, voir enonce
    end
    if isfield(options,'T_3')
        T3 = options.T_3;
    else
        T3 = 1050; % [°C] voir enonce  
    end
    if isfield(options,'eta_PiC')
        eta_PiC = options.eta_PiC;
    else
        eta_PiC = 0.9; % [-]  
    end
    if isfield(options,'eta_PiT')
        eta_PiT = options.eta_PiT;
    else
        eta_PiT = 0.9; % [-]  
    end
end


%% Allocation et Init. des variables globales
T = zeros(1,4);
p = zeros(1,4);
h = zeros(1,4);
s = zeros(1,4);
e = zeros(1,4);

LHV = 50150; % Pour methane [kJ/kg] venant de slides S3
y = 4;
x = 0;
ma1 = ((32 + 3.76*28)*(1+(y-2*x)/4)) / (12 + y + 16*x); % LMECA1855 p.133
% kg/mol
MmO2 = 2*16e-3;
MmCO2 = (2*16+12)*1e-3;
MmN2 = 2*14e-3;
MmH2O = (2*1+16)*1e-3;

%% Etat de reference et calculs de base
% Tres bonne approximation que de tout considerer comme gaz ideal ! (p.116)
% Ainsi, Delta_h = Cp|^T2_T1*DeltaT [kJ/kg]
% Et delta_s = Cp|^T2_T1*log(T2/T1) - R*log(p2/p1) [kJ/kg/K]

T0 = T0+273.15;

% Etat de reference
[Cpa_0,Ra] = Cpa(T0);
h0 = Cpa_0*(T0-273.15); % Par rapport a 0°C (p.118)
s0 = Cpa_0*log(T0/273.15);

%% Etat 1 : Arrivee de l'air
T(1) = T_ext+T0; % [K]
%p(1) = 1.01325; % [bar] - pression atmospherique
p(1) = 1;


h(1) = Cpa(T(1))*(T(1)-T0);
s(1) = Cpa(T(1))*log(T(1)/T0);
e(1) = (h(1)-h0) - T0*(s(1)-s0);

%% Etat 1-2 : Compression polytropique
p(2) = r*p(1);

T12 = [T(1) zeros(1,19)];
p12 = linspace(p(1),p(2),20);
h12 = [h(1) zeros(1,19)];
s12 = [s(1) zeros(1,19)];
for i=2:length(T12)
    f = @(T2) T(1)*(p12(i)/p(1))^(Ra / ( eta_PiC*Cpa(T2,T(1)) )) - T2;
    T12(i) = fsolve(f,400,opt);
    h12(i) = h12(i-1) + Cpa(T12(i))*(T12(i)-T12(i-1));
    s12(i) = s12(i-1) + Cpa(T12(i))*log(T12(i)/T12(i-1)) - Ra*log(p12(i)/p12(i-1));
end
T(2) = T12(length(T12));
h(2) = h12(length(h12));
s(2) = s12(length(s12));
e(2) = (h(2)-h0) - T0*(s(2)-s0);
WmC = h(2)-h(1);

%% Etat 2-3 : Combustion (quasi isobare)
p(3) = k_cc*p(2);
T(3) = T3+273.15; % [K]

%Generalement, en cas machine stationnaire, methane (CH4) utilise -- p.111

%exp = r^(Ra / ( eta_PiC*Cpa(T(2),T(1)) ));

% On tente de trouver la valeur de lambda
% On egalise l'equation 3.9 (p115) a l'equation 3.26 (p119)
% f = @(L) Cpa(T(2),T(1))*T(1)*((1+1/(L*ma1)) * T(3)/T(1) * ...
%     Cpg(L,T(3),T(2))/Cpa(T(2),T(1)) - exp) - LHV/(L*ma1);
%f = @(L) (1+LHV/((L*ma1+1)*Cpg(L,T(3),T(2))*T(2)))-(T(3)/T(2));
lambda = fsolve(@ComputeLambda,2.5,opt);
[MO2, MCO2, MN2, MH2O, Rg] = ComputeMassFraction(lambda);
MF = [MO2 MCO2 MN2 MH2O];

% On calcule le reste de l'etat
T23 = linspace(T(2),T(3),100);
p23 = linspace(p(2),p(3),100);
Cpg23 = Cpg(MF,T23);
h23 = [h(2) zeros(1,99)];
s23 = [s(2) zeros(1,99)];
for i=2:length(T23)
    h23(i) = h23(i-1) + Cpg23(i)*(T23(i)-T23(i-1));
    s23(i) = s23(i-1) + Cpg23(i)*log(T23(i)/T23(i-1)) - Rg*log(p23(i)/p23(i-1));
end
h(3) = h23(length(h23));
s(3) = s23(length(s23));
e(3) = (h(3)-h0) - T0*(s(3)-s0);

%% Etat 3-4 : Detente polytropique
p(4) = p(3)/(k_cc*r);

T34 = [T(3) zeros(1,19)];
p34 = linspace(p(3),p(4),20);
h34 = [h(3) zeros(1,19)];
s34 = [s(3) zeros(1,19)];
for i=2:length(T34)
    f = @(T4) T(3)*(p34(i)/p(3))^(eta_PiT*Rg/Cpg(MF,T4,T(3))) - T4; % p.119
    T34(i) = fsolve(f,400,opt);
    h34(i) = h34(i-1) + Cpg(MF,T34(i),T34(i-1))*(T34(i)-T34(i-1));
    s34(i) = s34(i-1) + Cpg(MF,T34(i),T34(i-1))*log(T34(i)/T34(i-1)) - Rg*log(p34(i)/p34(i-1));
end
T(4) = T34(length(T34));
h(4) = h34(length(h34));
s(4) = s34(length(s34));
e(4) = (h(4)-h0) - T0*(s(4)-s0);
WmT = h(3)-h(4);

%% Calcul des debits
%Systeme lineaire 3x3 a partir des eq. de la page 115:
% X = [m_a m_c m_g]'
MAT = [1 -lambda*ma1 0  ;  1 1 -1  ;  -abs(WmC)*(1-k_mec) 0 abs(WmT)*(1-k_mec)];
VEC = [0 0 P_e]';
X = MAT\VEC;
m_a = X(1);
m_c = X(2);
m_g = X(3);
PmT = WmT*m_g;
PmC = WmC*m_a;

%% Rendements et pertes
S288 = 288*(183.1/16); % Entropie standard du combustible [kJ/kg/K]
HHV = 55695; % [kJ/kg] venant de slides S3
ec = HHV-S288;
Pmec = k_mec*(PmT+PmC);
Qcomb = (1+1/(lambda*ma1))*h(3) - h(2);
Wm = (1+1/(lambda*ma1))*(WmT) - (WmC);

eta_cyclen = (Wm) / (Qcomb);
eta_toten = (PmT-PmC)/(m_c*LHV);
eta_cyclex = (PmT-PmC) / (m_g*e(3)-m_a*e(2));
eta_totex = P_e/(m_c*ec);
% Turbine et compresseur
eta_cex = (e(2)-e(1))/(h(2)-h(1));
eta_tex = (h(3)-h(4))/(e(3)-e(4));
eta_rotex = (PmT-PmC)/(m_g*(e(3)-e(4)) - m_a*(e(2)-e(1)));
eta_combex = (m_g*e(3)-m_a*e(2))/(m_c*ec);


Pechen = (Qcomb-Wm)*m_g;
Protex = (1-eta_cex)*PmC + (1-eta_tex)*PmT;
Pcombex = ec*m_c - (m_g*e(3)-m_a*e(2));
Pechex = e(4)*m_g;
%% Calcul des masses
fum = MF*m_g;

%% Renvoi des datas
DAT = [T;p;h;s;e];
ETA = [eta_cyclen eta_toten eta_cyclex eta_totex eta_rotex eta_combex];
DATEN = [Pmec Pechen];
DATEX = [Pmec Protex Pcombex Pechex];
MASSFLOW = X;
COMBUSTION = struct('LHV',LHV,'e_c',ec,'lambda',lambda,'Cp_g',Cpg(MF,400),'fum',fum);

%% Verif'
% T-T0
% p
% lambda
% h
% s
% e
% PmT
% PmC
% X
% COMBUSTION.fum/m_g

%% PLOT
if display == 1
    hdraw = [h12 h23 h34];
    sdraw = [s12 s23 s34];
    Tdraw = [T12 T23 T34];
    FIG(1) = figure; % h-s
    plot(sdraw,hdraw);
    title('Diagramme h-s de la turbine à gaz')
    xlabel('Entropie [kJ/kg/K]')
    ylabel('Enthalpie [kJ/kg]')
    FIG(2) = figure;
    plot(sdraw, Tdraw); % T-s
    title('Diagramme T-s de la turbine à gaz')
    xlabel('Entropie [kJ/kg/K]')
    ylabel('Température [°K]')
    FIG(3) = figure; % Energie
    legend = {sprintf('Pertes a l''échappement \n %.1f MW',Pechen/1000), ...
        sprintf('Pertes mécaniques \n %.1f MW',Pmec/1000), ...
        sprintf('Puissance effective \n %.1f MW',P_e/1000)};
    pie([Pechen Pmec P_e],legend);
    title('Flux d''énergie de la turbine à gaz');
    FIG(4) = figure; % Exergie
    legend = {sprintf('Irréversibilité de la combustion \n %.1f MW',Pcombex/1000),...
        sprintf('Pertes mécaniques \n %.1f MW',Pmec/1000), ...
        sprintf('Pertes à l''échappement \n %.1f MW',Pechex/1000), ...
        sprintf('Irréversibilités turb. et compr. \n %.1f MW',Protex/1000)...
        sprintf('Puissance effective \n %.1f MW',P_e/1000)};
    pie([Pcombex Pmec Pechex Protex P_e],legend);  
    title('Flux d''exergie de la turbine à gaz');
end


%% NESTED FUNCTIONS
% Cpa : Fonction qui renvoie une chaleur specifique de l'air a la temp. T
% et peut aussi donner sa constante des gazs, R*.
% Si nargout == 2, retourne aussi la constante R* de l'air [kJ/kg/K]
% Si nargin == 1, retourne la chaleur specifique a la T° demandee
% Si nargin == 2, retourne la chaleur specifique moyennee entre T1 et T2
function [Cpa_m, Ra] = Cpa(T2,T1)

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

% Renvoie la chaleur specifique des fumees a T2 ou moyennee entre T2 et T1.
% M est un vecteur : M(1) = MO2 // M(2) = MCO2 // M(3) = MN2 // M(4) = MH2O
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

function f = ComputeLambda(L)
    [MO2, MCO2, MN2, MH2O, Rg] = ComputeMassFraction(L);
    MassFr = [MO2 MCO2 MN2 MH2O];
    Cpg = Cpg(MassFr,T(3),T(2));
    f = Cpg*(T(3)-T(2))+(1/(L*ma1))*(Cpg*(T(3)-T0)-LHV);
end

function [MO2, MCO2, MN2, MH2O, Rg] = ComputeMassFraction(L)
    
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
end
