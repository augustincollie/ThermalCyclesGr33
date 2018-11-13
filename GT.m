function [ETA, DATEN, DATEX, DAT, MASSFLOW, COMBUSTION] = GT(P_e,options,display)
% GT Gas turbine modelisation
% GT(P_e,options,display) compute the thermodynamics states for a Gas
% turbine based on several inputs (given in OPTION) and based on a given 
% electricity production P_e. It returns the main results. It can as well
% plots graphs if input argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_E = electrical power output target [kW]
% OPTIONS is a structure containing :
%   -options.k_mec [-] : Shaft losses 
%   -options.T_0   [°C] : Reference temperature
%   -options.T_ext [°C] : External temperature
%   -options.r     [-] : Comperssion ratio
%   -options.k_cc  [-] : Coefficient of pressure losses due to combustion
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
%   -eta(1) : eta_cyclen, cycle energy efficiency
%   -eta(2) : eta_toten, overall energy efficiency
%   -eta(3) : eta_cyclex, cycle exegy efficiency
%   -eta(4) : eta_totex, overall exergie efficiency
%   -eta(5) : eta_rotex, compressor-turbine exergy efficiency
%   -eta(6) : eta_combex, Combustion exergy efficiency
%   FYI : eta(i) \in [0;1] [-]
% DATEN is a vector with : 
%   -daten(1) : perte_mec [kW]
%   -daten(2) : perte_ech [kW]
% DATEX is a vector with :
%   -datex(1) : perte_mec [kW]
%   -datex(2) : perte_rotex [kW]
%   -datex(3) : perte_combex [kW]
%   -datex(4) : perte_echex  [kW]
% DAT is a matrix containing :
% dat = {T_1       , T_2       , T_3       , T_4; [°C]
%        p_1       , p_2       , p_3       , p_4; [bar]
%        h_1       , h_2       , h_3       , h_4; [kJ/kg]
%        s_1       , s_2       , s_3       , s_4; [kJ/kg/K]
%        e_1       , e_2       , e_3       , e_4;};[kJ/kg]
% MASSFLOW is a vector containing : 
%   -massflow(1) = m_a, air massflow [kg/s]
%   -massflow(2) = m_c, combustible massflow [kg/s] 
%   -massflow(3) = m_f, exhaust gas massflow [kg/s]
% 
% COMBUSTION is a structure with :
%   -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
%   -combustion.e_c    : the combuistible exergie         [kJ/kg]
%   -combustion.lambda : the air excess                   [-]
%   -combustion.Cp_g   : heat capacity of exhaust gas at 400 K [kJ/kg/K]
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
%


%% Your Work
% https://che.utah.edu/~geoff/air_poll/lectures_handouts/combustion_basics.pdf
% Slide 13.
opt=optimset('display','off');
% Exemple of how to use 'nargin' to check your number of inputs
if nargin<3
    display=1;
   if nargin<2
       options=struct('k_mec',0.015,'T_0',15,'T_ext',15,'r',10,'k_cc',0.95, ...
           'T_3',1050,'eta_PiC',0.9,'eta_PiT',0.9);
       if nargin<1
           P_e=250e3;%250MW
       end
   end
end


% Prise des donnees dans la structure
for d=1
    if isfield(options,'T_0')
        T0 = options.T_0;
    else
        T0 = 15; % [°C]  
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


%% Allocation
DAT = zeros(5,4); % [T,p,h,s,e]' 
T = zeros(1,4);
p = zeros(1,4);
h = zeros(1,4);
s = zeros(1,4);
e = zeros(1,4);

%% Etat de reference et calculs de base
% Tres bonne approximation que de tout considerer comme gaz ideal ! (p.116)
% Ainsi, Delta_h = Cp|^T2_T1*DeltaT [kJ/kg]
% Et delta_s = Cp|^T2_T1*log(T2/T1) - R*log(p2/p1) [kJ/kg/K]

T0 = T0+273.15;

% Etat de reference
[Cpa_0,Ra] = Cpa(T0);
p0 = 1.01325;
h0 = Cpa_0*(T0-273.15); % Par rapport a 0°C
s0 = Cpa_0*log(T0/273.15);

%% Etat 1 : Arrivee de l'air
T(1) = T_ext+273.15; % [K]
p(1) = 1.01325; % [bar] - pression atmospherique
%p(1) = 1;

Cpa_1 = Cpa(T(1));
h(1) = Cpa_1*(T(1)-273.15);
s(1) = Cpa_1*log(T(1)/273.15);
e(1) = (h(1)-h0) - T0*(s(1)-s0);

%% Etat 2 : Compression polytropique
p(2) = r*p(1);

f = @(T2) T(1).*r.^(Ra ./ ( eta_PiC.*Cpa(T2,T(1)) )) - T2;
T(2) = fsolve(f,400,optimoptions('fsolve','display','iter'));
h(2) = Cpa(T(2))*(T(2)-273.15);
s(2) = Cpa(T(2),T0)*log(T(2)/T0) - Ra*log(p(2)/p0);
e(2) = (h(2)-h0) - T0*(s(2)-s0);
WmC = h(2)-h(1);

%% Etat 3 : Combustion (quasi isobare)
p(3) = k_cc*p(2);
T(3) = T3+273.15; % [K]

% On determine le fuel utilise - type CHyOx
%Generelement, en cas stationnaire, methane (CH4) utilise -- p.111
LHV = 50009; % Pour methane [kJ/kg] trouve dans la litterature
y = 4;
x = 0;
ma1 = ((2*32 + 3.76*2*28)*(1+(y-2*x)/4)) / (12 + y + 32*x); % LMECA1855 p.133
exp = r^(Ra / (eta_PiC*Cpa(T(2),T(1))));

% On tente de trouver la valeur de lambda
% On egalise l'equation 3.9 (p115) a l'equation 3.26 (p119)
f = @(L) Cpa(T(2),T(1)).*T(1).*((1+1./(L.*ma1)) .* T(3)./T(1) .* ...
    Cpg(L,T(3),T(2))./Cpa(T(2),T(1)) - exp) - LHV./(L.*ma1);
lambda = fsolve(f,2,optimoptions('fsolve','display','iter'));

% On calcule le reste de l'etat
[~, Rg] = Cpg(lambda,T(3),T0)
h(3) = Cpg(lambda,T(3))*(T(3)-273.15);
s(3) = Cpg(lambda,T(3))*log(T(3)/T0) - Rg*log(p(3)/p0);
e(3) = (h(3)-h0) - T0*(s(3)-s0);


%% Etat 4 : Detente polytropique
p(4) = p(3)/(k_cc*r);

f = @(T4) T(3)*(p(4)/p(3))^(eta_PiT*Rg/Cpg(lambda,T(3),T4)) - T4; % p.119
T(4) = fsolve(f,T(1),optimoptions('fsolve','display','iter'));
T-273.15
p
lambda
h(4) = Cpg(lambda,T(4))*(T(4)-273.15);
s(4) = Cpg(lambda,T(4),T0)*log(T(4)/T0) - Rg*log(p(4)/p0);
e(4) = (h(4)-h0) - T0*(s(4)-s0);
WmT = h(3)-h(4);


%% Calcul des debits
%Systeme lineaire 3x3 a partir des eq. de la page 115:
% X = [m_a m_c m_g]'
MAT = [1 -lambda*ma1 0  ;  1 1 -1  ;  h(2) LHV -h(3)];
VEC = [0 0 0]';
X = MAT\VEC
m_a = X(1);
m_c = X(2);
m_g = X(3);
end

% Retourne la chaleur specifique moyenne de l'air, entre T1 et T2, sur base
% de janaf. Si nargout == 2, retourne la constante R* de l'air [kJ/kg/K]
% Si nargin == 1, retourne la chaleur specifique a la T° demandee
function [Cpa_m, Ra] = Cpa(T2,T1)
    % Masse molaire des composantes [kg/mol]
    MmO2 = 2*16e-3;
    MmN2 = 2*14.01e-3;

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
    Cpa_m = sum(Cpa_vec)/length(Tvec);
    
    if nargout > 1
        R = 8.314*1e-3; % [kJ/mol/K]
        Ra = R/Ma;
    end
end

% Retourne la chaleur specifique moyenne des fumees, entre T1 et T2, sur base
% de janaf. Si nargout == 2, retourne la constante R* des fumees [kJ/kg/K]
% Si nargin == 1, retourne la chaleur specifique a la T° demandee
function [Cpg_m, Rg] = Cpg(L,T2,T1)
    y = 4;
    x = 0;

    % kg/mol
    MmO2 = 2*16e-3;
    MmCO2 = (2*16+12.01)*1e-3;
    MmN2 = 2*14.01e-3;
    MmH2O = (2*1.01+2*16)*1e-3;
    
    % LMECA2160 : p.31 eq 3.33 avec k=0
    % fraction molaire [-]
    A = (y-2*x)/4;
    %denom = y/2 + (L-1)*(1+A)+3.76*L*(1+A);
    denom = 1;
    O2 = (L-1)*(1+A)/denom;
    CO2 = 1/denom;
    N2 = 3.76*L*(1+A)/denom;
    H2O = (y/2)/denom;
    
    M_g = MmO2*O2+MmCO2*CO2+MmN2*N2+MmH2O*H2O; % kg/mol * mol
    % kg/mol * mol / kg -> fraction massique
    MO2 = MmO2*O2/M_g;
    MCO2 = MmO2*CO2/M_g;
    MN2 = MmO2*N2/M_g;
    MH2O = MmO2*H2O/M_g;
    if nargin ==3
        Tvec = linspace(T1,T2,100);
    elseif nargin ==2
        Tvec = T2;
    end
    CpO2 = janaf('c','O2',Tvec);
    CpCO2 = janaf('c','CO2',Tvec);
    CpN2 = janaf('c','N2',Tvec);
    CpH2O = janaf('c','H2O',Tvec);
    Cpg_vec = MO2*CpO2+MCO2*CpCO2+MN2*CpN2+MH2O*CpH2O; % Vecteur selon T
    Cpg_m = sum(Cpg_vec)/length(Tvec);
    
    if nargout > 1
        R = 8.314*1e-3; % [kJ/mol/K]
        Rg = R/M_g * (y/2 + (L-1)*(1+A)+3.76*L*(1+A)); % [kJ/kg/K]
    end
end
