function [DAT_WATER DAT_AIR MASSFLOW] = CoolingTower(P_w,options)
% COOLINGTOWER is a cooling tower 0D modelisation
% COOLINGTOWER(P_w,options) compute the thermodynamics states for a Cooling
% tower based on several inputs (given in OPTION) and based on a given 
% water power to dissipate P_w.
% It returns the main results. 
%
% INPUTS :
% P_W = Heat power output at the condenser [kW]
% OPTIONS is a structure containing :
%   -options.Tcond  [°C]: Temperature in the condenser
%   -options.Tpinch [°C]: Minimum tempearture pinch between Tw_out and the
%                         condenser temperature.
%   -options.Tw_out [°C]: Cooling water temperature at the condenser outlet
%   -options.Tw_in  [°C]: Cooling water temperature at the condenser inlet
%   -options.Triver [°C]: River temperature 
%   -options.Ta_in  [°C]: Atmospheric air temperature 
%   -options.Ta_out [°C]: Air outlet temperature of cooling tower 
%   -options.Phi_atm [-]: Relative humidity of atmospheric air
%   -options.Phi_out [-]: Maximum relative humidity of air at the cooling 
%                         tower outlet.
%
% OUTPUT :
% MassFlow [kg/s]: Vector containing the different massflow :
%   -massflow(1) : water massflow at the condenser
%   -massflow(2) : additionnal water massflow = water flow evaporated
%   -massflow(3) : air massflow at the cooling tower   
%
%  dat_water = [T_e1       , T_e2       , T_e3       , T_e4;  %[°C]
%               h_e1       , h_e2       , h_e3       , h_e4;  %[kJ/kg]
%               m_e1       , m_e2       , m_e3       , m_e4]; %[kg/s]
% 
%  dat_air   = [Ta_in       , Ta_out  ;  %[°C]
%               ha_in       , ha_out  ;  %[kJ/kg]
%               xa_in       , xa_out  ;  %[kg_water/kg_dry_air]
%               Phia_in     , Phia_out]; %[-] relative humidity
%  
%
% ADDITIONNAL INFORMATIONS
% Water points : 
%       1 : water outlet of cooling tower
%       2 : water just before condenser
%       3 : water just after  condenser
%       4 : water from the river (coming between 1 & 2)
%
% Air points :
%       a_in : air at the cooling tower inlet
%       a_out : air at the cooling tower outlet
%

%% YOUR WORK
if nargin<2
    options=struct();
    if nargin<1
        P_w=200e3;%200MW_heat
    end
end

if isfield(options,'Tpinch')
    Tpinch = options.Tpinch;
else
    Tpinch = 4; %[K]
end
end

%%Importation des données
for d=1 %Only to allow collapsing of the code...
    if isfield(options,'Tcond')
        Tcond = options.Tcond;
    else
        Tcond = 46.8;  % [C]
    end
   
    if isfield(options,'Tw_out')
        if (options.Tw_out > (Tcond-Tpinch))
        Tw_out = options.Tw_out;
        end
    else
        Tw_out = Tcond-Tpinch; %42.8;  % [C]
    end
    if isfield(options,'Tw_in')
        Tw_in = options.Tw_in;
    else
        Tw_in = 20;  % [C]
    end
    if isfield(options,'Triver')
        Triver = options.Triver;
    else
        Triver = 15;  % [C]
    end
    if isfield(options,'Ta_in')
        Ta_in = options.Ta_in;
    else
        Ta_in = 15;  % [C]
    end
    if isfield(options,'Ta_out')
        Ta_out = options.Ta_out;
    else
        Ta_out = 25;  % [-]
    end
    if isfield(options,'Phi_atm')
        Phi_atm= options.Phi_atm;
    else
        Phi_atm = 0.7;  % [-]
    end
    if isfield(options,'Phi_out')
        Phi_out = options.Phi_out;
    else
        Phi_out = 1;  % [-]
    end
end
    
%Allocation Outputs
    MASSFLOW = zeros(1,3); % [kg/s]
    DAT_WATER = zeros(3,4);
    DAT_AIR = zeros(4,2);

%Donnees Utiles
    Cpa = 1.009; % chaleur specifique air [kJ/kgK]
    Cpv = 1.854; % chaleur specifique vapeur d'eau [kJ/kgK]
    Cpl = 4.1868; % chaleur specifique eau liquide [kJ/kgK]
    hLV = 2.5016; % enthalpie d'evaporation de l'eau [kJ/kg]
    
    
%% LAST TRY ! 
    % Given outputs
    T_e2 = Tw_in; %[C]
    T_e3 = Tw_out;
    T_e4 = Triver;
    
    h_e2 = Cpl*T_e2; %[kJ/kg]
    h_e3 = Cpl*T_e3;
    h_e4 = Cpl*T_e4;
    
    m_e3 = P_w/(h_e3-h_e2); %[kg/s]
    m_e2 = m_e3;
    
    Ta_1 = Ta_in; %[C]
    Ta_2 = Ta_out;
    
    Phia_in = Phi_atm; %[-]
    Phia_out = Phi_out;
    
    %States of cooling Air - Psychrometrics
    [~, xa_in, ~, ha_in, ~, ~, ~] = Psychrometrics('Tdb',Ta_1,'phi',Phia_in);
    [~, xa_out, ~, ha_out, ~, ~, ~] = Psychrometrics('Tdb',Ta_2,'phi',Phia_out);
    ha_in = ha_in*1e-03;
    ha_out = ha_out*1e-03;
    
    %System 3 Eq - 3 Unkn.
    
    function Inc = FunS(X)
       Inc(1) = X(3)-(X(1)*(xa_out-xa_in)/(1+xa_in)) %Mass Conservation
       
       Term21 = m_e3*(h_e3-Cpl*X(2));
       Term22 = X(1)*((ha_out-ha_in)-(Cpl*X(2)*(xa_out-xa_in)))/(1+xa_in);
       Inc(2) = Term21 - Term22 %Enthalpy conservation at cooling tower
       
       Inc(3) = X(2)-((m_e3*h_e2-X(3)*h_e4)/((m_e3*X(3))*Cpl)) %Enthalpy conservation in water between 3 and 1
    end

    
    Inc = ones(1,3); %Vectors of 3 unknowns of the system; [mDotAir, T_e1, m_e4]
    Init = 1e04*Inc; %Initial iteration values of Inc
    
%     options = optimoptions('fsolve','OptimalityTolerance',1e-04);
    Ans = fsolve(@FunS, Init); %,options
    
    mDotAir = Ans(1);
    T_e1 = Ans(2);
    m_e4 = Ans(3);
    
    m_e1 = m_e3-m_e4;
    h_e1 = Cpl*T_e1;
    
    %Preparing Outputs
    
    MASSFLOW = [m_e3, m_e4, mDotAir] %[kg/s]
          
    DAT_WATER = [ T_e1       , T_e2       , T_e3       , T_e4;      %[°C]
                  h_e1       , h_e2       , h_e3       , h_e4;      %[kJ/kg]
                  m_e1       , m_e2       , m_e3       , m_e4 ]    %[kg/s]
    DAT_AIR   = [ Ta_in       , Ta_out ;     %[°C]
                  ha_in       , ha_out ;     %[kJ/kg]
                  xa_in       , xa_out ;     %[kg_water/kg_dry_air]
                  Phia_in     , Phia_out ]  %[-] relative humidity

    
    
%% Anciennes méthodes.. A SUPPRIMER    
%Evaluation du flux massique d'eau
    dhCond = Cpl*(Tw_out - Tw_in); % [kJ/kg]
    MASSFLOW(1) = P_w/dhCond % [kg/s]
    
%     plot(Cp_e,'ok'); %Vérifier la linéarité de Cp_eau
    
%Calcul des etats d'entree et sortie de l'air atmospherique
%     p_sat1 = (XSteam('psat_T',Ta_in)) * 1e06; % en Pa
%     p_sat2 = (XSteam('psat_T',Ta_out))* 1e06; % en Pa

%     x1 = (1/0622)*((Phi_atm*p_sat1)/(1e05-p_sat1)); 
%     x2 = (1/0622)*((Phi_out*p_sat2)/(1e05-p_sat2)); 
    [Tdb1, x1, Phi1, h1, Tdp1, v1, Twb1] = Psychrometrics('Tdb',(Ta_in),'phi',Phi_atm);
    [Tdb2, x2, Phi2, h2, Tdp2, v2, Twb2] = Psychrometrics('Tdb',(Ta_out),'phi',Phi_out);
    
    h1 = h1*1e-03; %[kJ/kg_dryAir]
    h2 = h2*1e-03; %[kJ/kg_dryAir]
    
DAT_AIR = [Tdb1, Tdb2; NaN, NaN; x1, x2; Phi1, Phi2]

%Calcul des etats de l'eau du circuit de refroidissement
    

    DAT_WATER(3,2) = MASSFLOW(1);
    DAT_WATER(3,3) = MASSFLOW(1);
    DAT_WATER(1,:) = [NaN Tw_in Tw_out Triver];
    DAT_WATER(2,:) = Cpl*DAT_WATER(1,:);
    
%     mDotEvap = sqrt(((h2-h1)/DAT_WATER(2,3))*(x2-x1)); %calcul de la masse d'eau evaporee
%     MASSFLOW(2) = mDotEvap;

    he_evap = (h2-h1)/(x2-x1);
    
    mDotDryAir = (MASSFLOW(1)*he_evap)/(h2-h1);
    MASSFLOW(3) = (1+x1)*mDotDryAir;
    MASSFLOW(2) = MASSFLOW(3)*(x2-x1)
    
    AirFlowRatio = mDotDryAir/MASSFLOW(3);
    DAT_AIR(2,:) = AirFlowRatio*[h1, h2]
    
    DAT_WATER(3,4) = MASSFLOW(2);
    DAT_WATER(3,1) = MASSFLOW(1)-MASSFLOW(2);
    dhAir = h2-h1;
    DAT_WATER(1,1) = (dhAir-DAT_WATER(3,3)*DAT_WATER(2,3))/(-DAT_WATER(3,1)*Cpl); %Evaluation de la temperature au point 1
    DAT_WATER(2,1) = DAT_WATER(1,1)*Cpl;
    

    
    MASSFLOW
    DAT_WATER
    DAT_AIR
    
end
