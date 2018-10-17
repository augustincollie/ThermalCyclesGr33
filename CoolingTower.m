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