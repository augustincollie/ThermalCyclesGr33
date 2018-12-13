% VARIABLES
% Tdb (dry bulb temperature) and Tdp(dew point temperature) in C
% w (humidity ratio) in kg/kg of dry air
% phi (relative humidity) in %
% h (enthalpy) in J/kg of dry air
% v (specific volume) in m3/kg of dry air
% Twb (wet bulb temperature) in C
% P (atmospheric pressure) in kPa

% The following cases are present:
% Tdb, w; Tdb, phi; Tdb, h; w, phi; w, h; phi, h; Tdb, Twb; w, Twb; phi, Twb;  

% Following ASHRAE 2013 Fundamentals SI Psychrometrics chapter equations are used:
% Eq6:Pws=f(Tdb); Eq22: w=f(Tdb, phi, p); Eq24: phi=f(Tdb, w, p); Eq28:v=f(Tdb, w and p); Eq32:h=f(Tdb, w and p); 
% Eq35:Twb=f(Tdb,w); Eq39:Tdp=f(Tdb, p);
%
% EXAMPLE, I know the relative humidity and the dry bulb temperature : 
% [Tdb, w, phi, h, Tdp, v, Twb] = Psychometrics('Tdb',300,'phi',0.8);
function [Tdb, w, phi, h, Tdp, v, Twb] = Psychrometrics(varargin)

if length(varargin)<4
    display('Need four inputs:''prop1'',value1,''prop2'',value2''');
    Tdb=[];w=[];phi=[];h=[];Tdp=[];v=[];Twb=[];    
    return
elseif length(varargin)>4 && length(varargin)<6
    display('Need six inputs:''prop1'',value1,''prop2'',value2'',,''Pamb'',value in kPa''');
    Tdb=[];w=[];phi=[];h=[];Tdp=[];v=[];Twb=[];    
    return
elseif length(varargin)==4
Tdb_in=[];w_in=[];phi_in=[];h_in=[];Twb_in=[];    
prop(1) = {lower(char(varargin(1)))};
prop(2) = {lower(char(varargin(3)))};
propVal(1) = cell2mat(varargin(2));
propVal(2) = cell2mat(varargin(4));
P = 101.325; %kPa, atmospheric pressure value
elseif length(varargin)==6
Tdb_in=[];w_in=[];phi_in=[];h_in=[];Twb_in=[];    
prop(1) = {lower(char(varargin(1)))};
prop(2) = {lower(char(varargin(3)))};
propVal(1) = cell2mat(varargin(2));
propVal(2) = cell2mat(varargin(4));
P = cell2mat(varargin(6));
end

for i=1:2    
switch prop{i}
    case 'tdb'
        Tdb_in=propVal(i);
    case 'w'
        w_in=propVal(i);
    case 'phi'
        phi_in=propVal(i);
    case 'h'
        h_in=propVal(i);
    case 'twb'
        Twb_in=propVal(i);
        
end
end
if (~isempty(Twb_in) && ~isempty(h_in))
    display('function not available');
    Tdb=[];w=[];phi=[];h=[];Tdp=[];v=[];Twb=[];    
    return
end
c_air = 1006; %J/kg, value from ASHRAE 2013 Fundamentals eq. 32
hlg = 2501000; %J/kg, value from ASHRAE 2013 Fundamentals eq. 32
cw  = 1860; %J/kg, value from ASHRAE 2013 Fundamentals eq. 32

%++++++++++++++++++++++++++++++++++++++
if (~isempty(Tdb_in) && ~isempty(w_in))
    Tdb=Tdb_in;w=w_in;
    
    % phi calculation from Tdb and w
    Pw=w*P/(0.621945+w); %partial pressure of water wapor
    Pws=Saturation_pressure(Tdb);
    phi=Pw/Pws*100;
    
    % h calculation from Tdb and w
    h=c_air*Tdb+w*(hlg+cw*Tdb); %ASHRAE 2013 fundamentals eq. 32
    
    % v calculation from Tdb and w
    v=0.287042*(Tdb+273.15)*(1+1.607858*w)/P; %ASHRAE 2013 fundamentals eq. 28
end
%++++++++++++++++++++++++++++++++++++++
if (~isempty(Tdb_in) && ~isempty(phi_in))
    Tdb=Tdb_in;phi=phi_in;
    
    % w calculation from Tdb and phi
    Pws=Saturation_pressure(Tdb);
    Pw=phi/100*Pws;
    w=0.621945*Pw/(P-Pw);
    
    % h calculation from Tdb and w
    h=c_air*Tdb+w*(hlg+cw*Tdb);
    
    % v calculation from Tdb and w
    v=0.287042*(Tdb+273.15)*(1+1.607858*w)/P;
end
%++++++++++++++++++++++++++++++++++++++
if (~isempty(Tdb_in) && ~isempty(h_in))
    Tdb=Tdb_in;h=h_in;
    
    % w calculation from Tdb and h
    w=(h - c_air*Tdb)/(hlg+cw*Tdb);
    
    % phi calculation from Tdb and w
    Pw=w*P/(0.621945+w); %partial pressure of water wapor
    Pws=Saturation_pressure(Tdb);
    phi=Pw/Pws*100;
    
    % v calculation from Tdb and w
    v=0.287042*(Tdb+273.15)*(1+1.607858*w)/P;
end
%++++++++++++++++++++++++++++++++++++++
if (~isempty(w_in) && ~isempty(h_in))
    w=w_in;h=h_in;
    
    % Tdb calculation from w and h
    Tdb=(h - w*hlg)/(c_air+w*cw);
    
    % phi calculation from Tdb and w
    Pw=w*P/(0.621945+w); %partial pressure of water wapor 
    Pws=Saturation_pressure(Tdb);
    phi=Pw/Pws*100;
    
    % v calculation from Tdb and w
    v=0.287042*(Tdb+273.15)*(1+1.607858*w)/P;
end
%++++++++++++++++++++++++++++++++++++++
if (~isempty(w_in) && ~isempty(phi_in))
    w=w_in;phi=phi_in;
    
    % Tdb calculation from phi and w
    Pw=w*P/(0.621945+w); %partial pressure of water wapor
    Pws=Pw/phi*100;
    options=optimset('LargeScale','off','Display','off');
    [y,val,exitflag]=fsolve(@Iteration_function_1, 20,options);Tdb =y(1);
    if exitflag<1
        disp('Iteration error')
    end
    
    % h calculation from Tdb and w
    h=c_air*Tdb+w*(hlg+cw*Tdb);
    
    % v calculation from Tdb and w
    v=0.287042*(Tdb+273.15)*(1+1.607858*w)/P;
end
%++++++++++++++++++++++++++++++++++++++
if (~isempty(phi_in) && ~isempty(h_in))
    phi=phi_in;h=h_in;
    
    % Tdb calculation from phi and h    
    options=optimset('LargeScale','off','Display','off');
    [y,val,exitflag]=fsolve(@Iteration_function_2, 20,options);Tdb =y(1);
    if exitflag<1
        disp('Iteration error')
    end
    
    % w calculation from Tdb and phi
    Pws=Saturation_pressure(Tdb);
    Pw=phi/100*Pws;
    w=0.621945*Pw/(P-Pw);
    
    % h calculation from Tdb and w
    h=c_air*Tdb+w*(hlg+cw*Tdb);
    
    % v calculation from Tdb and w
    v=0.287042*(Tdb+273.15)*(1+1.607858*w)/P;
end
%++++++++++++++++++++++++++++++++++++++
if (~isempty(Tdb_in) && ~isempty(Twb_in))
    Tdb=Tdb_in;Twb=Twb_in;
    
    % w calculation from Tdb and Twb
    Pws=Saturation_pressure(Tdb);
    Pwsasterik=Saturation_pressure(Twb);
    ws=0.621945*Pwsasterik/(P-Pwsasterik);
    w= ((hlg-2.326e3*Twb)*ws-c_air*(Tdb-Twb))/(hlg+cw*Tdb-4.186e3*Twb);
    
    % phi calculation from Tdb and w
    Pw=w*P/(0.621945+w); %partial pressure of water wapor
    phi=Pw*100/Pws;
    
    % h calculation from Tdb and w
    h=c_air*Tdb_in+w*(hlg+cw*Tdb_in);
    
    % v calculation from Tdb and w
    v=0.287042*(Tdb_in+273.15)*(1+1.607858*w)/P;
end
%++++++++++++++++++++++++++++++++++++++
if (~isempty(w_in) && ~isempty(Twb_in))
    w=w_in;Twb=Twb_in;
    
    % Tdb calculation from Twb and w
    Pwsasterik=Saturation_pressure(Twb);
    ws=0.621945*Pwsasterik/(P-Pwsasterik);
    options=optimset('LargeScale','off','Display','off');
    [y,val,exitflag]=fsolve(@Iteration_function_4, Twb,options);Tdb =y(1);
    if exitflag<1
        disp('Iteration error')
    end

    % phi calculation from Tdb and w
    Pws=Saturation_pressure(Tdb);
    Pw=w*P/(0.621945+w); %partial pressure of water wapor
    phi=Pw*100/Pws;
    
    % h calculation from Tdb and w
    h=c_air*Tdb_in+w*(hlg+cw*Tdb_in);
    
    % v calculation from Tdb and w
    v=0.287042*(Tdb_in+273.15)*(1+1.607858*w)/P;    
end
%++++++++++++++++++++++++++++++++++++++
if (~isempty(phi_in) && ~isempty(Twb_in))
    phi=phi_in;Twb=Twb_in;
    
    % Tdb calculation from phi and Twb
    Pwsasterik=Saturation_pressure(Twb);
    ws=0.621945*Pwsasterik/(P-Pwsasterik);
    options=optimset('LargeScale','off','Display','off');  
    [y,val,exitflag]=fsolve(@Iteration_function_5, Twb,options);Tdb =y(1);
    if exitflag<1
        disp('Iteration error')
    end

    % w calculation from Tdb and phi
    Pws=Saturation_pressure(Tdb);
    Pw=phi/100*Pws;
    w=0.621945*Pw/(P-Pw);
    
    % h calculation from Tdb and w
    h=c_air*Tdb+w*(hlg+cw*Tdb);
    
    % v calculation from Tdb and w
    v=0.287042*(Tdb+273.15)*(1+1.607858*w)/P;
end

% dew point calculation from w
% pw=(P*w)/(0.621945+w); % water vapor partial pressure in kPa
alpha=log(Pw);
Tdp=6.54 + 14.526*alpha+0.7389*(alpha^2)+0.09486*(alpha^3)+0.4569*(Pw^0.1984); % valid for Tdp between 0 C and 93 C

%++++++++++++++++++++++++++++++++++++++
if nargout>6 && isempty(Twb_in)
% Note: this Twb calc. equations are good for patm=101325 Pa only. 
if abs(Tdb - Tdp) < .001, Twb=Tdb;return;end
options=optimset('LargeScale','off','Display','off');
[y,val,exitflag]=fsolve(@Iteration_function_3, Tdb,options);Twb=y(1);
if Twb > Tdb,Twb=Tdb;end
if Twb < Tdp,Twb=Tdp;end   
end
% if phi>100   
%     Tdb = NaN;
%     w   = NaN;
%     phi = NaN;
%     h   =  NaN;
%     Tdp = NaN;
%     v  =NaN; 
%     
%     disp('ERROR: Point is outside the chart')
% end


    function [Pws] = Saturation_pressure(Tdb) %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6
        T=Tdb+273.15;
        Pws=exp(-(5.8002206e3)/T+1.3914993+-(4.8640239e-2)*T+(4.1764768e-5)*(T^2)-(1.4452093e-8)*(T^3)+6.5459673*log(T)); %in Pa valid for 0 to 200C
        Pws=Pws/1000; % in kPa
    end


    function result = Iteration_function_1(y) %calc Tdb from phi and w
        Tdb_as=y(1);
        Pws=Saturation_pressure(Tdb_as);
        phi_as=Pw/Pws*100; %ASHRAE 2013 fundamentals eq. 24
        % equation to satisfy
        result=phi_as-phi;
    end


    function result = Iteration_function_2(y) %calc Tdb from phi and h
        Tdb_as=y(1);
        % w calculation from Tdb and phi
        Pws=Saturation_pressure(Tdb_as);
        Pw=phi/100*Pws;
        w_as=0.621945*Pw/(P-Pw); %ASHRAE 2013 fundamentals eq. 22
        % h calculation from Tdb and w
        h_as=c_air*Tdb_as+w_as*(hlg+cw*Tdb_as);
        % equation to satisfy
        result=h_as-h;
    end
 
    function result = Iteration_function_3(y) %calc Twb from Tdb and w using ASHRAE 2013 fundamentals eq. 35
        Twb_as=y(1);
        Pws_as=Saturation_pressure(Twb_as);
        ws=0.621945*Pws_as/(P-Pws_as);
        w_as= ((hlg-2.326e3*Twb_as)*ws-c_air*(Tdb-Twb_as))/(hlg+cw*Tdb-4.186e3*Twb_as);
        result=(w-w_as)*1000;
    end
        
    function result = Iteration_function_4(y) %calc Tdb from Twb and w Tdp using ASHRAE 2013 fundamentals eq. 35
        Tdb_as=y(1);
        w_as= ((hlg-2.326e3*Twb)*ws-c_air*(Tdb_as-Twb))/(hlg+cw*Tdb_as-4.186e3*Twb);
        result=(w-w_as)*1000;
    end

    function result = Iteration_function_5(y) %calc Tdb from Twb and phi Tdp using ASHRAE 2013 fundamentals eq. 35
        Tdb_as=y(1);
        w_as= ((hlg-2.326e3*Twb)*ws-c_air*(Tdb_as-Twb))/(hlg+cw*Tdb_as-4.186e3*Twb);
        Pw_as=w_as*P/(0.621945+w_as); %partial pressure of water wapor
        Pws_as=Saturation_pressure(Tdb_as);
        phi_as=Pw_as*100/Pws_as;
        result=phi-phi_as;
    end

    function result = Iteration_function_6(y) %calc Pw from Tdb and Tdp using ASHRAE 2013 fundamentals eq. 39
        Pw_as=y(1);
        Tdp_as= 6.54+14.526*log(Pw_as)+0.7389*(log(Pw_as))^2+0.09486*(log(Pw_as))^3+0.4569*(Pw_as^0.1984); % valid for Tdp between 0 C and 93 C
        result=Tdp-Tdp_as;
    end
end