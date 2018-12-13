function [FIG] = DrawST(DAT, options, eta_SiT, tbleed, pbleed, hbleed, sbleed, DAT9, DATVEC, REHEAT)
% Trace les diagrammes h-s et t-s de la Steam Turbine

%% Extraction des donnees utiles
t = DATVEC(1,:);
p = DATVEC(2,:);
h = DATVEC(3,:);
s = DATVEC(4,:);
nsout = options.nsout;
reheat = options.reheat;
if nsout > 0
    t9 = DAT9(:,1);
    h9 = DAT9(:,2);
    s9 = DAT9(:,3);
end

if reheat == 2
    t3_2 = REHEAT(1,1);     t4_2 = REHEAT(1,2);
    p3_2 = REHEAT(2,1);     p4_2 = REHEAT(2,2);
    h3_2 = REHEAT(3,1);     h4_2 = REHEAT(3,2);
    s3_2 = REHEAT(4,1);     s4_2 = REHEAT(4,2);
end

%% Plots
%%%%%%%%%%%%%%%%%%% linearisation entre les etats
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

        % 8 - 1 : Parcours des echanges de chaleur
        ts89 = fliplr(t9');
        hs89 = fliplr(h9');
        s89 = fliplr(s9');
        
        % Soutirage de surchauffe a vapeur sature
        tsout = zeros(50,nsout);
        hsout = zeros(50,nsout);
        ssout = zeros(50,nsout);
        for k=1:nsout
            ssout(:,k) = linspace(sbleed(1,k),sbleed(2,k),50);
            tsout(:,k) = arrayfun( @(s) XSteam('T_ps',pbleed(1,k),s), ssout(:,k));
            hsout(:,k) = arrayfun( @(s) XSteam('h_ps',pbleed(1,k),s), ssout(:,k));
        end
        ssout = [ssout ; sbleed(3:4,:)];
        tsout = [tsout ; tbleed(3:4,:)];
        hsout = [hsout ; hbleed(3:4,:)];
    else
        s78 = [];
        s89 = [];
        ts78 = [];
        ts89 = [];
        hs78 = [];
        hs89 = [];
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
    plot([s12 s23 s34 s43p s3p4p s4p5 s45 s56 s67 s78 s89] , [ts12 ts23 ts34 ts43p ts3p4p ts4p5 ts45 ts56 ts67 ts78 ts89],'r');    
    if nsout > 0
        plot(ssout,tsout,'--')
    end
    plot(DAT(4,:),DAT(1,:),'*')
    legend('Courbe de saturation','Fluide principal','Soutirages','Location','NorthWest')
    title('Diagramme T-s de l''installation')
    xlabel('Entropie [kJ/kg/K]')
    ylabel('Température [°C]')
    
    FIG(2) = figure; % diagramme h-s
    plot([sliq svap],[hliq hvap],'-k'); % Cloche de saturation
    hold on
    plot([s12 s23 s34 s43p s3p4p s4p5 s45 s56 s67 s78 s89] , [hs12 hs23 hs34 hs43p hs3p4p hs4p5 hs45 hs56 hs67 hs78 hs89],'r');
    if nsout > 0
        plot(ssout,hsout,'--')
    end
    plot(DAT(4,:),DAT(3,:),'*')
    legend('Courbe de saturation','Fluide principal','Soutirages','Location','SouthEast')    
    title('Diagramme h-s de l''installation')
    xlabel('Entropie [kJ/kg/K]')
    ylabel('Enthalpie [kJ/kg]')
    
    
end

