function [STATES,tbleed,pbleed,hbleed,sbleed,xbleed,ebleed,DAT9,Xflow,bache,DHbleed,DHRES] = Soutirages(STATES,options,eta_SiT,REHEAT)
% Calcule les soutirages et leurs etats intermediaires, ainsi que les
% fractions massiques de vapeur extraites

% On considere quelques hypotheses :
% 1) Il y a toujours un soutirage en sortie de HP (si nsout > 0, reheat > 0).
% 2) Pour les autres soutirages, ils sont repartis de manière
% "equidistants" au niveau enthalpique.
% 3) La détente dans les turbines est isentropique.
% 4) on peut considerer que p_6is = p_6i en sortie des bleeders.
% 5) La sortie d'un échangeur est à l'état de liquide saturé (avant la vanne).

%% Unfold des donnees
nsout = options.nsout;
Tdrum = options.Tdrum;
reheat = options.reheat;
drumFlag = options.drumFlag;
TpinchEx = options.TpinchEx;
TpinchSub = options.TpinchSub;
T_0 = options.T_0;
eta_SiC = options.eta_SiC;

t = STATES(1,:);
p = STATES(2,:);
h = STATES(3,:);
s = STATES(4,:);
x = STATES(5,:);
e = STATES(6,:);

if reheat == 2
    t4_2 = REHEAT(1);
    p4_2 = REHEAT(2);
    h4_2 = REHEAT(3);
    s4_2 = REHEAT(4);
    x4_2 = REHEAT(5);
    e4_2 = REHEAT(6);
end


if nsout == 0 % Pas de soutirage - court-circuitage comme cycle Rankin-Hirn
    p(1) = p(7);
    t(1) = t(7);
    x(1) = x(7);
    s(1) = s(7);
    h(1) = h(7);
    e(1) = e(7);
    tbleed = [];    pbleed = [];    hbleed = [];
    sbleed = [];    xbleed = [];    ebleed = [];
    DAT9 = [];      DHbleed = [];   DHRES = [];
    bache = 0;
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
        h(9) = XSteam('hL_p',p(9));
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
    DHres = DHres'; % Residus d'enthalpie "post" soutirage
    h9 = zeros(n+1,1); % Enthalpie du flux principal
    t9 = tbleed(3,:) - TpinchEx;
    if drumFlag == 1
        t9(b) = t(9);
    end
    t9 = [t9 0];
    DHbleed = hbleed(1,:) - hbleed(3,:);
    for i=1:nsout            
        if i <= b % apres bache
            h9(i) = XSteam('h_pT',p(10),t9(i));
        else % avant bache
            h9(i) = XSteam('h_pT',p(8),t9(i));
        end
    end
    h90 = h(8); % Valeur de depart pour les iterations
    while abs(h90-h9(n+1)) > 1e-4
        h9(n+1) = h90;
        for i=1:n
            if i == b
                DHliq(i) = h9(i)-h9(i+1) - (h(10)-h(9));
            elseif i == n
                DHliq(i) = h9(i)-h(8);
            else
                DHliq(i) = h9(i)-h9(i+1);
            end
        end

        %Creation de la matrice necessaire au systeme lineaire
        DHRES = DHres*ones(1,n);
        DHRES = tril(DHRES,-1);
        DHRES(n,:) = DHRES(n,:) + (hbleed(3,nsout)-hbleed(4,nsout))*ones(1,n);
        DHLIQ = DHliq*ones(1,n);
        if drumFlag == 1
            DHRES(b+1:n,1:b) = 0;
            DHRES(:,b) = 0;
            DHLIQ(b+1:n,1:b) = 0;
            DHLIQ(b,1:b) = 0;
        end
        DHBLEED = diag(DHbleed); 
        A = DHBLEED + DHRES - DHLIQ;
        B = DHliq;
        Xflow = A\B;
        
        h90 = sum(Xflow(b+1:n))*(hbleed(3,n)-hbleed(4,n))/(1+sum(Xflow(b+1:n))) + h(8);
    end
    t9(n+1) = XSteam('T_ph',p(8),h90);    
    % Reste de l'etat 9    
    t9 = t9';
    if drumFlag == 1
        p9 = [p(10)*ones(1,b-1) p(8)*ones(1,nsout-b+2)]';
    else
        p9 = [p(8)*ones(1,nsout-b+1)]';
    end
    x9 = arrayfun( @(p,h) XSteam('x_ph',p,h),p9,h9);
    s9 = arrayfun( @(p,h) XSteam('s_ph',p,h),p9,h9);
    e9 = arrayfun( @(h,s) Exergie(h,s),h9,s9);

elseif nsout ==1
    Xflow = (h(2)-h(7))/(hbleed(1,1)-h(2));
else 
    %Rankin-Hirn
    Xflow = 0;
end


%% Renvoi des donnees
if nsout > 0
    DAT9 = [t9 p9 h9 s9 x9 e9];
end
STATES(1,:) = t;
STATES(2,:) = p;
STATES(3,:) = h;
STATES(4,:) = s;
STATES(5,:) = x;
STATES(6,:) = e;


%% Fonctions annexes
% Calcule l'exergie
function e = Exergie(h , s)
    h_0= XSteam('hL_T',T_0);
    s_0= XSteam('sL_T',T_0);
    e = (h-h_0) - (273.15+T_0)*(s-s_0); 
end

% Trouve la pression en debut de soutirage
function f = pressionBleeders(ptest)
    [~,hbl,~,~,~,~] = detenteTurb(t(5),p(5),h(5),s(5),ptest,eta_SiT(2));
    f = hbl - hbleed(1,i);
end
end

