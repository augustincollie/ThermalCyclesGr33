function [ ] = AnalyseST( )
% Fait des analyses de rendement en fonction de differents parametres

options = struct('nsout',8,'reheat',1,'T_max',525,'T_cond_out',30,'p3_hp',200, ...
            'drumFlag',1,'eta_mec',0.98,'comb',struct('Tmax',1200,'x',0,'y',4) ...
            ,'T_exhaust',120,'x4',0.89,'T_0',15,'TpinchSub',4,'TpinchEx',15, ...
            'TpinchCond',15,'Tdrum',120,'eta_SiC',0.9,'eta_SiT',[0.9 0.9]);
P_e = 225e3;  

% cyclen//toten//cyclex//totex//gen//gex//combex//chemex//transex
% ETA = zeros(2,9);
% for i=1:2
%     options.reheat = i;
%     [ETA(i,:),~,~,~,~,~,~] = ST(P_e,options,0);    
% end
% ETA

ETA = zeros(11,9);
for i=0:10
    options.p3_hp = 200+3*i;
    [ETA(i+1,:),~,~,~,~,~,~] = ST(P_e,options,0);
end
figure;
k=200:3:230;
plot(k,ETA(:,3));
hold on;
plot(k,ETA(:,4));
plot(k,ETA(:,6));
plot(k,ETA(:,7));
plot(k,ETA(:,8));
plot(k,ETA(:,9));
hold off
legend('cyclex','totex','gex','combex','chemex','transex')