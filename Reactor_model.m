function m = Reactor_model(t_final,dt_plot,P_stable,PF_retarded)
tic
%%

%%%%%%%%%%%%%%%%%%%%%%%
%Donnees et hypotheses%
%%%%%%%%%%%%%%%%%%%%%%%

NA = 6.02214076e23; %avogadro
V = 10; %volume du reacteur [m^3]

m_Utot = 1000; %masse d'uranium [kg]
m_U235 = m_Utot*0.03; %masse d'uranium 235
m_U238 = m_Utot*0.97; %masse d'uranium 238

N_U235 = m_U235/molarMass('U235'); %nombre de moles d'U235 [mol]
N_U238 = m_U238/molarMass('U238'); %nombre de moles d'U238 [mol]

E_thermique = 0.025; %energie neutron thermique [eV]
E_rapide = 1e6; %energie neutron rapide [eV]

v_thermique = 10; %vitesse neutron thermique [m/s]
v_rapide = 1000; %vitesse neutron rapide [m/s]

n_0 = 1e10; %nombre initial de neutrons (thermiques)

lambda_PF = log(2)/1; %lambda pour PF* = PF + n + 5 MeV : log(2)/1
lambda_rt = log(2)/(5*10^-4); %lambda transition rapide->thermique : log(2)/(5*10^-4)
lambda_perte_th = 12; %lambda neutron thermique pour fuites et barres de controle : 800
lambda_perte_rap = 65; %lambda neutron rapide pour fuites et barres de controle : 2000

%--------------------------------------------------------------------------
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculs des sections efficaces et des temps de demi-vie%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%U235
fis_U235_th = 1e-28*Section_efficace('U235','Fission',E_thermique,'DATABASE');
fis_U235_rap = 1e-28*Section_efficace('U235','Fission',E_rapide,'DATABASE');
    
%U238
fis_U238_th = 1e-28*Section_efficace('U238','Fission',E_thermique,'DATABASE');
fis_U238_rap = 1e-28*Section_efficace('U238','Fission',E_rapide,'DATABASE');
cap_U238_th = 1e-28*Section_efficace('U238','Capture',E_thermique,'DATABASE');
cap_U238_rap = 1e-28*Section_efficace('U238','Capture',E_rapide,'DATABASE');

%U239
fis_U239_th = 1e-28*Section_efficace('U239','Fission',E_thermique,'DATABASE');
fis_U239_rap = 1e-28*Section_efficace('U239','Fission',E_rapide,'DATABASE');
demi_U239 = Demi_vie('U239','BetaMinus');
    
%Np239
fis_Np239_th = 1e-28*Section_efficace('Np239','Fission',E_thermique,'DATABASE');
fis_Np239_rap = 1e-28*Section_efficace('Np239','Fission',E_rapide,'DATABASE');
demi_Np239 = Demi_vie('Np239','BetaMinus');
    
%Pu239
fis_Pu239_th = 1e-28*Section_efficace('Pu239','Fission',E_thermique,'DATABASE');
fis_Pu239_rap = 1e-28*Section_efficace('Pu239','Fission',E_rapide,'DATABASE');

%Xe135
cap_Xe135_th = 1e-28*Section_efficace('Xe135','Capture',E_thermique,'DATABASE');
cap_Xe135_rap = 1e-28*Section_efficace('Xe135','Capture',E_rapide,'DATABASE');
demi_Xe135 = Demi_vie('Xe135','BetaMinus');

%--------------------------------------------------------------------------
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialisation des vecteurs%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_final = 10; %[s]
dt_gen = 10^-4;
T = [0:dt_gen:t_final];
Y = zeros(length(T),8); %U235,U238,U239,Np239,Pu239,PF*,Xe135,PF
Y(1,:) = [N_U235 N_U238 0 0 0 0 0 0]; %quantites initiales [mol]
n_th = zeros(length(T),1); %neutrons thermiques 
n_th(1,1) = n_0; %nombre initial [#]
n_rap = zeros(length(T),1); %neutrons rapides
n_rap(1,1) = 0; %nombre initial [#]
phi_th = zeros(length(T),1); %flux de neutrons thermiques
phi_th(1,1) = n_th(1,1)*v_thermique/V; %flux initial [#/m^2.s]
phi_rap = zeros(length(T),1); %flux de neutrons rapides
phi_rap(1,1) = n_rap(1,1)*v_rapide/V; %flux initial [#/m^2.s]

% n = zeros(length(T),1); %tous les neutrons (thermique+rapide)
% n(1,1) = n_th(1,1) + n_rap(1,1);
% n_ret = zeros(length(T),1); %neutrons retardes
% n_ret(1,1) = 0;

%--------------------------------------------------------------------------
%%

%%%%%%%%%%%%%%%%%%%%%%
%Equations cinetiques%
%%%%%%%%%%%%%%%%%%%%%%

for i = 2:length(T)
    
    Y(i,1) = Y(i-1,1) + (- Y(i-1,1)*fis_U235_th*phi_th(i-1,1) - Y(i-1,1)*fis_U235_rap*phi_rap(i-1,1))*dt_gen;
    Y(i,2) = Y(i-1,2) + (- Y(i-1,2)*fis_U238_th*phi_th(i-1,1) - Y(i-1,2)*fis_U238_rap*phi_rap(i-1,1) - Y(i-1,2)*cap_U238_th*phi_th(i-1,1) - Y(i-1,2)*cap_U238_rap*phi_rap(i-1,1))*dt_gen;
    Y(i,3) = Y(i-1,3) + (Y(i-1,2)*cap_U238_th*phi_th(i-1,1) + Y(i-1,2)*cap_U238_rap*phi_rap(i-1,1) - Y(i-1,3)*fis_U239_th*phi_th(i-1,1) - Y(i-1,3)*fis_U239_rap*phi_rap(i-1,1) - Y(i-1,3)*log(2)/demi_U239)*dt_gen;
    Y(i,4) = Y(i-1,4) + (Y(i-1,3)*log(2)/demi_U239 - Y(i-1,4)*fis_Np239_th*phi_th(i-1,1) - Y(i-1,4)*fis_Np239_rap*phi_rap(i-1,1) - Y(i-1,4)*log(2)/demi_Np239)*dt_gen; 
    Y(i,5) = Y(i-1,5) + (Y(i-1,4)*log(2)/demi_Np239 - Y(i-1,5)*fis_Pu239_th*phi_th(i-1,1) - Y(i-1,5)*fis_Pu239_rap*phi_rap(i-1,1))*dt_gen;
    Y(i,6) = Y(i-1,6) + (Y(i-1,1)*fis_U235_th*phi_th(i-1,1) + Y(i-1,2)*fis_U238_th*phi_th(i-1,1) + Y(i-1,3)*fis_U239_th*phi_th(i-1,1) + Y(i-1,4)*fis_Np239_th*phi_th(i-1,1) + Y(i-1,5)*fis_Pu239_th*phi_th(i-1,1) + Y(i-1,1)*fis_U235_rap*phi_rap(i-1,1) + Y(i-1,2)*fis_U238_rap*phi_rap(i-1,1)+ Y(i-1,3)*fis_U239_rap*phi_rap(i-1,1)+ Y(i-1,4)*fis_Np239_rap*phi_rap(i-1,1) + Y(i-1,5)*fis_Pu239_rap*phi_rap(i-1,1))*2*0.95*dt_gen - Y(i-1,6)*lambda_PF*dt_gen;
    Y(i,7) = Y(i-1,7) + (Y(i-1,1)*fis_U235_th*phi_th(i-1,1) + Y(i-1,2)*fis_U238_th*phi_th(i-1,1) + Y(i-1,3)*fis_U239_th*phi_th(i-1,1) + Y(i-1,4)*fis_Np239_th*phi_th(i-1,1) + Y(i-1,5)*fis_Pu239_th*phi_th(i-1,1) + Y(i-1,1)*fis_U235_rap*phi_rap(i-1,1) + Y(i-1,2)*fis_U238_rap*phi_rap(i-1,1)+ Y(i-1,3)*fis_U239_rap*phi_rap(i-1,1)+ Y(i-1,4)*fis_Np239_rap*phi_rap(i-1,1) + Y(i-1,5)*fis_Pu239_rap*phi_rap(i-1,1))*2*0.05*dt_gen + (- Y(i-1,7)*cap_Xe135_th*phi_th(i-1,1) - Y(i-1,7)*cap_Xe135_rap*phi_rap(i-1,1) - Y(i-1,7)*log(2)/demi_Xe135)*dt_gen;
    Y(i,8) = Y(i-1,8) + (Y(i-1,6)*lambda_PF + Y(i-1,7)*cap_Xe135_th*phi_th(i-1,1) + Y(i-1,7)*cap_Xe135_rap*phi_rap(i-1,1) + Y(i-1,7)*log(2)/demi_Xe135)*dt_gen;
    
    n_th(i,1) = n_th(i-1,1) + (- Y(i-1,1)*fis_U235_th*phi_th(i-1,1) - Y(i-1,2)*fis_U238_th*phi_th(i-1,1) - Y(i-1,3)*fis_U239_th*phi_th(i-1,1) - Y(i-1,4)*fis_Np239_th*phi_th(i-1,1) - Y(i-1,5)*fis_Pu239_th*phi_th(i-1,1) - Y(i-1,2)*cap_U238_th*phi_th(i-1,1) + Y(i-1,6)*lambda_PF)*NA*dt_gen + (n_rap(i-1,1)*lambda_rt - n_th(i-1,1)*lambda_perte_th)*dt_gen;
    phi_th(i,1) = n_th(i,1)*v_thermique/V;
    n_rap(i,1) = n_rap(i-1,1) + (Y(i-1,1)*fis_U235_th*phi_th(i-1,1) + Y(i-1,2)*fis_U238_th*phi_th(i-1,1) + Y(i-1,3)*fis_U239_th*phi_th(i-1,1) + Y(i-1,4)*fis_Np239_th*phi_th(i-1,1) + Y(i-1,5)*fis_Pu239_th*phi_th(i-1,1))*2*NA*dt_gen + (Y(i-1,1)*fis_U235_rap*phi_rap(i-1,1) + Y(i-1,2)*fis_U238_rap*phi_rap(i-1,1) + Y(i-1,3)*fis_U239_rap*phi_rap(i-1,1) + Y(i-1,4)*fis_Np239_rap*phi_rap(i-1,1) + Y(i-1,5)*fis_Pu239_rap*phi_rap(i-1,1) - Y(i-1,2)*cap_U238_rap*phi_rap(i-1,1))*NA*dt_gen + (- n_rap(i-1,1)*lambda_rt - n_rap(i-1,1)*lambda_perte_rap)*dt_gen;
    phi_rap(i,1) = n_rap(i,1)*v_rapide/V;
    
%     n(i,1) = n_th(i,1) + n_rap(i,1); %tous les neutrons (thermique+rapide)
%     n_ret(i,1) = n_ret(i-1,1) + Y(i-1,6)*lambda_PF*NA*dt_gen; %neutrons retardes
  
end

%--------------------------------------------------------------------------
%%

%%%%%%%%%
%GRAPHES%
%%%%%%%%%

%especes
figure;
loglog(T,Y(:,1));
hold on;
loglog(T,Y(:,2));
hold on;
loglog(T,Y(:,3));
hold on;
loglog(T,Y(:,4));
hold on;
loglog(T,Y(:,5));
hold on;
loglog(T,Y(:,6),'--');
hold on;
loglog(T,Y(:,7),'--');
hold on;
loglog(T,Y(:,8),'--');
xlabel('Temps [s]');
ylabel('Espèces [mol]');
legend('U235','U238','U239','Np239','Pu239','PF*','Xe135','PF'); %,'Location','southeast'
hold off;

%neutrons
figure;
loglog(T,n_th(:,1));
hold on;
loglog(T,n_rap(:,1));
% hold on;
% loglog(T,n_ret(:,1));
xlabel('Temps [s]');
ylabel('Population de neutrons [#]');
legend('Thermique','Rapide');
hold off;

%flux de neutrons
figure;
loglog(T,phi_th(:,1));
hold on;
loglog(T,phi_rap(:,1));
xlabel('Temps [s]');
ylabel('Flux de neutrons [#/m^2.s]');
legend('Thermique','Rapide');
hold off;

%puissance


%--------------------------------------------------------------------------
%%
toc
end