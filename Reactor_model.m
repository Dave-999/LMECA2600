function [] = Reactor_model(t_final,dt_plot,P_stable,PF_retarded)
tic
close all;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Donnees et hypotheses %
%%%%%%%%%%%%%%%%%%%%%%%%%

NA = 6.02214076e23; %avogadro
V = 10; %volume du reacteur [m^3]

m_Utot = 1000; %masse d'uranium [kg]
m_U235 = m_Utot*0.03; %masse d'uranium 235
m_U238 = m_Utot*0.97; %masse d'uranium 238

N_U235 = m_U235/molarMass('U235'); %nombre de moles d'U235 [mol]
N_U238 = m_U238/molarMass('U238'); %nombre de moles d'U238 [mol]

n_0 = 1e15; %nombre initial de neutrons (thermiques)

E_n_th = 0.025; %energie neutron thermique [eV]
E_n_rap = 1e6; %energie neutron rapide [eV]

v_n_th = 10; %vitesse neutron thermique [m/s]
v_n_rap = 1000; %vitesse neutron rapide [m/s]

lambda_PF = log(2)/1; %lambda pour PF* = PF + n + 5 MeV
lambda_rt = log(2)/(5*10^-4); %lambda transition rapide->thermique
lambda_perte_th = 800; %lambda thermique pour fuites et barres de controle
lambda_perte_rap = 2000; %lambda rapide pour fuites et barres de controle

E_fis = 200e6 * 1.602176565e-19 * 1e-9; %energie fission [GJ]
E_PF = 5e6 * 1.602176565e-19 * 1e-9; %energie PF* = PF + n [GJ]
E_rt = 1e6 * 1.602176565e-19 * 1e-9; %energie transition rapide->thermique [GJ]

P_min = 1; %puissance min [GW]
P_max = 3; %puissance max [GW]
P_objectif = 2; %puissance désirée [GW]

Time = 1; %modification de lambda_perte après "Time" seconde(s)

%--------------------------------------------------------------------------
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculs des sections efficaces et des temps de demi-vie %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%U235
U235_sig_fis_th = 1e-28*Section_efficace('U235','Fission',E_n_th,'DATABASE');
U235_sig_fis_rap = 1e-28*Section_efficace('U235','Fission',E_n_rap,'DATABASE');
    
%U238
U238_sig_fis_th = 1e-28*Section_efficace('U238','Fission',E_n_th,'DATABASE');
U238_sig_fis_rap = 1e-28*Section_efficace('U238','Fission',E_n_rap,'DATABASE');
U238_sig_cap_th = 1e-28*Section_efficace('U238','Capture',E_n_th,'DATABASE');
U238_sig_cap_rap = 1e-28*Section_efficace('U238','Capture',E_n_rap,'DATABASE');

%U239
U239_sig_fis_th = 1e-28*Section_efficace('U239','Fission',E_n_th,'DATABASE');
U239_sig_fis_rap = 1e-28*Section_efficace('U239','Fission',E_n_rap,'DATABASE');
U239_demi_vie = Demi_vie('U239','BetaMinus');
    
%Np239
Np239_sig_fis_th = 1e-28*Section_efficace('Np239','Fission',E_n_th,'DATABASE');
Np239_sig_fis_rap = 1e-28*Section_efficace('Np239','Fission',E_n_rap,'DATABASE');
Np239_demi_vie = Demi_vie('Np239','BetaMinus');
    
%Pu239
Pu239_sig_fis_th = 1e-28*Section_efficace('Pu239','Fission',E_n_th,'DATABASE');
Pu239_sig_fis_rap = 1e-28*Section_efficace('Pu239','Fission',E_n_rap,'DATABASE');

%Xe135
Xe135_sig_cap_th = 1e-28*Section_efficace('Xe135','Capture',E_n_th,'DATABASE');
Xe135_sig_cap_rap = 1e-28*Section_efficace('Xe135','Capture',E_n_rap,'DATABASE');
Xe135_demi_vie = Demi_vie('Xe135','BetaMinus');

%--------------------------------------------------------------------------
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation des vecteurs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_final = 100; %[s]
dt_gen = 1e-4;
T = [0:dt_gen:t_final];

Y = zeros(length(T),8); %U235,U238,U239,Np239,Pu239,PF*,Xe135,PF
Y(1,:) = [N_U235 N_U238 0 0 0 0 0 0]; %quantites initiales [mol]

n_th = zeros(length(T),1); %neutrons thermiques 
n_th(1,1) = n_0; %nombre initial [#]
n_rap = zeros(length(T),1); %neutrons rapides
n_rap(1,1) = 0; %nombre initial [#]

phi_th = zeros(length(T),1); %flux de neutrons thermiques
phi_th(1,1) = n_th(1,1)*v_n_th/V; %flux initial [#/m^2.s]
phi_rap = zeros(length(T),1); %flux de neutrons rapides
phi_rap(1,1) = n_rap(1,1)*v_n_rap/V; %flux initial [#/m^2.s]

Power = zeros(length(T),1); %puissance [GW]
Power(1,1) = 0;

%--------------------------------------------------------------------------
%%

%%%%%%%%%%%%%%%%%%%%%%%%
% Equations cinetiques %
%%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:length(T)
    
    Y(i,1) = Y(i-1,1) + (- Y(i-1,1)*U235_sig_fis_th*phi_th(i-1,1) - Y(i-1,1)*U235_sig_fis_rap*phi_rap(i-1,1))*dt_gen;
    Y(i,2) = Y(i-1,2) + (- Y(i-1,2)*U238_sig_fis_th*phi_th(i-1,1) - Y(i-1,2)*U238_sig_fis_rap*phi_rap(i-1,1) - Y(i-1,2)*U238_sig_cap_th*phi_th(i-1,1) - Y(i-1,2)*U238_sig_cap_rap*phi_rap(i-1,1))*dt_gen;
    Y(i,3) = Y(i-1,3) + (Y(i-1,2)*U238_sig_cap_th*phi_th(i-1,1) + Y(i-1,2)*U238_sig_cap_rap*phi_rap(i-1,1) - Y(i-1,3)*U239_sig_fis_th*phi_th(i-1,1) - Y(i-1,3)*U239_sig_fis_rap*phi_rap(i-1,1) - Y(i-1,3)*log(2)/U239_demi_vie)*dt_gen;
    Y(i,4) = Y(i-1,4) + (Y(i-1,3)*log(2)/U239_demi_vie - Y(i-1,4)*Np239_sig_fis_th*phi_th(i-1,1) - Y(i-1,4)*Np239_sig_fis_rap*phi_rap(i-1,1) - Y(i-1,4)*log(2)/Np239_demi_vie)*dt_gen; 
    Y(i,5) = Y(i-1,5) + (Y(i-1,4)*log(2)/Np239_demi_vie - Y(i-1,5)*Pu239_sig_fis_th*phi_th(i-1,1) - Y(i-1,5)*Pu239_sig_fis_rap*phi_rap(i-1,1))*dt_gen;
    Y(i,6) = Y(i-1,6) + (Y(i-1,1)*U235_sig_fis_th*phi_th(i-1,1) + Y(i-1,2)*U238_sig_fis_th*phi_th(i-1,1) + Y(i-1,3)*U239_sig_fis_th*phi_th(i-1,1) + Y(i-1,4)*Np239_sig_fis_th*phi_th(i-1,1) + Y(i-1,5)*Pu239_sig_fis_th*phi_th(i-1,1) + Y(i-1,1)*U235_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,2)*U238_sig_fis_rap*phi_rap(i-1,1)+ Y(i-1,3)*U239_sig_fis_rap*phi_rap(i-1,1)+ Y(i-1,4)*Np239_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,5)*Pu239_sig_fis_rap*phi_rap(i-1,1))*2*0.95*dt_gen - Y(i-1,6)*lambda_PF*dt_gen;
    Y(i,7) = Y(i-1,7) + (Y(i-1,1)*U235_sig_fis_th*phi_th(i-1,1) + Y(i-1,2)*U238_sig_fis_th*phi_th(i-1,1) + Y(i-1,3)*U239_sig_fis_th*phi_th(i-1,1) + Y(i-1,4)*Np239_sig_fis_th*phi_th(i-1,1) + Y(i-1,5)*Pu239_sig_fis_th*phi_th(i-1,1) + Y(i-1,1)*U235_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,2)*U238_sig_fis_rap*phi_rap(i-1,1)+ Y(i-1,3)*U239_sig_fis_rap*phi_rap(i-1,1)+ Y(i-1,4)*Np239_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,5)*Pu239_sig_fis_rap*phi_rap(i-1,1))*2*0.05*dt_gen + (- Y(i-1,7)*Xe135_sig_cap_th*phi_th(i-1,1) - Y(i-1,7)*Xe135_sig_cap_rap*phi_rap(i-1,1) - Y(i-1,7)*log(2)/Xe135_demi_vie)*dt_gen;
    Y(i,8) = Y(i-1,8) + (Y(i-1,6)*lambda_PF + Y(i-1,7)*Xe135_sig_cap_th*phi_th(i-1,1) + Y(i-1,7)*Xe135_sig_cap_rap*phi_rap(i-1,1) + Y(i-1,7)*log(2)/Xe135_demi_vie)*dt_gen;
    
    n_th(i,1) = n_th(i-1,1) + (- Y(i-1,1)*U235_sig_fis_th*phi_th(i-1,1) - Y(i-1,2)*U238_sig_fis_th*phi_th(i-1,1) - Y(i-1,3)*U239_sig_fis_th*phi_th(i-1,1) - Y(i-1,4)*Np239_sig_fis_th*phi_th(i-1,1) - Y(i-1,5)*Pu239_sig_fis_th*phi_th(i-1,1) - Y(i-1,2)*U238_sig_cap_th*phi_th(i-1,1) + Y(i-1,6)*lambda_PF)*NA*dt_gen + (n_rap(i-1,1)*lambda_rt - n_th(i-1,1)*lambda_perte_th)*dt_gen;
    phi_th(i,1) = n_th(i,1)*v_n_th/V;
    n_rap(i,1) = n_rap(i-1,1) + (Y(i-1,1)*U235_sig_fis_th*phi_th(i-1,1) + Y(i-1,2)*U238_sig_fis_th*phi_th(i-1,1) + Y(i-1,3)*U239_sig_fis_th*phi_th(i-1,1) + Y(i-1,4)*Np239_sig_fis_th*phi_th(i-1,1) + Y(i-1,5)*Pu239_sig_fis_th*phi_th(i-1,1))*2*NA*dt_gen + (Y(i-1,1)*U235_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,2)*U238_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,3)*U239_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,4)*Np239_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,5)*Pu239_sig_fis_rap*phi_rap(i-1,1) - Y(i-1,2)*U238_sig_cap_rap*phi_rap(i-1,1))*NA*dt_gen + (- n_rap(i-1,1)*lambda_rt - n_rap(i-1,1)*lambda_perte_rap)*dt_gen;
    phi_rap(i,1) = n_rap(i,1)*v_n_rap/V;
    
    Power(i,1) = (Y(i-1,1)*U235_sig_fis_th*phi_th(i-1,1) + Y(i-1,2)*U238_sig_fis_th*phi_th(i-1,1) + Y(i-1,3)*U239_sig_fis_th*phi_th(i-1,1) + Y(i-1,4)*Np239_sig_fis_th*phi_th(i-1,1) + Y(i-1,5)*Pu239_sig_fis_th*phi_th(i-1,1) + Y(i-1,1)*U235_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,2)*U238_sig_fis_rap*phi_rap(i-1,1)+ Y(i-1,3)*U239_sig_fis_rap*phi_rap(i-1,1)+ Y(i-1,4)*Np239_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,5)*Pu239_sig_fis_rap*phi_rap(i-1,1))*NA*E_fis + Y(i-1,6)*lambda_PF*NA*E_PF + n_rap(i-1,1)*lambda_rt*E_rt;
    
    if mod((i-1)*dt_gen/Time,1) == 0 %on verifie si le temps requis s'est ecoule
        if Power(i,1) < P_min
            lambda_perte_th = lambda_perte_th * 0.95;
            lambda_perte_rap = lambda_perte_rap * 0.95;
        elseif Power(i,1) >= P_min && Power(i,1) < P_max
            
            if Power(i,1) < P_objectif
                
                if (Power(i,1) - Power(i-1/dt_gen*Time,1)) > 0
                    lambda_perte_th = lambda_perte_th / (0.95 * (P_objectif-Power(i,1))/Power(i,1));
                    lambda_perte_rap = lambda_perte_rap / (0.95 * (P_objectif-Power(i,1))/Power(i,1));
                elseif (Power(i,1) - Power(i-1/dt_gen*Time,1)) <= 0
                    lambda_perte_th = lambda_perte_th * (0.95 * (P_objectif-Power(i,1))/Power(i,1));
                    lambda_perte_rap = lambda_perte_rap * (0.95 * (P_objectif-Power(i,1))/Power(i,1));
                end
                
            elseif Power(i,1) > P_objectif
                
                if (Power(i,1) - Power(i-1/dt_gen*Time,1)) > 0
                    lambda_perte_th = lambda_perte_th / (0.95 * (P_objectif-Power(i,1))/Power(i,1));
                    lambda_perte_rap = lambda_perte_rap / (0.95 * (P_objectif-Power(i,1))/Power(i,1));
                elseif (Power(i,1) - Power(i-1/dt_gen*Time,1)) <= 0
                    lambda_perte_th = lambda_perte_th * (0.95 * (P_objectif-Power(i,1))/Power(i,1));
                    lambda_perte_rap = lambda_perte_rap * (0.95 * (P_objectif-Power(i,1))/Power(i,1));
                end
                
            end
            
        elseif Power(i,1) >= P_max
            lambda_perte_th = lambda_perte_th / 0.95;
            lambda_perte_rap = lambda_perte_rap / 0.95;
        end
    end
    
end

%--------------------------------------------------------------------------
%%

%%%%%%%%%%%
% Graphes %
%%%%%%%%%%%

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
loglog(T,Y(:,7),'k--');
hold on;
loglog(T,Y(:,8),'--');
xlabel('Temps [s]');
ylabel('Espèces [mol]');
legend('U235','U238','U239','Np239','Pu239','PF*','Xe135','PF');
hold off;

%neutrons
figure;
loglog(T,n_th(:,1),'b');
hold on;
loglog(T,n_rap(:,1),'r');
xlabel('Temps [s]');
ylabel('Population de neutrons [#]');
legend('Thermique','Rapide');
hold off;

%flux de neutrons
figure;
loglog(T,phi_th(:,1),'b');
hold on;
loglog(T,phi_rap(:,1),'r');
xlabel('Temps [s]');
ylabel('Flux de neutrons [#/m^2.s]');
legend('Thermique','Rapide');
hold off;

%puissance
figure;
plot([0 t_final],[P_min P_min],'b');
hold on;
plot([0 t_final],[P_max P_max],'r');
hold on;
plot(T,Power(:,1),'g');
xlabel('Temps [s]');
ylabel('Puissance libérée [GW]');
legend('P_{MIN}','P_{MAX}','POWER');
hold off;

%--------------------------------------------------------------------------
%%
toc
end