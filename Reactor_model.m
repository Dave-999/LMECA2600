function [Lambda_BC_thermal Lambda_BC_fast U5_burning_rate] = Reactor_model(t_final,n_th_init ,n_fa_init ,m_tot,U5_pour, U8_pour,Pu9_pour,Poison_pour)
tic
%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Donnees et hypotheses %
%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    
    t_final = 100; %temps de simulation [s]
    n_th_init = 1e15; %nombre initial de neutrons thermiques
    n_fa_init = 0; %nombre intial de neutrons rapides
    m_tot = 25000; %masse totale de combustible [kg]
    U5_pour = 3; %pourcentage massique d'U235 dans le combustible [%]
    U8_pour = 97; %pourcentage massique d'U238 dans le combustible [%]
    Pu9_pour = 0; %pourcentage massique de Pu239 dans le combustible [%]
    Poison_pour = 5; %pourcentage molaire de PF*_poison dans les produits de fission [%]
    
end

NA = 6.02214076e23; %nombre d'Avogadro
m_neutron = 1.67493e-27; %masse neutron [kg]
eV_Joule = 1.60218e-19; %conversion eV en J
W_GW = 1e-9; %conversion W en GW
barn_m2 = 1e-28; %conversion barn en m2
V = 10; %volume du reacteur [m^3]

U235_init = m_tot*U5_pour/100/molarMass('U235'); %nombre initial de moles de U235 [mol/l]
U238_init = m_tot*U8_pour/100/molarMass('U238'); %nombre initial de moles de U238 [mol/l]
Pu239_init = m_tot*Pu9_pour/100/molarMass('Pu239'); %nombre initial de moles de Pu239 [mol/l]
U239_init = 0; %nombre initial de moles de Pu239 [mol/l]
Np239_init = 0; %nombre initial de moles de Np239 [mol/l]
PF_star_init = 0; %nombre initial de moles de PF* [mol/l]
PF_star_poison_init = 0; %nombre initial de moles de PF*_poison [mol/l]
PF_init = 0; %nombre initial de moles de PF [mol/l]

E_th = 0.025; %energie neutron thermique [eV]
E_rap = 1e6; %energie neutron rapide [eV]
v_th = sqrt(E_th*eV_Joule*2/m_neutron); %vitesse neutron thermique [m/s]
v_rap = sqrt(E_rap*eV_Joule*2/m_neutron); %vitesse neutron rapide [m/s]

T_PF = 7; %temps de demi-vie pour PF* = PF + n + 5 MeV
lambda_PF = log(2)/T_PF; %lambda pour PF* = PF + n + 5 MeV
lambda_rt = log(2)/(5e-4); %lambda transition rapide->thermique
Lambda_BC_thermal = 130; %lambda thermique pour les barres de controle
Lambda_BC_fast = 260; %lambda rapide pour les barres de controle

E_fis = 200e6 * eV_Joule * W_GW; %energie fission [GJ]
E_PF = 5e6 * eV_Joule * W_GW; %energie PF* = PF + n [GJ]
E_rt = 1e6 * eV_Joule * W_GW; %energie transition rapide->thermique [GJ]

Power_init = 0; %puissance initiale [GW]
P_min = 1; %puissance min [GW]
P_max = 3; %puissance max [GW]

%--------------------------------------------------------------------------
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculs des sections efficaces et des temps de demi-vie %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%U235
U235_sig_fis_th = barn_m2*Section_efficace('U235','Fission',E_th,'DATABASE');
U235_sig_fis_rap = barn_m2*Section_efficace('U235','Fission',E_rap,'DATABASE');
    
%U238
U238_sig_fis_th = barn_m2*Section_efficace('U238','Fission',E_th,'DATABASE');
U238_sig_fis_rap = barn_m2*Section_efficace('U238','Fission',E_rap,'DATABASE');
U238_sig_cap_th = barn_m2*Section_efficace('U238','Capture',E_th,'DATABASE');
U238_sig_cap_rap = barn_m2*Section_efficace('U238','Capture',E_rap,'DATABASE');

%U239
U239_sig_fis_th = barn_m2*Section_efficace('U239','Fission',E_th,'DATABASE');
U239_sig_fis_rap = barn_m2*Section_efficace('U239','Fission',E_rap,'DATABASE');
U239_demi_vie = Demi_vie('U239','BetaMinus');
    
%Np239
Np239_sig_fis_th = barn_m2*Section_efficace('Np239','Fission',E_th,'DATABASE');
Np239_sig_fis_rap = barn_m2*Section_efficace('Np239','Fission',E_rap,'DATABASE');
Np239_demi_vie = Demi_vie('Np239','BetaMinus');
    
%Pu239
Pu239_sig_fis_th = barn_m2*Section_efficace('Pu239','Fission',E_th,'DATABASE');
Pu239_sig_fis_rap = barn_m2*Section_efficace('Pu239','Fission',E_rap,'DATABASE');

%Xe135
Xe135_sig_cap_th = barn_m2*Section_efficace('Xe135','Capture',E_th,'DATABASE');
Xe135_sig_cap_rap = barn_m2*Section_efficace('Xe135','Capture',E_rap,'DATABASE');
Xe135_demi_vie = Demi_vie('Xe135','BetaMinus');

%--------------------------------------------------------------------------
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation des vecteurs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt_gen = 1e-4;
T = (0:dt_gen:t_final);

Y = zeros(length(T),8); %U235,U238,U239,Np239,Pu239,PF*,PF*_poison,PF
Y(1,:) = [U235_init U238_init U239_init Np239_init Pu239_init PF_star_init PF_star_poison_init PF_init]; %quantites initiales [mol]

n_th = zeros(length(T),1); %neutrons thermiques 
n_th(1,1) = n_th_init; %nombre initial [#]
n_rap = zeros(length(T),1); %neutrons rapides
n_rap(1,1) = n_fa_init; %nombre initial [#]

phi_th = zeros(length(T),1); %flux de neutrons thermiques
phi_th(1,1) = n_th(1,1)*v_th/V; %flux initial [#/m^2.s]
phi_rap = zeros(length(T),1); %flux de neutrons rapides
phi_rap(1,1) = n_rap(1,1)*v_rap/V; %flux initial [#/m^2.s]

Power = zeros(length(T),1); %puissance [GW]
Power(1,1) = Power_init; %puissance initiale [GW]

%--------------------------------------------------------------------------
%%

%%%%%%%%%%%%%%%%%%%%%%%%
% Equations cinetiques %
%%%%%%%%%%%%%%%%%%%%%%%%

%Equations detaillees dans le rapport

for i = 2:length(T)
    
    for j = 1:1
        
        Y(i,1) = Y(i-1,1) + (- Y(i-1,1)*U235_sig_fis_th*phi_th(i-1,1) - Y(i-1,1)*U235_sig_fis_rap*phi_rap(i-1,1))*dt_gen;
        Y(i,2) = Y(i-1,2) + (- Y(i-1,2)*U238_sig_fis_th*phi_th(i-1,1) - Y(i-1,2)*U238_sig_fis_rap*phi_rap(i-1,1) - Y(i-1,2)*U238_sig_cap_th*phi_th(i-1,1) - Y(i-1,2)*U238_sig_cap_rap*phi_rap(i-1,1))*dt_gen;
        Y(i,3) = Y(i-1,3) + (Y(i-1,2)*U238_sig_cap_th*phi_th(i-1,1) + Y(i-1,2)*U238_sig_cap_rap*phi_rap(i-1,1) - Y(i-1,3)*U239_sig_fis_th*phi_th(i-1,1) - Y(i-1,3)*U239_sig_fis_rap*phi_rap(i-1,1) - Y(i-1,3)*log(2)/U239_demi_vie)*dt_gen;
        Y(i,4) = Y(i-1,4) + (Y(i-1,3)*log(2)/U239_demi_vie - Y(i-1,4)*Np239_sig_fis_th*phi_th(i-1,1) - Y(i-1,4)*Np239_sig_fis_rap*phi_rap(i-1,1) - Y(i-1,4)*log(2)/Np239_demi_vie)*dt_gen;
        Y(i,5) = Y(i-1,5) + (Y(i-1,4)*log(2)/Np239_demi_vie - Y(i-1,5)*Pu239_sig_fis_th*phi_th(i-1,1) - Y(i-1,5)*Pu239_sig_fis_rap*phi_rap(i-1,1))*dt_gen;
        Y(i,6) = Y(i-1,6) + (Y(i-1,1)*U235_sig_fis_th*phi_th(i-1,1) + Y(i-1,2)*U238_sig_fis_th*phi_th(i-1,1) + Y(i-1,3)*U239_sig_fis_th*phi_th(i-1,1) + Y(i-1,4)*Np239_sig_fis_th*phi_th(i-1,1) + Y(i-1,5)*Pu239_sig_fis_th*phi_th(i-1,1) + Y(i-1,1)*U235_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,2)*U238_sig_fis_rap*phi_rap(i-1,1)+ Y(i-1,3)*U239_sig_fis_rap*phi_rap(i-1,1)+ Y(i-1,4)*Np239_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,5)*Pu239_sig_fis_rap*phi_rap(i-1,1))*2*(1-Poison_pour/100)*dt_gen - Y(i-1,6)*lambda_PF*dt_gen;
        Y(i,7) = Y(i-1,7) + (Y(i-1,1)*U235_sig_fis_th*phi_th(i-1,1) + Y(i-1,2)*U238_sig_fis_th*phi_th(i-1,1) + Y(i-1,3)*U239_sig_fis_th*phi_th(i-1,1) + Y(i-1,4)*Np239_sig_fis_th*phi_th(i-1,1) + Y(i-1,5)*Pu239_sig_fis_th*phi_th(i-1,1) + Y(i-1,1)*U235_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,2)*U238_sig_fis_rap*phi_rap(i-1,1)+ Y(i-1,3)*U239_sig_fis_rap*phi_rap(i-1,1)+ Y(i-1,4)*Np239_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,5)*Pu239_sig_fis_rap*phi_rap(i-1,1))*2*Poison_pour/100*dt_gen + (- Y(i-1,7)*Xe135_sig_cap_th*phi_th(i-1,1) - Y(i-1,7)*Xe135_sig_cap_rap*phi_rap(i-1,1) - Y(i-1,7)*log(2)/Xe135_demi_vie)*dt_gen;
        Y(i,8) = Y(i-1,8) + (Y(i-1,6)*lambda_PF + Y(i-1,7)*Xe135_sig_cap_th*phi_th(i-1,1) + Y(i-1,7)*Xe135_sig_cap_rap*phi_rap(i-1,1) + Y(i-1,7)*log(2)/Xe135_demi_vie)*dt_gen;
        
        n_th(i,1) = n_th(i-1,1) + (- Y(i-1,1)*U235_sig_fis_th*phi_th(i-1,1) - Y(i-1,2)*U238_sig_fis_th*phi_th(i-1,1) - Y(i-1,3)*U239_sig_fis_th*phi_th(i-1,1) - Y(i-1,4)*Np239_sig_fis_th*phi_th(i-1,1) - Y(i-1,5)*Pu239_sig_fis_th*phi_th(i-1,1) - Y(i-1,2)*U238_sig_cap_th*phi_th(i-1,1) + Y(i-1,6)*lambda_PF)*NA*dt_gen + (n_rap(i-1,1)*lambda_rt - n_th(i-1,1)*Lambda_BC_thermal)*dt_gen;
        phi_th(i,1) = n_th(i,1)*v_th/V;
        n_rap(i,1) = n_rap(i-1,1) + (Y(i-1,1)*U235_sig_fis_th*phi_th(i-1,1) + Y(i-1,2)*U238_sig_fis_th*phi_th(i-1,1) + Y(i-1,3)*U239_sig_fis_th*phi_th(i-1,1) + Y(i-1,4)*Np239_sig_fis_th*phi_th(i-1,1) + Y(i-1,5)*Pu239_sig_fis_th*phi_th(i-1,1))*2*NA*dt_gen + (Y(i-1,1)*U235_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,2)*U238_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,3)*U239_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,4)*Np239_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,5)*Pu239_sig_fis_rap*phi_rap(i-1,1) - Y(i-1,2)*U238_sig_cap_rap*phi_rap(i-1,1))*NA*dt_gen + (- n_rap(i-1,1)*lambda_rt - n_rap(i-1,1)*Lambda_BC_fast)*dt_gen;
        phi_rap(i,1) = n_rap(i,1)*v_rap/V;
        
        Power(i,1) = (Y(i-1,1)*U235_sig_fis_th*phi_th(i-1,1) + Y(i-1,2)*U238_sig_fis_th*phi_th(i-1,1) + Y(i-1,3)*U239_sig_fis_th*phi_th(i-1,1) + Y(i-1,4)*Np239_sig_fis_th*phi_th(i-1,1) + Y(i-1,5)*Pu239_sig_fis_th*phi_th(i-1,1) + Y(i-1,1)*U235_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,2)*U238_sig_fis_rap*phi_rap(i-1,1)+ Y(i-1,3)*U239_sig_fis_rap*phi_rap(i-1,1)+ Y(i-1,4)*Np239_sig_fis_rap*phi_rap(i-1,1) + Y(i-1,5)*Pu239_sig_fis_rap*phi_rap(i-1,1))*NA*E_fis + Y(i-1,6)*lambda_PF*NA*E_PF + n_rap(i-1,1)*lambda_rt*E_rt;
    
%         phi_th(i-1,1) = (phi_th(i-1,1) + phi_th(i,1))/2;
%         phi_rap(i-1,1) = (phi_rap(i-1,1) + phi_rap(i,1))/2;
        
        if mod(T(i),1) == 0
            if Power(i,1) < P_min
                Lambda_BC_thermal = Lambda_BC_thermal * 0.95;
                Lambda_BC_fast = Lambda_BC_fast * 0.95;
            end
            if Power(i,1) >= P_min
                if Power(i,1)-Power(i-1e4,1) > 0
                    Lambda_BC_thermal = Lambda_BC_thermal * (1 + 0.6*(Power(i,1)-Power(i-1e4,1))/Power(i,1));
                    Lambda_BC_fast = Lambda_BC_fast * (1 + 0.6*(Power(i,1)-Power(i-1e4,1))/Power(i,1));
                elseif Power(i,1)-Power(i-1e4,1) < 0
                    Lambda_BC_thermal = Lambda_BC_thermal * (1 - 0.35*(Power(i-1e4,1)-Power(i,1))/Power(i-1e4,1));
                    Lambda_BC_fast = Lambda_BC_fast * (1 - 0.35*(Power(i-1e4,1)-Power(i,1))/Power(i-1e4,1));
                end
            end
        end
        
    end %fin boucle it�ration flux
    
end %fin boucle lenght(T)



%--------------------------------------------------------------------------
%%

%%%%%%%%%%%
% Graphes %
%%%%%%%%%%%

%especes
figure;
loglog(T,Y(:,1)/V);
hold on;
loglog(T,Y(:,2)/V);
hold on;
loglog(T,Y(:,3)/V);
hold on;
loglog(T,Y(:,4)/V);
hold on;
loglog(T,Y(:,5)/V);
hold on;
loglog(T,Y(:,6)/V,'--');
hold on;
loglog(T,Y(:,7)/V,'k--');
hold on;
loglog(T,Y(:,8)/V,'--');
xlabel('Temps [s]');
ylabel('Concentration molaire des esp�ces [mol/l]');
legend('U235','U238','U239','Np239','Pu239','PF^{*}','PF^{*}_poison','PF');
hold off;

% %neutrons
% figure;
% loglog(T,n_th(:,1),'b');
% hold on;
% loglog(T,n_rap(:,1),'r');
% xlabel('Temps [s]');
% ylabel('Population de neutrons [#]');
% legend('Thermique','Rapide');
% hold off;

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
ylabel('Puissance lib�r�e [GW]');
legend('P_{MIN}','P_{MAX}','POWER');
hold off;

%--------------------------------------------------------------------------
%%
toc
end