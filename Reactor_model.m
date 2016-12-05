function m = Reactor_model(t_final,dt_plot,P_stable,PF_retarded)
tic
%Donnees

NA = 6.022e23;
V = 10; %Volume du reacteur en [m^3]

m_Utot = 1000; %[kg]
m_U235 = m_Utot*0.03;
m_U238 = m_Utot*0.97;

N_U235 = m_U235/molarMass('U235'); %Nombre de moles d'U235 [mol]
N_U238 = m_U238/molarMass('U238'); %Nombre de moles d'U238 [mol]

E_thermal = 0.025; %[eV]
E_fast = 1e6; %[eV]

v_thermal = 10; %[m/s]
v_fast= 1000; %[m/s]

n_thermal = 1e10; %Nombre de neutrons thermiques en t=0
flux_thermal = n_thermal*v_thermal/V; %Flux de neutrons thermiques en t=0 [#/m^2.s]

lambda_neutron = 5; %Lambda pour les fuites de neutrons et barres de controle
T_PF = 1;
lambda_PF = log(2)/T_PF; %Lambda pour PF* = PF + n + 5 MeV

%--------------------------------------------------------------------------
%%

fis_U235 = 1e-28*Section_efficace('U235','Fission',E_thermal,'DATABASE');
    
fis_U238 = 1e-28*Section_efficace('U238','Fission',E_thermal,'DATABASE');
cap_U238 = 1e-28*Section_efficace('U238','Capture',E_thermal,'DATABASE');

fis_U239 = 1e-28*Section_efficace('U239','Fission',E_thermal,'DATABASE');
demi_U239 = Demi_vie('U239','BetaMinus');
    
fis_Np239 = 1e-28*Section_efficace('Np239','Fission',E_thermal,'DATABASE');
demi_Np239 = Demi_vie('Np239','BetaMinus');
    
fis_Pu239 = 1e-28*Section_efficace('Pu239','Fission',E_thermal,'DATABASE');


t_final = 10; %[s]
dt_gen = 10^-4;
T = [0:dt_gen:t_final];
Y = zeros(length(T),7); %U235,U238,U239,Np239,Pu239,PF*,PF
Y(1,:) = [N_U235 N_U238 0 0 0 0 0]; %Quantites initiales [mol]
N = zeros(length(T),1); %Flux de neutrons thermiques
N(1,1) = n_thermal; %Nb de neutrons initial [#]
Phi = zeros(length(T),1);
Phi(1,1) = flux_thermal; %Flux initial [#/m^2.s]

for i = 2:length(T)
    Y(i,1) = Y(i-1,1) + (- Y(i-1,1)*fis_U235*Phi(i-1,1))*dt_gen; %U235
    Y(i,2) = Y(i-1,2) + (- Y(i-1,2)*fis_U238*Phi(i-1,1) - Y(i-1,2)*cap_U238*Phi(i-1,1))*dt_gen; %U238
    Y(i,3) = Y(i-1,3) + (Y(i-1,2)*cap_U238*Phi(i-1,1) - Y(i-1,3)*fis_U239*Phi(i-1,1) - Y(i-1,3)*log(2)/demi_U239)*dt_gen; %U239
    Y(i,4) = Y(i-1,4) + (Y(i-1,3)*log(2)/demi_U239 - Y(i-1,4)*fis_Np239*Phi(i-1,1) - Y(i-1,4)*log(2)/demi_Np239)*dt_gen; %Np239
    Y(i,5) = Y(i-1,5) + (Y(i-1,4)*log(2)/demi_Np239 - Y(i-1,5)*fis_Pu239*Phi(i-1,1))*dt_gen; %Pu239
    Y(i,6) = Y(i-1,6) + (Y(i-1,1)*fis_U235*Phi(i-1,1) + Y(i-1,2)*fis_U238*Phi(i-1,1) + Y(i-1,3)*fis_U239*Phi(i-1,1) + Y(i-1,4)*fis_Np239*Phi(i-1,1) + Y(i-1,5)*fis_Pu239*Phi(i-1,1))*2*dt_gen - Y(i-1,6)*lambda_PF*dt_gen; %PF*
    Y(i,7) = Y(i-1,7) + Y(i-1,6)*lambda_PF*dt_gen;
    
    N(i,1) = N(i-1,1) + (Y(i-1,1)*fis_U235*Phi(i-1,1) + Y(i-1,2)*fis_U238*Phi(i-1,1) - Y(i-1,2)*cap_U238*Phi(i-1,1) + Y(i-1,3)*fis_U239*Phi(i-1,1) + Y(i-1,4)*fis_Np239*Phi(i-1,1) + Y(i-1,5)*fis_Pu239*Phi(i-1,1) + Y(i-1,6)*lambda_PF)*NA*dt_gen - N(i-1,1)*lambda_neutron*dt_gen; %Nombres de neutrons (thermiques)
    Phi(i,1) = N(i,1)*v_thermal/V; %flux thermique
end

%GRAPHES
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
xlabel('Temps [s]');
ylabel('Espèces [#]');
legend('U235','U238','U239','Np239','Pu239','PF*','PF','Location','southeast');
hold off;

figure;
loglog(T,Phi(:,1));
xlabel('Temps [s]');
ylabel('Flux de neutrons thermiques [#/m^2.s]');


toc
end