function m = Reactor_model(t_final,dt_plot,P_stable,PF_retarded)
%Donnees
tic
V = 30; %Volume du reacteur en [m^3]
m_Utot = 1000; %Masse d'Uranium [kg]

m_U235 = m_Utot*0.07; %Masse d'U235 [kg]
N_U235 = m_U235/molarMass('U235'); %Nombre de moles d'U235 [mol]

m_U238 = m_Utot*0.93; %Masse d'U238 [kg]
N_U238 = m_U238/molarMass('U238'); %Nombre de moles d'U238 [mol]

E_thermal = 0.025; %[eV]
E_fast = 1e6; %[eV]

v_thermal = 10; %[m/s]
v_fast= 1000; %[m/s]

n_thermal = 1e10; %Nombre de neutrons thermiques en t=0
flux_thermal = n_thermal*v_thermal/V; %Flux de neutrons thermiques en t=0 [#/m^2/s]

%--------------------------------------------------------------------------
%%

%ALTERNATIVE RAPIDE, INCORRECTE OU NON?
sigma_U235 = Section_efficace('U235','Fission',E_thermal,'DATABASE');
    demi_U235 = Demi_vie('U235','Alpha');
    
    sigma_U238 = Section_efficace('U238','Capture',E_thermal,'DATABASE');
    demi_U238 = Demi_vie('U238','Alpha');
    
    demi_U239 = Demi_vie('U239','BetaMinus');
    
    demi_Np239 = Demi_vie('Np239','BetaMinus');
    
    sigma_Pu239 = Section_efficace('Pu239','Fission',E_thermal,'DATABASE');
    demi_Pu239 = Demi_vie('Pu239','Alpha');


t_final = 100; %[s]
dt=10^-4;
T=linspace(0,t_final,t_final/dt);
Y=zeros(length(T),6);
Y(1,:)=[N_U235 N_U238 0 0 0 0];

for i= 2:length(T)
    Y(i,1) = Y(i-1,1)+(- Y(i-1,1)*sigma_U235*10e-28*flux_thermal - Y(i-1,1)*log(2)/demi_U235)*dt; %U235
    Y(i,2) = Y(i-1,2)+(- Y(i-1,2)*sigma_U238*10e-28*flux_thermal - Y(i-1,2)*log(2)/demi_U238)*dt; %U238
    Y(i,3) = Y(i-1,3)+(Y(i-1,2)*sigma_U238*10e-28*flux_thermal - Y(i-1,3)*log(2)/demi_U239)*dt; %U239
    Y(i,4) = Y(i-1,4)+(Y(i-1,3)*log(2)/demi_U239 - Y(i-1,4)*log(2)/demi_Np239)*dt; %Np239
    Y(i,5) = Y(i-1,5)+(Y(i-1,4)*log(2)/demi_Np239 - Y(i-1,5)*sigma_U238*10e-28*flux_thermal - Y(i-1,5)*log(2)/demi_Pu239)*dt; %Pu239
    Y(i,6) = Y(i-1,6)+(Y(i-1,5)*sigma_U238*10e-28*flux_thermal)*dt; %PF*
end


%[T,Y] = ode45(@fun,[0,t_final],[N_U235,N_U238,0,0,0,0]);

figure;
semilogy(T,Y(:,1));
hold on;
semilogy(T,Y(:,2));
hold on;
semilogy(T,Y(:,3));
hold on;
semilogy(T,Y(:,4));
hold on;
semilogy(T,Y(:,5));
hold on;
semilogy(T,Y(:,6));
hold on;
xlabel('Temps [s]');
ylabel('Espèces [mol]');
legend('U235','U238','U239','Np239','Pu239','PF*');
hold off;

toc
end