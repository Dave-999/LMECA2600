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
flux_thermal = n_thermal*v_thermal/V; %Flux de neutrons thermiques en t=0

%--------------------------------------------------------------------------
%%

t_final = 100; %[s]

[T,Y] = ode45(@fun,[0,t_final],[N_U235,N_U238,0,0,0,0]);

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
loglog(T,Y(:,6));
hold on;
xlabel('Temps [s]');
ylabel('Espèces [mol]');
legend('U235','U238','U239','Np239','Pu239','PF*');
hold off;

toc
end