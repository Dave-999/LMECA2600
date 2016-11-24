function m = Reactor_model(t_final,dt_plot,P_stable,PF_retarded)
%Donnees

V = 30; %Volume du reacteur en [m^3]
m_Utot = 1; %Masse d'Uranium [kg]

m_U235 = m_Utot*0.07; %Masse d'U235 [kg]
N_U235 = m_U235/molarMass('U235'); %Nombre de moles d'U235 [mol]

m_U238 = m_Utot*0.93; %Masse d'U238 [kg]
N_U238 = m_U238/molarMass('U238'); %Nombre de moles d'U238 [mol]

n_thermal = 10^10; %Nombre de neutrons thermiques en t=0
flux_thermal = n_thermal/V; %Flux de neutrons thermiques en t=0

E_thermal = 0.025; %eV
E_fast = 1e6; %eV

%--------------------------------------------------------------------------

t_final = 10;

[T,Y] = ode45(@f1,[0,t_final],N_U235);

figure;
plot(T,N_U235-Y);
xlabel('t [s]');
ylabel('Species [mol]');

end

function dydt = f1(t,y)
    V = 30;
    n_thermal = 10^10;
    E_thermal = 0.025;
    flux_thermal = n_thermal/V;
    dydt = y*flux_thermal*10e-28*Section_efficace('U235','Fission',E_thermal,'DATABASE');
end