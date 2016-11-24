function m = Reactor_model(t_final,dt_plot,P_stable,PF_retarded)

m_Utot = 1000; %[kg]
m_U235 = m_Utot*0.07;
m_U238 = m_Utot*0.93;

n_initial = 10^10;

%--------------------------------------------------

Section_efficace('U235','Fission',30000,'DATABASE')

end