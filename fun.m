function dydt = fun(t,y)

    E_thermal = 0.025;
    n_thermal = 10^10;
    v_thermal = 10;
    V = 30;
    flux_thermal = n_thermal*v_thermal/V;
    
    fis_U235 = Section_efficace('U235','Fission',E_thermal,'DATABASE');
    cap_U235 = Section_efficace('U235','Capture',E_thermal,'DATABASE');
    demi_U235 = Demi_vie('U235','Alpha');
    
    cap_U238 = Section_efficace('U238','Capture',E_thermal,'DATABASE');
    demi_U238 = Demi_vie('U238','Alpha');
    
    fis_U239 = Section_efficace('U239','Fission',E_thermal,'DATABASE');
    demi_U239 = Demi_vie('U239','BetaMinus');
    
    fis_Np239 = Section_efficace('Np239','Fission',E_thermal,'DATABASE');
    demi_Np239 = Demi_vie('Np239','BetaMinus');
    
    fis_Pu239 = Section_efficace('Pu239','Fission',E_thermal,'DATABASE');
    demi_Pu239 = Demi_vie('Pu239','Alpha');
    
    %----------------------------------------------------------------------
    
    dydt = zeros(6,1);
    dydt(1) = - y(1)*fis_U235*10e-28*flux_thermal - y(1)*cap_U235*10e-28*flux_thermal - y(1)*log(2)/demi_U235; %U235
    dydt(2) = - y(2)*cap_U238*10e-28*flux_thermal - y(2)*log(2)/demi_U238; %U238
    dydt(3) = y(2)*cap_U238*10e-28*flux_thermal - y(3)*fis_U239*10e-28*flux_thermal - y(3)*log(2)/demi_U239; %U239
    dydt(4) = y(3)*log(2)/demi_U239 - y(4)*fis_Np239*10e-28*flux_thermal - y(4)*log(2)/demi_Np239; %Np239
    dydt(5) =  y(4)*log(2)/demi_Np239 - y(5)*fis_Pu239*10e-28*flux_thermal - y(5)*log(2)/demi_Pu239; %Pu239
    dydt(6) = y(1)*fis_U235*10e-28*flux_thermal + y(3)*fis_U239*10e-28*flux_thermal + y(4)*fis_Np239*10e-28*flux_thermal + y(5)*fis_Pu239*10e-28*flux_thermal; %PF*
    
end

