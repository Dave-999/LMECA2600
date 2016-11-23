function M = molarMass(X)
% M = molarMass(X) returns the molar mass of element X
%
% INPUT : 
% [X] : 'U235' , 'U238' , 'U239', 'Np239', 'Pu239', 'n'
%
% OUTPUT : 
% M : Molar mass [kg/mol]
% FROM : http://wwwndc.jaea.go.jp/NuC/

if nargin == 0
    X='U235';
end

switch X
    case 'U235'
        M=235.0439299/1000;%kg/mol
    case 'U238'
        M=238.0507826/1000;%kg/mol
    case 'U239'
        M=239.054294518/1000;%kg/mol
    case 'Np239'
        M=239.052940487/1000;%kg/mol
    case 'Pu239'
        M=239.052164844/1000;%kg/mol
%     case 'PFstar' %assimilé à du Sn116 pour un souci de cohérence energétique (voir exo 2.4)
%         M=116.902956/1000;%kg/mol
%     case 'PF' %assimilé à du Pd115
%         M=114.913658506/1000;%kg/mol
    case 'n'
        M=1.00866491578/1000;%kg/mol
        
    otherwise
        fprintf('\n ---------------- WARNING ----------------- \n No molar mass for species (%s)',X);
        M=0;
end



end