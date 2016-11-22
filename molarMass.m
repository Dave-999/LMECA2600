function [M] = molarMass(X)
%Source : http://wwwndc.jaea.go.jp/NuC/
%masse molaire en [kg/kmole]

if strcmp('U235',X)
    M = 235.043931368;
elseif strcmp('U238',X)
    M = 238.050789466;
elseif strcmp('U239',X)
    M = 239.054294518;
elseif strcmp('Np239',X)
    M = 239.052940487;
elseif strcmp('Pu239',X)
    M = 239.052164844;
end

end