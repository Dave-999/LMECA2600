function [demi_vie] = Demi_vie(X,Decay)
%FROM : http://wwwndc.jaea.go.jp/NuC/
%Plus precisement : http://wwwndc.jaea.go.jp/CN14/index.html
%Input
%   [X] : 'U235' , 'U238' , 'U239', 'Np239', 'Pu239' or 'Xe135'
%   [Transfo] : Alpha, BetaMinus or BetaPlus
%Output
%   [demi_vie] : half-life expressed in seconds

if strcmp('U235',X)
    if strcmp('Alpha',Decay)
        demi_vie = 2.21950368*10^16;
    elseif strcmp('Beta-',Decay)
        demi_vie = 10^99;
        fprintf('\n The chemical %s is not implemented, demi_vie = 10^99 <=> never happend',X);
    elseif strcmp('Beta+',Decay)
        demi_vie = 10^99;
        fprintf('\n The chemical %s is not implemented, demi_vie = 10^99 <=> never happend',X);
    end
elseif strcmp('U238',X)
    if strcmp('Alpha',Decay)
        demi_vie = 1.40902848*10^17;
    elseif strcmp('Beta-',Decay)
        demi_vie = 10^99;
        fprintf('\n The chemical %s is not implemented, demi_vie = 10^99 <=> never happend',X);
    elseif strcmp('Beta+',Decay)
        demi_vie = 10^99;
        fprintf('\n The chemical %s is not implemented, demi_vie = 10^99 <=> never happend',X);
    end
elseif strcmp('U239',X)
    if strcmp('Alpha',Decay)
        demi_vie = 10^99;
        fprintf('\n The chemical %s is not implemented, demi_vie = 10^99 <=> never happend',X);
    elseif strcmp('Beta-',Decay)
        demi_vie = 1425;
    elseif strcmp('Beta+',Decay)
        demi_vie = 10^99;
        fprintf('\n The chemical %s is not implemented, demi_vie = 10^99 <=> never happend',X);
    end
elseif strcmp('Np239',X)
    if strcmp('Alpha',Decay)
        demi_vie = 10^99;
        fprintf('\n The chemical %s is not implemented, demi_vie = 10^99 <=> never happend',X);
    elseif strcmp('Beta-',Decay)
        demi_vie = 203558.4;
    elseif strcmp('Beta+',Decay)
        demi_vie = 10^99;
        fprintf('\n The chemical %s is not implemented, demi_vie = 10^99 <=> never happend',X);
    end
elseif strcmp('Pu239',X)
    if strcmp('Alpha',Decay)
        demi_vie = 7.6033296*10^11;
    elseif strcmp('Beta-',Decay)
        demi_vie = 10^99;
        fprintf('\n The chemical %s is not implemented, demi_vie = 10^99 <=> never happend',X);
    elseif strcmp('Beta+',Decay)
        demi_vie = 10^99;
        fprintf('\n The chemical %s is not implemented, demi_vie = 10^99 <=> never happend',X);
    end
end

end