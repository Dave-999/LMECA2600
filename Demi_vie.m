function [demi_vie] = Demi_vie(X,Decay)
%Source : http://wwwndc.jaea.go.jp/NuC/
%demi-vie en [seconde]

if strcmp('U235',X)
    if strcmp('Alpha',Decay)
        demi_vie = 2.21950368*10^16;
    elseif strcmp('Beta-',Decay)
        demi_vie = inf;
        fprintf('Erreur : pas de %s pour %s, demi-vie "infinie"', Decay, X);
    elseif strcmp('Beta+',Decay)
        demi_vie = inf;
        fprintf('Erreur : pas de %s pour %s, demi-vie "infinie"', Decay, X);
    end
elseif strcmp('U238',X)
    if strcmp('Alpha',Decay)
        demi_vie = 1.40902848*10^17;
    elseif strcmp('Beta-',Decay)
        demi_vie = inf;
        fprintf('Erreur : pas de %s pour %s, demi-vie "infinie"', Decay, X);
    elseif strcmp('Beta+',Decay)
        demi_vie = inf;
        fprintf('Erreur : pas de %s pour %s, demi-vie "infinie"', Decay, X);
    end
elseif strcmp('U239',X)
    if strcmp('Alpha',Decay)
        demi_vie = inf;
        fprintf('Erreur : pas de %s pour %s, demi-vie "infinie"', Decay, X);
    elseif strcmp('Beta-',Decay)
        demi_vie = 1425;
    elseif strcmp('Beta+',Decay)
        demi_vie = inf;
        fprintf('Erreur : pas de %s pour %s, demi-vie "infinie"', Decay, X);
    end
elseif strcmp('Np239',X)
    if strcmp('Alpha',Decay)
        demi_vie = inf;
        fprintf('Erreur : pas de %s pour %s, demi-vie "infinie"', Decay, X);
    elseif strcmp('Beta-',Decay)
        demi_vie = 203558.4;
    elseif strcmp('Beta+',Decay)
        demi_vie = inf;
        fprintf('Erreur : pas de %s pour %s, demi-vie "infinie"', Decay, X);
    end
elseif strcmp('Pu239',X)
    if strcmp('Alpha',Decay)
        demi_vie = 7.6033296*10^11;
    elseif strcmp('Beta-',Decay)
        demi_vie = inf;
        fprintf('Erreur : pas de %s pour %s, demi-vie "infinie"', Decay, X);
    elseif strcmp('Beta+',Decay)
        demi_vie = inf;
        fprintf('Erreur : pas de %s pour %s, demi-vie "infinie"', Decay, X);
    end
end

end