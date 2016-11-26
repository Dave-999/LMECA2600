function demi_vie = Demi_vie(X,Transfo)
%FROM : http://wwwndc.jaea.go.jp/NuC/
% Plus précisement : http://wwwndc.jaea.go.jp/CN14/index.html
% Input
%   [X] : 'U235' , 'U238' , 'U239', 'Np239', 'Pu239' or 'Xe135'
%   [Transfo] : Alpha, BetaMinus or BetaPlus
% Output
%   [demi_vie] : half-life expressed in seconds
%
% WARNING : version 1.0_2016/10/20. Only BetaMinus has been implemented for
% U239, Np239 and Xe135

%% Check argument IN : 
%On vérifie que l'espece chimique est dans nos database

switch X
    case 'U235'
    case 'U236'
    case 'U237'
    case 'U238'
    case 'U239'
    case 'Np239'
    case 'Pu239'
    case 'Xe135'
    otherwise
        fprintf('\n WARNING : There is no database for element %s. \n Please check function information',X);
        return;
end

%On vérifie que le type de transformation existe :
switch Transfo
    case 'Alpha'
    case 'BetaMinus'
    case 'BetaPlus'
    case 'Gamma' % Se révèle très courte
    otherwise
        fprintf('\n WARNING : These transformation are not implemented %s. \n Please check function information',Transfo);
end


%% Reading information in file
switch Transfo
    case 'Gamma'% A IMPLEMENTER
        
    case 'Alpha'
        switch X
            case 'U235'
                demi_vie = 2.21950368*10^16; %[s]
            case 'U238'
                demi_vie = 1.40902848*10^17; %[s]
            case 'Pu239'
                demi_vie = 7.6033296*10^11; %[s]
            otherwise
                fprintf('\n The chemical %s is not implemented. half_life=10^99 <=> never happend',X);
                demi_vie=10^99;
        end
        
    case 'BetaMinus'
        switch X
            case 'U239'
                demi_vie=23.45*60;%23.45min [s]
            case 'Np239'
                demi_vie=2.3565*24*3600;%2.3565 j [s]
            case 'Xe135'
                demi_vie=9.17*3600;%9.17 h [s]
            otherwise
                fprintf('\n The chemical %s is not implemented. half_life=10^99 <=> never happend',X);
                demi_vie=10^99;
        end
        
        
    case 'BetaPlus'% A IMPLEMENTER
        
        
    otherwise 
        fprintf('\n The chemical %s is not implemented. half_life=10^99 <=> never happend',X);
        demi_vie=10^99;
end


end