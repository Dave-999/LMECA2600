function [sigma]=Section_efficace(X,Transfo,n_eV,Path)
% Section_efficace [barn]  Section efficace d'une transformation pour 1 composant
%   [SIGMA] = Section_efficace(X,TRANSFO,N_EV) donne la section 
%   efficace SIGMA de la transformation TRANSFO de l'element chimique X
%   lorsque le neutron incident a une énergie de N_EV
%
%   X : Espèce chimique.
%       ATTENTION : ici on travail avec un nombre limité d'espèces
%       chimiques
%
%   TRANSFO : seules les transformations ci dessous sont utilisées :
%       Fission : Probabilite qu'un noyau absorbe un neutron et fissionne
%       Capture : Probabilite qu'un noyau absorbe un neutron et fissionne 
%
%   n_eV : energie du neutron incident. Peut etre un vecteur
%       ATTENTION : On suppose une energie comprise entre 1e-5 et 2e7 [eV]
%
%   Path : Adress of the data base
%
%   ETAPES : etapes intermediaires pour construire la data base :
%       1°) Construction d'une data base pour les elements Ux et Np9 et Pu9
%       où Ux peuvent soit fissionner soit capturer jusqu'à U9. Puis U9 à
%       Np9 peuvent soit fissionner soit beta -. Pu peut juste fissionner
%
%   SOURCES : fichiers viennent de
%   https://www-nds.iaea.org/exfor/endf.htm
%   C'est la base ENDF qui a généralement été utilisé

%% Check argument IN : 

%On vérifie que l'espece chimique est dans nos database

if nargin==0
    n_eV=logspace(-5,6,10000);
    Path='C:/Users/ZINNIA/Desktop/GL/Cours/MECA2600-Génie des réacteurs nucléaires/Projet 2016-2017/Code/DATABASE';
    X='U235';
    Transfo='Fission';
elseif(n_eV~=sort(n_eV))
    error('The asked cross sections vector must be sorted!');
end

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
        sigma=0;
        fprintf('\n WARNING : There is no database for element %s. \n Please check function information',X);
        return;
end

%On vérifie que le type de transformation existe :
switch Transfo
    case 'Capture'
    case 'Fission'
    otherwise
        sigma=0;
        fprintf('\n WARNING : There is no database for Transformation %s. \n Please check function information',Transfo);
        return;
end

if (min(n_eV)<0.99e-5||max(n_eV)>2e7)
    error('\n Neutron energy of %d in data base of %s is out of range : [1e-5 2e7]',n_eV,X);
end

%% Reading information in file

switch Transfo
    case 'Capture'
        file_name=sprintf('%s/%s_CAPTURE.txt',Path,X);
        [E,Sig,~]=textread(file_name,'%f%f%s');
        
    case 'Fission'
        file_name=sprintf('%s/%s_FISSION.txt',Path,X);
        [E,Sig,~]=textread(file_name,'%f%f%s');        

    otherwise
        error('Trou dans la raquette :/');
end


counter=1;
sigma=0*n_eV;%Allocation taille
for ii=1:length(n_eV)% Il existe des abscisses à plusiuuers ordonnees :/
    for jj=counter:(length(E))
        if(((n_eV(ii)-E(jj))/n_eV(ii))<10^-6)% Différence relative négligeable. Ceci a cause d'erreurs d'arrondis :/
            sigma(ii)=Sig(jj);
            counter=jj;
            break;
        end
        if(E(jj)>n_eV(ii))
            sigma(ii)= (n_eV(ii)-E(jj-1))/(E(jj)-E(jj-1))*Sig(jj-1)+...
                       (E(jj)-n_eV(ii))/(E(jj)-E(jj-1))*Sig(jj); 
            counter=jj;
            break;
        end
        
    end
end

if nargin==0 
    loglog(n_eV,sigma,'k') 
end



end