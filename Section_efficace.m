function [sigma] = Section_efficace(X,Transfo,n_eV,Path)
%Source : https://www-nds.iaea.org/exfor/endf.htm
%sigma en [barn]

if strcmp('Fission',Transfo)
    if strcmp('U235',X)
    fileID = fopen(fullfile(Path,'U235_Fission.txt'));
    elseif strcmp('U238',X)
    fileID = fopen(fullfile(Path,'U238_Fission.txt'));
    elseif strcmp('U239',X)
    fileID = fopen(fullfile(Path,'U239_Fission.txt'));
    elseif strcmp('Np239',X)
    fileID = fopen(fullfile(Path,'Np239_Fission.txt'));
    elseif strcmp('Pu239',X)
    fileID = fopen(fullfile(Path,'Pu239_Fission.txt'));
    end
elseif strcmp('Capture',Transfo)
    if strcmp('U235',X)
    fileID = fopen(fullfile(Path,'U235_Capture.txt'));
    elseif strcmp('U238',X)
    fileID = fopen(fullfile(Path,'U238_Capture.txt'));
    elseif strcmp('U239',X)
    fileID = fopen(fullfile(Path,'U239_Capture.txt'));
    elseif strcmp('Np239',X)
    fileID = fopen(fullfile(Path,'Np239_Capture.txt'));
    elseif strcmp('Pu239',X)
    fileID = fopen(fullfile(Path,'Pu239_Capture.txt'));
    end
end

data = fscanf(fileID,'%f',[2 inf]);
data = data';

% loglog(data(:,1),data(:,2));
% title('Section efficace');
% xlabel('E [eV]');
% ylabel('Sigma [barn]');

i = 1;
while data(i,1) < n_eV
    i = i+1;
end
if n_eV == data(i,1)
    sigma = data(i,2);
else
    sigma = (n_eV-data(i-1,1))/(data(i,1)-data(i-1,1))*(data(i,2)-data(i-1,2)) + data(i-1,2);
end

fclose(fileID);


end