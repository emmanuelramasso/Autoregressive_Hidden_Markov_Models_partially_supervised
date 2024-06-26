function [lesreponsesLOGL, lesreponsesMSE, ...
    lesreponsesNDEI, lesreponsesMAXERR, ...
    lesreponsesPSNR, lesreponsesL2RAT] = ...
    appliquerModelesSurSalves(models, salvesOut)

% ici on applique les modeles
lesreponsesLOGL = zeros(size(salvesOut,2),length(models));
lesreponsesMSE = zeros(size(salvesOut,2),length(models));
lesreponsesNDEI = zeros(size(salvesOut,2),length(models));
lesreponsesMAXERR = zeros(size(salvesOut,2),length(models));
lesreponsesPSNR = zeros(size(salvesOut,2),length(models));
lesreponsesL2RAT = zeros(size(salvesOut,2),length(models));

disp(sprintf('Process %d salves...',length(salvesOut)))

k=0;
for lasalve=1:length(salvesOut)

    k=k+1;
    if mod(k,length(salvesOut)*5/100)==0
        disp(sprintf('%d/%d',lasalve,length(salvesOut)));
        k=0;
    end
        
    [lesreponsesLOGL(lasalve,:),...
    lesreponsesMSE(lasalve,:), ...
    lesreponsesNDEI(lasalve,:), ...
    lesreponsesMAXERR(lasalve,:), ...
    lesreponsesPSNR(lasalve,:), ...
    lesreponsesL2RAT(lasalve,:)] = appliquerModeles(models, salvesOut{lasalve}(:));
    
end
disp('Inference sur salves ok.')


