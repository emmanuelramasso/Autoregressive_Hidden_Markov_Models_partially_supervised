function perf = testARPHMM_salves_silmulees(lessalves, labelInit) 

addpath PHMM/
addpath PHMM/utils/
addpath utils/

S={}; for i=1:length(lessalves), for j=1:size(lessalves{i},1), 
        S=[S; (lessalves{i}(j,:)')/max(lessalves{i}(j,:))]; 
    end, 
end
salvesOut = S;% lesSalves 
% for i=1:length(salvesOut)
%     salvesOut{i}=salvesOut{i}+randn(size(salvesOut{i}))*0.1;
% end

nbAvant = 200000
nbApres = 100000 % MC1 LOC

for i=1:length(salvesOut)
    [a b]=max(salvesOut{i});
    salvesOut{i} = salvesOut{i}(max(1,b-nbAvant+1):min(b+nbApres,length(salvesOut{i})));
    salvesOut{i} = salvesOut{i}(:);
end

% figure, for i=1:length(salvesOut), clf, plot(salvesOut{i}), i, pause(1), end

% numero des salves choisies pour apprendre les modeles
% referencesSalves=[26 60 130] % MC2 LOC
referencesSalves=[1]% 107 205 306]

salvesModeles = salvesOut(referencesSalves);
% salvesModeles{1}=salvesModeles{1}(345:704);
for i=1:length(referencesSalves)
    figure,plot(salvesModeles{i})
    title(sprintf('Salve de référence %d (%d-ieme de la liste)',i,referencesSalves(i)))
end


clear models
k = 1 % ne pas toucher
predc = 10; % nb de predecesseurs dans le AR
Q = 5; % nb d'etats dans la chaine de Markov
D = 5 % pas important, pour affichage

% fichier temporaire pour l'apprentissage
nomTemp = 'tmptmtptmp';
s=num2str(clock);
f=strfind(s,' '); s(f)=[]; 
f=strfind(s,'.'); s(f)=[];
nomTemp = [nomTemp s];

for lsa=1:length(salvesModeles)% parcours les salves de référence et apprend un modele pour chaque salve
    
    z=salvesModeles{lsa};%(1:400)';% un autre
    
    figure,plot(z)
    
    y = z;
    my = mean(y); stdy = std(y);
    y=(y-my)/stdy;
    
    % particularité, à voir plus tard
    if 0
        gtruth=randi(Q,[length(y),1]);
        pl0=labels2matrix(gtruth, Q);
    else
        pl0 = ones(length(y),Q);
    end
    
    figure, hold on,
    qualitymf_supZ = -inf;
    for l=1:3
        [gamma2_sup, NDEI_sup, qualitymf_sup, signalt_sup, ...
            signalcomb_sup, B_sup, sigma_sup, A_sup, Pi_sup, ~,~, ~, ~, ~, LL_sup] =...
            arphmm(y, pl0,Q,predc);
        
        stem(y), plot(signalcomb_sup,'ro')
        
        if LL_sup(end) > qualitymf_supZ %qualitymf_sup<qualitymf_supZ
            save(nomTemp,'LL_sup','gamma2_sup','NDEI_sup','qualitymf_sup','signalt_sup','signalcomb_sup','B_sup','sigma_sup','A_sup','Pi_sup');
            qualitymf_supZ = LL_sup(end);%qualitymf_sup;
        end
    end
    clear('LL_sup','gamma2_sup','NDEI_sup','qualitymf_sup','signalt_sup','signalcomb_sup','B_sup','sigma_sup','A_sup','Pi_sup')
    load(nomTemp)
    
    mean((y(D:end)-signalcomb_sup(D:end)).^2)
    [a b]=max(gamma2_sup,[],2);
    %[AR,~,~,~, C]=valid_RandIndex(gtruth, b)
    
    figure, hold on, stem(y), plot(signalcomb_sup,'ro')
    %figure,plot(signalcomb_sup-y,'ro')
    
    models(k).B = B_sup;
    models(k).S = sigma_sup;
    models(k).A = A_sup;
    models(k).Pi = Pi_sup;
    models(k).gamma = gamma2_sup;
    models(k).pl = pl0;
    models(k).predc = predc;
    models(k).Q = Q;
    models(k).T = length(y);
    models(k).my = my;
    models(k).stdy = stdy;
    models(k).NDEI = NDEI_sup;
    models(k).qualitymf = qualitymf_sup;
    models(k).LL = LL_sup;
    
    k=k+1;
    close all
end

if 0
    Fe=200;
    close all; x=salvesModeles{lsa}; figure(2);
    [P,w]=periodogram(x,[],[],Fe); plot(w,10*log10(P)); s={'period. x'}; hold on, for i=1:size(B_sup,1)
        ff=filter(B_sup(i,:),1,x);[P,w]=periodogram(ff,[],[],Fe); plot(w,10*log10(P)); s{end+1}=sprintf('period %d',i);end,
    xlabel('Frequency (Hz)');
    ylabel('One-sided PSD (dB/Hz)');
    hold on, [P,w]=periodogram(signalcomb_sup,[],[],Fe); plot(w,10*log10(P)); s{end+1}='period. final'; legend(s);
end

% Fe=200;
% close all; x=salvesModeles{lsa}; 
% for i=1:size(B_sup,1)
%     figure, [P,w]=periodogram(x,[],[],Fe); plot(w,10*log10(P)); s={'period. x'}; hold on, 
% ff=filter(B_sup(i,:),1,x);[P,w]=periodogram(ff,[],[],Fe); plot(w,10*log10(P)); s{end+1}=sprintf('period %d',i); legend(s), end, 
% xlabel('Frequency (Hz)');
% ylabel('One-sided PSD (dB/Hz)');

% figure,hold on, for i=1:size(B_sup,1), impulse(tf(1,B_sup(i,:),1/Fe)), end


% ici on applique les modeles
lesreponsesLOGL = zeros(size(salvesOut,2),length(models));
lesreponsesMSE = zeros(size(salvesOut,2),length(models));
lesreponsesNDEI = zeros(size(salvesOut,2),length(models));
lesreponsesMAXERR = zeros(size(salvesOut,2),length(models));
lesreponsesPSNR = zeros(size(salvesOut,2),length(models));
lesreponsesL2RAT = zeros(size(salvesOut,2),length(models));

for lasalve=1:1:length(salvesOut)
    
    z = salvesOut{lasalve}(:);
    
    for lemodel=1:length(models)
        
        try % cas ou la taille de la salve est trop courte wrt apprentissage
            
            y = z;
            
            %%% INFERE
            y = y - models(lemodel).my;
            y = y / models(lemodel).stdy;
            
            T =length(y);
            predata=cell(T,1);
            datmat=zeros(T+models(lemodel).predc,1);
            datmat(1:models(lemodel).predc,:)=y(1:models(lemodel).predc,:);
            datmat(models(lemodel).predc+1:T+models(lemodel).predc,:)=y;
            for t=1:T
                for pp=1:models(lemodel).predc
                    predata{t}(pp,:)=datmat(t+models(lemodel).predc-pp,:);
                end
            end
            
            % Calcul of b
            obslik=zeros(T,Q);
            for i=1:Q
                for t=1:T
                    obslik(t,i)=mvnpdf(y(t,:)+(models(lemodel).B(i,:)*predata{t}),0,models(lemodel).S{i});
                end
            end
            
            if isempty(find(sum(obslik,2)==0,1))
                
                if 1
                    [~,~, gamma, loglik] = fwdback_phmm_mix(models(lemodel).Pi, ...
                        models(lemodel).A, obslik, ones(T,models(lemodel).Q));
                else
                    %                     [~,~, gamma, loglik] = fwdback_phmm_mix(models(lemodel).Pi, ...
                    %                         models(lemodel).A, obslik, models(lemodel).gamma);
                end
                
                signal=autoregres(predata,models(lemodel).B,models(lemodel).S);
                ssmat=cell2mat(signal);
                signalt=reshape(ssmat,T,Q);
                
                % le signal de sortie est donné par signal{k} avec le max de gamma
                signalcomb=zeros(T,1);
                for t=1:T
                    for j=1:Q
                        signalcomb(t)=signalcomb(t) + gamma(t,j)*signalt(t,j);
                    end
                end
                %figure,plot(signalcomb,'ro'), hold on, stem(y)
                
                [PSNR,MSE,MAXERR,L2RAT] = measerr(y,signalcomb);
                
                lesreponsesLOGL(lasalve,lemodel) = loglik / T;
                lesreponsesMSE(lasalve,lemodel) = MSE; %mean((signalcomb-y).^2);
                lesreponsesNDEI(lasalve,lemodel) = sqrt(lesreponsesMSE(lasalve,lemodel))/std(y);
                lesreponsesMAXERR(lasalve,lemodel) = PSNR;
                lesreponsesPSNR(lasalve,lemodel) = MAXERR;
                lesreponsesL2RAT(lasalve,lemodel) = L2RAT;
                
            else
                lesreponsesLOGL(lasalve,lemodel) = -inf;
                lesreponsesMSE(lasalve,lemodel) = inf;
                lesreponsesNDEI(lasalve,lemodel) = inf;
                lesreponsesMAXERR(lasalve,lemodel) = inf; 
                lesreponsesPSNR(lasalve,lemodel) = inf;
                lesreponsesL2RAT(lasalve,lemodel) = inf;
            end
        catch
            lesreponsesLOGL(lasalve,lemodel) = -inf;
            lesreponsesMSE(lasalve,lemodel) = inf;
            lesreponsesNDEI(lasalve,lemodel) = inf;
            lesreponsesMAXERR(lasalve,lemodel) = inf;
            lesreponsesPSNR(lasalve,lemodel) = inf;
            lesreponsesL2RAT(lasalve,lemodel) = inf;
        end
        [lasalve, length(salvesOut), lemodel, length(models)]
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % analyse des resultats
% on veut segmenter l'indicateur


% [A B]=max(lesreponsesLOGL,[],2); figure,plot(B)
% [A B]=min(lesreponsesMSE,[],2);figure,plot(B)
% [A B]=min(lesreponsesNDEI,[],2); figure,plot(B)
% figure,imagesc(lesreponsesLOGL)
% figure,imagesc(lesreponsesMSE)
% figure,imagesc(lesreponsesNDEI)

% data
noms2={'LOGL','MSE','NDEI'};
% global noms, noms={'LOGL','NDEI'};
X=[lesreponsesLOGL(:,:) lesreponsesNDEI(:,:)];%lesreponsesMSE(:,:) 

% ############### POINT IMPORTANT
% on ramene tt a la reference
% par exemple la premiere reference
% lsa = referencesSalves(1)
% 
% ref1=X(lsa,1); X(:,1)=X(:,1)./ref1;
% ref2=X(lsa,2); X(:,2)=(X(:,2)-ref2)./ref2;
% ref3=X(lsa,3); X(:,3)=(X(:,3)-ref3)./ref3;
% 
% X(lsa,:) % 1ere colonne => =1 ??
% figure,plot(X), legend(noms2)


% X=[lesreponsesLOGL(:,1) lesreponsesNDEI(:,1)];
aretirer = find(isinf(sum(X,2)) | isnan(sum(X,2)) );
garde = 1:size(X,1); garde(aretirer)=[];
disp(sprintf('Retire %d pts / %d',length(aretirer),size(X,1)))

X(aretirer,:)=[];
labelInit2=labelInit; labelInit2(aretirer)=[];


KKK=4;

% que la likelihood
% global normaliser, normaliser=true;
% c1 = GK_multi_init(X, 1, KKK, 1, [] ,[], 'gmm1');%, 'fcm');
% [c1 c2]=kmeans(zscore(X), KKK, 'replicates', 100);

% V1
% [c1 c2]=kmeans((X), KKK, 'replicates', 100);
% tabulate(c1)
% [perf.a perf.b perf.c perf.d perf.r]=valid_RandIndex(labelInit2,c1)

% V2
options = statset('MaxIter',1000)
Sigma = 'full';%'diagonal';%,'full'};
SharedCovariance = true; %{true,false};
gmfit = fitgmdist((X),KKK,'CovarianceType',Sigma,...
'SharedCovariance',SharedCovariance,'Options',options);
c1 = cluster(gmfit,X);
tabulate(c1)
[perf.a perf.b perf.c perf.d perf.r]=valid_RandIndex(labelInit2,c1)
perf.r


% peut on trouver automatiquement le nb de clusters ?

clear eva
k=1;
rng('default');  % For reproducibility
methods={'kmeans','gmdistribution','linkage'};
indices={'CalinskiHarabasz','silhouette','DaviesBouldin'};%,'gap'};% gap peut etre long!!!!
for i=1:length(methods)
    for j=1:length(indices)
        disp(sprintf('En cours : method %s - indic. %s',methods{i},indices{j}))
        eva{k} = evalclusters(zscore(X),methods{i},indices{j},'Klist',[3:8]);
        k=k+1;
    end
end


perf.eva = eva;

