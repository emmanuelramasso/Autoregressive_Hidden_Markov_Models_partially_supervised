function perf = testARPHMM_salves_silmulees(lessalves, labelInit) 
% ADD TO PATH WITH SUBFOLDERS
% /home/emmanuel.ramasso/Documents/PROGRAMMES/GITHUB/CHMM_with_partial_labels 
% AND CURRENT FOLDER
%
% Run [lessalves, lesentrees, labelInit] = test_generation 
% to generate signals in four classes
% Then run perf = testARPHMM_salves_silmulees(lessalves, labelInit)
% to evaluate the performance in terms of different criteria
%
% The algorithm works as follows, see [3]:
% Phase 1 : train on a few signals (in vector "reference", hard coded here)
% Phase 2 : apply the model on testing, generate anomaly values and cluster
% them.
%
% References:
% The idea of including partial knowledge on states is from [1] and [2]
% applied to prognostics, while its application to acoustic emission
% signals is from [3]. 
%
% [1] On partially supervised learning and inference in dynamic Bayesian 
% networks for prognostics with uncertain factual evidence: Illustration 
% with Markov switching models, Pablo Juesas, Emmanuel Ramasso, 
% Sébastien Drujont, Vincent Placet,  Proceedings of the European 
% Conference of the PHM Society, Vol. 3 No. 1 (2016) 
%
% [2] Autoregressive Hidden Markov Models with partial knowledge on latent, 
% space applied to aero-engines prognostics, Pablo Juesas, Emmanuel
% Ramasso, Sébastien Drujont, Vincent Placet,  arxiv https://arxiv.org/abs/2105.00211, 2024.
% 
% [3] Ramasso, E., Butaud, P., Jeannin, T., Sarasini, F., Placet, V., 
% Godin, N., ...Gabrion, X. (2020). Learning the representation of raw 
% acoustic emission signals by direct generative modelling and its use in
% chronology-based clusters identification. Eng. Appl. Artif. Intell., 90,
% 103478. doi: 10.1016/j.engappai.2020.103478
%
% E. Ramasso 2016-2024
%

 
% signals selected to train the model, can specify several of them
referencesSalves = [1] % 107 205 306

% normalise signals
S={}; 
for i=1:length(lessalves), 
    for j=1:size(lessalves{i},1), 
        S=[S; (lessalves{i}(j,:)')/max(lessalves{i}(j,:))]; 
    end
end
salvesOut = S;

% add noise eventually 
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

% Training 
for lsa=1:length(salvesModeles)% parcours les salves de référence et apprend un modele pour chaque salve
    
    z=salvesModeles{lsa};%(1:400)';% un autre
    
    figure,plot(z)
    
    y = z;
    my = mean(y); stdy = std(y);
    y=(y-my)/stdy;
    
    % particularité, à voir plus tard
    if 0 
        % avec a priori
        gtruth=randi(Q,[length(y),1]);
        pl0=labels2matrix(gtruth, Q);
    else 
        % sans a priori
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

% Inference
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

% ############### 
% test : on ramene tt a la reference
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

% nb clusters expected ?
KKK=4;

% test
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
disp('adjusted Rand index, unadjusted Rand index, Mirkin''s index and Hubert''s index')
[perf.ari perf.ri perf.mi perf.hi perf.c]=valid_RandIndex(labelInit2,c1);
perf, perf.c

% Find the nb of clusters automatically
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

