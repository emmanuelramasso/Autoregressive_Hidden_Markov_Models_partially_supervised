
clear all; close all
% load '/home/emmanuel.ramasso/Documents/Dropbox/6_BENCHMARK/Biocomposites/FLAX_2017/P1_UD_07_TS0_10/file_2017_04_13_08_47_27953_Features_et_salves_0.mat'
%load('/home/emmanuel.ramasso/Documents/Dropbox/6_BENCHMARK/NATHALIE/minicomposite/MC2_locs/caracMC2_loc_NG.mat','carac','lesSalves');
% load('/home/emmanuel.ramasso/Documents/Dropbox/6_BENCHMARK/NATHALIE/minicomposite/MC1_locs/caracMC1_loc_NG.mat','carac','lesSalves');
load('/home/emmanuel.ramasso/Documents/Dropbox/6_BENCHMARK/NATHALIE/minicomposite/MC1/caracMC1_NG.mat','carac','lesSalves');
% load('/home/emmanuel.ramasso/Documents/Hanaa/P3_4/P3_4_salves.mat');
% load('/home/emmanuel.ramasso/Documents/Hanaa/P3_4/P3_4.mat'); salves(1:13)=[]; carac(1:13,:)=[]; lesSalves=salves; clear salves

addpath PHMM/
addpath PHMM/utils/
addpath utils/

% [Features,tempsSalves,ordre] = reordonneTri(Features);
% salvesOut=salvesOut(ordre);
salvesOut = lesSalves; clear lesSalves
Features = carac;
tempsSalves = carac(:,5); clear carac


MiniComp1 = lireMistrasMiniComp1;
figure,plot(MiniComp1(:,2),MiniComp1(:,4)*2)
tempsCarac=Features(:,5); tempsCarac=tempsCarac-tempsCarac(1); tempsCarac=tempsCarac/max(tempsCarac);
tempsDef=MiniComp1(:,2); tempsDef=tempsDef-tempsDef(1); tempsDef=tempsDef/max(tempsDef);
figure,plot(tempsCarac), hold on, plot(tempsDef)
defcarac = zeros(size(tempsCarac));

for i=1:length(tempsCarac)
    r=tempsDef;%(max(1,i-100):min(i+100,length(tempsDef)));
    u=MiniComp1(:,4)*2;%(max(1,i-100):min(i+100,length(tempsDef)),4)*2;
    d=(tempsCarac(i)-r).^2;
    [a b]=min(d);
    defcarac(i)=u(b);
end
figure,plot(tempsCarac,defcarac)


disp('Normalisation des salves...')
S={};
for i=1:length(salvesOut),
    S=[S; salvesOut{i}/max(salvesOut{i})];
end
salvesOut = S;%

% selectionne les pts avant et apres le max pour simplifier
% nbAvant = 50; nbApres = 150; % MC2 LOC
nbAvant = 100; nbApres = 300; % MC1 LOC

for i=1:length(salvesOut)
    [a b]=max(salvesOut{i});
    salvesOut{i} = salvesOut{i}(max(1,b-nbAvant+1):min(b+nbApres,length(salvesOut{i})));
    salvesOut{i} = salvesOut{i}(:);
end

% figure, for i=1:length(salvesOut), clf, plot(salvesOut{i}), i, pause(1), end

% numero des salves choisies pour apprendre les modeles
% referencesSalves=[26 60 130] % MC2 LOC
referencesSalves=[22 60 150 1500 6589 8000] % P3_4
referencesSalves=[742 1914 2820 8000 8200 8548 9042] %MC1 LOC
% CLASSE 2
% à 8000e => 2158.4229277   136.7796     0.3295*2
% à 8200e => 2561.8036060   141.8959     0.3569*2
% à 8540e => 3017.3715425   151.9470     0.4083*2
% à 9042e => 3038.3982967   161.4683     0.4438*2
% CLASSE 1
% à 742e  => 491.5214350    84.6320     0.0504*2
% à 1915e => 689.7546858    99.9811     0.0805*2
% à 2820e => 705.4806000   104.3557     0.1082*2


salvesModeles = salvesOut(referencesSalves);
% salvesModeles{1}=salvesModeles{1}(345:704);
for i=1:length(referencesSalves)
    figure,plot(salvesModeles{i})
    title(sprintf('Salve de référence %d (%d-ieme de la liste)',i,referencesSalves(i)))
end


clear models
k=1
predc = 20;
Q = 3;
D = 5

for lsa=1:length(salvesModeles)%    48;
    
    z=salvesModeles{lsa};%(1:400)';% un autre
    
    figure,plot(z)
    
    y = z;
    my = mean(y); stdy = std(y);
    y=(y-my)/stdy;
    
    if 0
        gtruth=randi(Q,[length(y),1]);
        pl0=labels2matrix(gtruth, Q);
    else
        pl0 = ones(length(y),Q);
    end
    
    figure, hold on,
    qualitymf_supZ = inf;
    for l=1:3
        [gamma2_sup, NDEI_sup, qualitymf_sup, signalt_sup, ...
            signalcomb_sup, B_sup, sigma_sup, A_sup, Pi_sup, ~,~, ~, ~, ~, LL_sup] =...
            arphmm(y, pl0,Q,predc);
        
        stem(y), plot(signalcomb_sup,'ro')
        
        if qualitymf_sup<qualitymf_supZ
            save tmptmp5454354454 LL_sup gamma2_sup NDEI_sup qualitymf_sup signalt_sup signalcomb_sup B_sup sigma_sup A_sup Pi_sup
            qualitymf_supZ = qualitymf_sup;
        end
    end
    clear gamma2_sup NDEI_sup qualitymf_sup signalt_sup signalcomb_sup B_sup sigma_sup A_sup Pi_sup
    load tmptmp5454354454
    
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

% test la reconnaissance 
[lesreponsesLOGL_train, lesreponsesMSE_train, ...
    lesreponsesNDEI_train, lesreponsesMAXERR_train, ...
    lesreponsesPSNR_train, lesreponsesL2RAT_train] = ...
    appliquerModelesSurSalves(models, salvesOut(referencesSalves));


figure,imagesc(lesreponsesLOGL_train)
lesreponsesLOGL_train


tic
[lesreponsesLOGL, lesreponsesMSE, ...
    lesreponsesNDEI, lesreponsesMAXERR, ...
    lesreponsesPSNR, lesreponsesL2RAT] = ...
    appliquerModelesSurSalves(models, salvesOut);
toc

figure
subplot(321),hold on,plot(lesreponsesL2RAT,'o'),title('L2RAT')
subplot(322),hold on,plot(lesreponsesLOGL,'o'),title('LOGL')
subplot(323),hold on,plot(lesreponsesMAXERR,'o'),title('MAXE')
subplot(324),hold on,plot(lesreponsesNDEI,'o'),title('NDEI')
subplot(325),hold on,plot(lesreponsesPSNR,'o'),title('PSNR')
subplot(326),hold on,plot(lesreponsesMSE,'o'),title('MSE')

X=[lesreponsesL2RAT lesreponsesLOGL lesreponsesMAXERR lesreponsesNDEI lesreponsesPSNR lesreponsesMSE defcarac];

figure
subplot(321),hold on,plot(diff(X(:,1)),'o'),title('L2RAT')
subplot(322),hold on,plot(diff(X(:,2)),'o'),title('LOGL')
subplot(323),hold on,plot(diff(X(:,3)),'o'),title('MAXE')
subplot(324),hold on,plot(diff(X(:,4)),'o'),title('NDEI')
subplot(325),hold on,plot(diff(X(:,5)),'o'),title('PSNR')
subplot(326),hold on,plot(diff(X(:,6)),'o'),title('MSE')

figure
subplot(321),hold on,plot(diff(medfilt1(X(:,1),15,'omitnan')),'o'),title('L2RAT')
subplot(322),hold on,plot(diff(medfilt1(X(:,2),15,'omitnan')),'o'),title('LOGL')
subplot(323),hold on,plot(diff(medfilt1(X(:,3),15,'omitnan')),'o'),title('MAXE')
subplot(324),hold on,plot(diff(medfilt1(X(:,4),15,'omitnan')),'o'),title('NDEI')
subplot(325),hold on,plot(diff(medfilt1(X(:,5),15,'omitnan')),'o'),title('PSNR')
subplot(326),hold on,plot(diff(medfilt1(X(:,6),15,'omitnan')),'o'),title('MSE')

% figure
% for i=1:length(models)
%     ll=diff(lesreponsesL2RAT);%  lesreponsesLOGL; 
%     %ll=ll/abs(ll(referencesSalves(i),i));
%     %ll(find(ll>0))=-realmin;
%     subplot(4,2,i),
%     plot(tempsSalves(2:end),exp(ll(:,i))./repmat(sum(exp(ll),2),1,7))
%     title(sprintf('class %d',i))
% end


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
% noms2={'LOGL','MSE','NDEI'};
% global noms, noms={'LOGL','NDEI'};
% X=[lesreponsesLOGL(:,:) lesreponsesNDEI(:,:)];% lesreponsesMSE(:,:) lesreponsesNDEI(:,:)];

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

% system.c = 68000000;
% system.type = 'barre'; % anneau ou barre
% system.nbCapteurs = 2;
% system.pos = [65 145];
% [caracLocalisees, indices, positionLocalisation] = localisation_barre_generique(Features, system);


T=tempsSalves;


[a b]=kmeans(lesreponsesL2RAT,5);
global reordonnerPourVisu, reordonnerPourVisu=true;
cc=duree_vs_amplitude(a, lesreponsesL2RAT, T, KKK, 1,1);
close(gcf)
visualiser_resultats_shm_clust_morceaux_2(cc,charge(:),KKK,[]);
superpose_energie(Features(:,22),charge(:),1);




indic_L2RAT=indicateursARPHMM(lesreponsesL2RAT,'min');
indic_LOGL=indicateursARPHMM(lesreponsesLOGL,'max');
indic_MAXERR=indicateursARPHMM(lesreponsesMAXERR,'min');
indic_MSE=indicateursARPHMM(lesreponsesMSE,'min');
indic_NDEI=indicateursARPHMM(lesreponsesNDEI,'min');
indic_PSNR=indicateursARPHMM(lesreponsesPSNR,'min');

% X=[indic_L2RAT indic_LOGL indic_MAXERR indic_MSE indic_NDEI indic_PSNR];
X=[indic_L2RAT indic_LOGL indic_NDEI indic_MAXERR indic_PSNR indic_MSE];

[L1, L2, F1, F2, doute]=labelsMajoritaire(X);
X=[X L1];
charge = defcarac;
if 1
    m=0;
    for i=2:length(charge)
        m=max(m,charge(i-1));
        if charge(i)<m, charge(i)=m; end
    end
end
    
for i=1:size(X,2)
    
    garde=1:size(X,1);
    garde(find(X(:,i)==0))=[];
   
    
    rr=1:size(X,2);
    rr(i)=[]; if isempty(rr), rr=i; end
    
    KKK=length(unique(X(:,i)));
    
    % ON AFFICHE LES RESULTATS SANS CONSIDERE Bbarre
    global reordonnerPourVisu, reordonnerPourVisu=true;
    cc=duree_vs_amplitude(X(garde,i), X(garde,:), T(garde), KKK, i, rr(1));
    tabulate(cc)
    close(gcf)
    
    %visualiser_resultats_shm_clust_morceaux_2(cc,T(garde),KKK,[]);
    %superpose_energie(Features(garde,22),T(garde),1);
    
    visualiser_resultats_shm_clust_morceaux_2(cc,charge(garde),KKK,[]);
    superpose_energie(Features(garde,22),charge(garde),1);
    
    %densites_amplitude(cc, Features(garde,13), 1) %F(garde,5), 1);
    
    %tracer_energie_dans_cluster(cc, Features(garde,11), T(garde));
    
    %tracer_energie_dans_cluster_normalise(cc, Features(garde,11), T(garde));
    
    %densites_amplitude(cc, Features(garde,29), 1) %F(garde,5), 1);
    
    %duree_vs_amplitude(cc, [Features(garde,13) log10(Features(garde,12))], T(garde), KKK, 1,2);
    
    %duree_vs_amplitude(cc, [Features(garde,29) Features(garde,24)], T(garde), KKK, 1,2);
   
    pause 
    close all
end



return




X = [lesreponsesL2RAT]; X(find(isinf(X)))=nan;
F1 = fillmissing(X,'spline');

X = [lesreponsesNDEI]; X(find(isinf(X)))=nan;
F2 = fillmissing(X,'spline');

% X = [lesreponsesMAXERR]; X(find(isinf(X)))=nan;
% F3 = fillmissing(X,'nearest');

X = [lesreponsesLOGL]; X(find(isinf(X)))=nan;
F4 = fillmissing(X,'spline');

X = [F1 F2  F4];

KKK=10;

options = statset('MaxIter',1000,'RobustWgtFun','cauchy')
Sigma = 'full';%'diagonal';%,'full'};
SharedCovariance = true; %{true,false};
gmfit = fitgmdist(X,KKK,'CovarianceType',Sigma,...
    'SharedCovariance',SharedCovariance,'Options',options);
c1 = cluster(gmfit,X);


% f=find(cc==6);
% xo=0; figure,for i=1:length(f), x=xo+linspace(0,1,length(salvesOut{f(i)})); hold on, plot(x,salvesOut{f(i)}), xo=x(end); pause, end


% ON AFFICHE LES RESULTATS SANS CONSIDERE Bbarre
global reordonnerPourVisu, reordonnerPourVisu=true;
cc=duree_vs_amplitude(c1, X, T(garde), KKK, 2, 1);
tabulate(cc)

visualiser_resultats_shm_clust_morceaux_2(cc,T,KKK,[]);
superpose_energie(Features(:,22),T,1);

visualiser_resultats_shm_clust_morceaux_2(cc,defcarac,KKK,[]);
superpose_energie(Features(:,22),defcarac,1);

densites_amplitude(cc, Features(garde,13), 1) %F(garde,5), 1);
    
    tracer_energie_dans_cluster(cc, Features(garde,11), T(garde));
    
    tracer_energie_dans_cluster_normalise(cc, Features(garde,11), T(garde));
    
    densites_amplitude(cc, Features(garde,29), 1) %F(garde,5), 1);
    
    duree_vs_amplitude(cc, [Features(garde,13) log10(Features(garde,12))], T(garde), KKK, 1,2);
    
    duree_vs_amplitude(cc, [Features(garde,29) Features(garde,24)], T(garde), KKK, 1,2);
    












visualiser_resultats_shm_clust_morceaux_2(cc,defCaracLoca,KKK,[]);
superpose_energie(Features(indices,11),defCaracLoca,1);


X=[lesreponsesL2RAT];
[a b]=min(exp(X)./repmat(sum(exp(X),2),1,size(X,2)),[],2);
X(find(isinf(X)))=-inf;

X=[lesreponsesLOGL];
[a b]=max(exp(X)./repmat(sum(exp(X),2),1,size(X,2)),[],2);
X(find(isinf(X)))=-inf;

X=[lesreponsesNDEI];
[a b]=min(exp(X)./repmat(sum(exp(X),2),1,size(X,2)),[],2);
X(find(isinf(X)))=-inf;


X=[lesreponsesL2RAT lesreponsesLOGL lesreponsesMAXERR lesreponsesNDEI lesreponsesPSNR lesreponsesMSE];

%X=[lesreponsesL2RAT(:,1) lesreponsesLOGL(:,1) lesreponsesMAXERR(:,1) lesreponsesNDEI(:,1) lesreponsesPSNR(:,1) lesreponsesMSE(:,1)];

% X=[lesreponsesLOGL(:,1) lesreponsesNDEI(:,1)];
aretirer = find(isinf(sum(X,2)) | isnan(sum(X,2)) | (abs(sum(abs(X),2)-eps)<=eps));
garde = 1:size(X,1); garde(aretirer)=[];
disp(sprintf('Retire %d pts / %d',length(aretirer),size(X,1)))

X(aretirer,:)=[];


KKK=4;

% que la likelihood
% global normaliser, normaliser=true;
% c1 = GK_multi_init(X, 1, KKK, 1, [] ,[], 'gmm1');%, 'fcm');
% [c1 c2]=kmeans(zscore(X), KKK, 'replicates', 100);
% [c1 c2]=kmeans(zscore(X), KKK, 'replicates', 100);


options = statset('MaxIter',1000)
Sigma = 'full';%'diagonal';%,'full'};
SharedCovariance = true; %{true,false};
gmfit = fitgmdist(X,KKK,'CovarianceType',Sigma,...
    'SharedCovariance',SharedCovariance,'Options',options);
c1 = cluster(gmfit,X);


% ts les criteres
% global normaliser, normaliser=true;
% c1 = GK_multi_init(X, 1:size(X,2), KKK, 1:size(X,2), [] ,[], 'kmeans');%, 'fcm');

% ON AFFICHE LES RESULTATS SANS CONSIDERE Bbarre
global reordonnerPourVisu, reordonnerPourVisu=true;
cc=duree_vs_amplitude(c1, X, T(garde), KKK, 2, 1);
tabulate(cc)

% clusterassocie = cc(referencesSalves)
% f=find(cc==clusterassocie)
% figure, for i=1:length(f), plot(salvesOut{f(i)}), pause, end

visualiser_resultats_shm_clust_morceaux_2(cc,T(garde),KKK,[]);
superpose_energie(Features(garde,11),T(garde),1);

densites_amplitude(cc, Features(garde,13), 1) %F(garde,5), 1);

figure, tracer_energie_dans_cluster(cc, Features(garde,11), T(garde));

tracer_energie_dans_cluster_normalise(cc, Features(garde,11), T(garde));

    
densites_amplitude(cc, Features(garde,29), 1) %F(garde,5), 1);

duree_vs_amplitude(cc, [Features(garde,13) log10(Features(garde,12))], T(garde), KKK, 1,2);

duree_vs_amplitude(cc, [Features(garde,29) Features(garde,24)], T(garde), KKK, 1,2);


return

k = 5;
Sigma = {'diagonal','full'};
nSigma = numel(Sigma);
SharedCovariance = {true,false};
SCtext = {'true','false'};
nSC = numel(SharedCovariance);
d = 500;
x1 = linspace(min(X(:,1)) - 2,max(X(:,1)) + 2,d);
x2 = linspace(min(X(:,2)) - 2,max(X(:,2)) + 2,d);
[x1grid,x2grid] = meshgrid(x1,x2);
X0 = X(1,:);%[x1grid(:) x2grid(:)];
threshold = sqrt(chi2inv(0.99,2));
options = statset('MaxIter',1000); % Increase number of EM iterations

figure;
c = 1;
for i = 1:nSigma;
    for j = 1:nSC;
        gmfit = fitgmdist(X,k,'CovarianceType',Sigma{i},...
            'SharedCovariance',SharedCovariance{j},'Options',options);
        clusterX = cluster(gmfit,X);
        mahalDist = mahal(gmfit,X0);
        subplot(2,2,c);
        h1 = gscatter(X(:,1),X(:,2),clusterX);
        hold on;
        for m = 1:k;
            idx = mahalDist(:,m)<=threshold;
            Color = h1(m).Color*0.75 + -0.5*(h1(m).Color - 1);
            h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
            uistack(h2,'bottom');
        end
        plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
        title(sprintf('Sigma is %s, SharedCovariance = %s',...
            Sigma{i},SCtext{j}),'FontSize',8)
        legend(h1,{'1','2','3'});
        hold off
        c = c + 1;
    end
end


return


% on segmente l'indicateur
X=lesreponsesLOGL;
aretirer = find(isinf(sum(X,2)) | isnan(sum(X,2)) );
garde = 1:size(X,1); garde(aretirer)=[];
disp(sprintf('Retire %d pts / %d',length(aretirer),size(X,1)))
X(aretirer,:)=[];

caracSelectionnes = Features(garde,:);

clear eva
k=1;
rng('default');  % For reproducibility
methods={'kmeans','gmdistribution','linkage'};
indices={'CalinskiHarabasz','silhouette','DaviesBouldin','gap'};% gap peut etre long!!!!
for i=1:length(methods)
    for j=1:length(indices)
        disp(sprintf('En cours : method %s - indic. %s',methods{i},indices{j}))
        eva{k} = evalclusters(zscore(X),methods{i},indices{j},'Klist',[3:8]);
        k=k+1;
    end
end

% trace clusters
for k=1:length(eva)
    
    close all
    
    %     %%%%%%%%%%%%%%%%%% trace logCSCA
    %     figure;
    %     hold on
    %     for i=1:eva{k}.OptimalK
    %         f=find(eva{k}.OptimalY==i);
    %         v=zeros(length(eva{k}.OptimalY),1);
    %         v(f)=1;
    %         plot(temps,log(cumsum(v)))
    %         title(sprintf('Cumul clusters vs tps -> methode = %s, critere = %s, K=%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK))
    %         xlabel ('temps(s)')
    %         ylabel ('logCSCA')
    %     end
    %
    %     %saveas(gcf,sprintf('cumul_clusters_vs_tps_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'fig')
    %     saveas(gcf,sprintf('cumul_clusters_vs_tps_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'png')
    
    %%%%%%%%%%%%%%%%%%%%%%%%% trace duree vs amplitude
    global reordonnerPourVisu, reordonnerPourVisu=true;
    cc=duree_vs_amplitude(eva{k}.OptimalY, [caracSelectionnes(:,13) log10(caracSelectionnes(:,12))],...
        caracSelectionnes(:,5), eva{k}.OptimalK, 1, 2);
    
    saveas(gcf,sprintf('dureeAmpls_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'png')
    
    %%%%%%%%%%%%%%%%%%%%%%%%% trace logCSCA superpose energie absolue
    visualiser_resultats_shm_clust_morceaux_2(cc,caracSelectionnes(:,5),eva{k}.OptimalK,[]);
    superpose_energie(caracSelectionnes(:,22),caracSelectionnes(:,5),1);
    
    saveas(gcf,sprintf('logCSCA_vs_tps_absE_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'png')
    
    %%%%%%%%%%%%%%%%%%%%%%%%% trace logCSCA superpose energie MARSE
    visualiser_resultats_shm_clust_morceaux_2(cc,caracSelectionnes(:,5),eva{k}.OptimalK,[]);
    superpose_energie(caracSelectionnes(:,11),caracSelectionnes(:,5),1);
    
    saveas(gcf,sprintf('logCSCA_vs_tps_Emarse_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'png')
    
    %%%%%%%%%%%%%%%%%%%%%%%%% trace densites amplitudes
    densites_amplitude(cc, caracSelectionnes(:,13), 1) %F(garde,5), 1);
    
    saveas(gcf,sprintf('densite_amplitude_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'png')
    
    %%%%%%%%%%%%%%%%%%%%%%%%% trace energie marse
    tracer_energie_dans_cluster(cc, caracSelectionnes(:,11), caracSelectionnes(:,5));
    
    saveas(gcf,sprintf('energieParCluster_marse_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'png')
    
    %%%%%%%%%%%%%%%%%%%%%%%%% trace energie abs
    tracer_energie_dans_cluster(cc, caracSelectionnes(:,22), caracSelectionnes(:,5));
    
    saveas(gcf,sprintf('energieParCluster_abs_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'png')
    
    %%%%%%%%%%%%%%%%%%%%%%%%% trace energie - strength
    tracer_energie_dans_cluster(cc, caracSelectionnes(:,21), caracSelectionnes(:,5));
    
    %saveas(gcf,sprintf('Duree_vs_amplitude_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'fig')
    saveas(gcf,sprintf('strength_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'png')
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% trace energie marse
    tracer_energie_dans_cluster_normalise(cc, caracSelectionnes(:,11), caracSelectionnes(:,5));
    
    saveas(gcf,sprintf('energieParCluster_marse_normalisee_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'png')
    
    %%%%%%%%%%%%%%%%%%%%%%%%% trace energie abs
    tracer_energie_dans_cluster_normalise(cc, caracSelectionnes(:,22), caracSelectionnes(:,5));
    
    saveas(gcf,sprintf('energieParCluster_abs_normalisee_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'png')
    
    %%%%%%%%%%%%%%%%%%%%%%%%% trace energie - strength
    tracer_energie_dans_cluster_normalise(cc, caracSelectionnes(:,21), caracSelectionnes(:,5));
    
    %saveas(gcf,sprintf('Duree_vs_amplitude_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'fig')
    saveas(gcf,sprintf('strength_normalisee_meth_%s_crit_%s_K%d',eva{k}.ClusteringFunction,eva{k}.CriterionName,eva{k}.OptimalK),'png')
    
    
end


return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dans aretirer on a les salves qui ne sont pas de tout dans le modele
c2=(KKK+1)*ones(size(Features,1),1);
c2(garde)=cc;
c2(aretirer)=KKK+1;
KKK=KKK+1;
tabulate(c2)

global reordonnerPourVisu, reordonnerPourVisu=false;
duree_vs_amplitude(c2, [Features(:,13) log10(Features(:,12))], T, KKK, 1,2);

visualiser_resultats_shm_clust_morceaux_2(c2,T,KKK,[]);
superpose_energie(Features(:,11),T,1);

densites_amplitude(c2, Features(:,13), 1) %F(garde,5), 1);
figure,tracer_energie_dans_cluster(c2, Features(:,11), T);

duree_vs_amplitude(c2, [Features(:,29) (Features(:,24))], T, KKK, 1,2);

f=find(c2==KKK);%rejet
figure, for i=f(:)',
    stem((salvesOut{1,lsa}-my)/stdy,'r'),
    hold on, plot((salvesOut{1,i}-my)/stdy); i, pause, clf,
end

f=find(c2==1);%rejet
figure, for i=f(:)',
    stem(salvesOut{1,lsa}/max(salvesOut{1,lsa}),'r'),
    hold on, plot(salvesOut{1,i}/max(salvesOut{1,i})); i, pause, clf,
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%


lsa=1;
z=salvesOut{1,lsa}';% un autre

figure,plot(z)
lemodel=1;
W=400;
L=nan*ones(length(z),1);
MS=nan*ones(length(z),1);
figure, hold on,

for Rt=W+1:30:length(z)
    
    y = z(Rt-W:Rt) - models(lemodel).my;
    y=y(:);
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
        [alp,~, gamma, loglik] = fwdback_phmm_TD(models(lemodel).Pi, ...
            models(lemodel).A, obslik, ones(T,models(lemodel).Q));
        
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
        
        U=loglik / models(lemodel).LL(end);
        if U>1, U=-inf; end
        L(Rt-W:Rt)=U;
        MS(Rt-W:Rt) = sum((y-signalcomb).^2);
        clf,
        subplot(411), plot(z),
        subplot(412),plot(L,'o'), hold on, plot(Rt,L(Rt),'ro','linewidth',3,'markersize',10), xlim([0 length(z)]),drawnow
        subplot(413), plot(y), hold on, plot(signalcomb,'ro')
        subplot(414), plot(MS), hold on, plot(Rt,MS(Rt),'ro','linewidth',3,'markersize',10), xlim([0 length(z)]),drawnow
        
    end
    [Rt, length(z)]
end
figure,plot(L)



% lit les fichiers 1a17 de FA21 traités par Sylavin
% clear all
% de=1;
% a=17;
% ch='/home/ramassoemmanuel/Documents/Dropbox/6_BENCHMARK/NouvellesDonneesDinh_x45_8plis/Fatigue/FA21_traiteParSylvainDec2015';
% H={}; D=[]; H2={};
% for i=de:1:a
%     clear salvesOut
%     load([ch '/salvesOut_PT' num2str(i)]);
%     for j=1:length(salvesOut(1,:))
%         H{1,end+1} = salvesOut{1,j};
%     end
%     clear carac
%     load([ch '/carac_PT' num2str(i)]);
%     f = find(carac(:,8)==1);% que la voie 1 qui recoit toujours
%     if f(1)==1, cpt=1; else cpt=0; end
%     for j=2:length(f)
%         if f(j)==f(j-1)+1, cpt=cpt+1; end
%         if cpt==3
%             break
%         end
%     end
%     if cpt~=3, error('??'); end
%     D = [D ; carac(f(j),:)];
%     H2{end+1} = salvesOut{1,f(j)};
%     % D = [D ; carac];
% end
% salvesOut = H2;

