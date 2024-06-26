
clear all; close all
% load '/home/emmanuel.ramasso/Documents/Dropbox/6_BENCHMARK/Biocomposites/FLAX_2017/P1_UD_07_TS0_10/file_2017_04_13_08_47_27953_Features_et_salves_0.mat'
addpath PHMM/
addpath PHMM/utils/
addpath utils/

[Features,tempsSalves,ordre] = reordonneTri(Features);
salvesOut=salvesOut(ordre);


% figure, for i=1:length(salvesOut), clf, plot(salvesOut{i}), i, pause(1), end

% numero des salves choisies pour apprendre les modeles
% referencesSalves=[169]%61 120];
% salvesModeles = salvesOut(referencesSalves);
referencesSalves=[6 61 120]%61 120];
salvesModeles = salvesOut(referencesSalves);
salvesModeles{1}=salvesModeles{1}(345:704);
figure,plot(salvesModeles{1})

clear models
k=1
predc = 4;
Q = 5;
D = 5

for lsa=1:length(salvesModeles)%    48;
    
    z=salvesModeles{lsa};%(1:400)';% un autre
    
    figure,plot(z)
    
    y = z;
    my = mean(y); stdy = std(y);
    y=(y-my)/stdy;
    
    if 1
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

% ici on applique les modeles
lesreponsesLOGL = zeros(size(salvesOut,2),length(models));
lesreponsesMSE = zeros(size(salvesOut,2),length(models));
lesreponsesNDEI = zeros(size(salvesOut,2),length(models));

for lasalve=1:1:length(salvesOut)
    
    y = salvesOut{lasalve}(:);
    z = y;
    
    for lemodel=1:length(models)        
        
        try % cas ou la taille de la salve est trop courte wrt apprentissage
            
            y =z;
            
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
                
                lesreponsesLOGL(lasalve,lemodel) = loglik / T;
                lesreponsesMSE(lasalve,lemodel) = mean((signalcomb-y).^2);
                lesreponsesNDEI(lasalve,lemodel) = sqrt(lesreponsesMSE(lasalve,lemodel))/std(z);
            else
                lesreponsesLOGL(lasalve,lemodel) = -inf;
                lesreponsesMSE(lasalve,lemodel) = inf;
                lesreponsesNDEI(lasalve,lemodel) = inf;
            end
        catch
                lesreponsesLOGL(lasalve,lemodel) = -inf;
                lesreponsesMSE(lasalve,lemodel) = inf;
                lesreponsesNDEI(lasalve,lemodel) = inf;
        end
        [lasalve, length(salvesOut), lemodel, length(models)]
    end
end


% % analyse des resultats 
% [A B]=max(lesreponsesLOGL,[],2); figure,plot(B)
% [A B]=min(lesreponsesMSE,[],2);figure,plot(B)
% [A B]=min(lesreponsesNDEI,[],2); figure,plot(B)
% figure,imagesc(lesreponsesLOGL)
% figure,imagesc(lesreponsesMSE)
% figure,imagesc(lesreponsesNDEI)

% data
noms2={'LOGL','MSE','NDEI'};
% global noms, noms={'LOGL','NDEI'};
X=[lesreponsesLOGL(:,:) lesreponsesMSE(:,:) lesreponsesNDEI(:,:)];

% ############### POINT IMPORTANT
% on ramene tt a la reference
% par exemple la premiere reference 
lsa = referencesSalves(1)

ref1=X(lsa,1); X(:,1)=X(:,1)./ref1;
ref2=X(lsa,2); X(:,2)=(X(:,2)-ref2)./ref2;
ref3=X(lsa,3); X(:,3)=(X(:,3)-ref3)./ref3;

X(lsa,:) % 1ere colonne => =1 ??
figure,plot(X), legend(noms2)


% X=[lesreponsesLOGL(:,1) lesreponsesNDEI(:,1)];
aretirer = find(isinf(sum(X,2)) | isnan(sum(X,2)) );
garde = 1:size(X,1); garde(aretirer)=[];
disp(sprintf('Retire %d pts / %d',length(aretirer),size(X,1)))

X(aretirer,:)=[];

T=tempsSalves;

KKK=10;

% que la likelihood 
% global normaliser, normaliser=true;
% c1 = GK_multi_init(X, 1, KKK, 1, [] ,[], 'gmm1');%, 'fcm');
% [c1 c2]=kmeans(zscore(X), KKK, 'replicates', 100);
[c1 c2]=kmeans(zscore(X(:,1)), KKK, 'replicates', 100);


% ts les criteres
% global normaliser, normaliser=true;
% c1 = GK_multi_init(X, 1:size(X,2), KKK, 1:size(X,2), [] ,[], 'kmeans');%, 'fcm');

% ON AFFICHE LES RESULTATS SANS CONSIDERE Bbarre
global reordonnerPourVisu, reordonnerPourVisu=true;
cc=duree_vs_amplitude(c1, X, T(garde), KKK, 3, 1); 
tabulate(cc)

% clusterassocie = cc(referencesSalves)
% f=find(cc==clusterassocie)
% figure, for i=1:length(f), plot(salvesOut{f(i)}), pause, end

visualiser_resultats_shm_clust_morceaux_2(cc,T(garde),KKK,[]);
superpose_energie(Features(garde,11),T(garde),1);

densites_amplitude(cc, Features(garde,13), 1) %F(garde,5), 1);

figure, tracer_energie_dans_cluster(cc, Features(garde,11), T(garde));

densites_amplitude(cc, Features(garde,29), 1) %F(garde,5), 1);

duree_vs_amplitude(cc, [Features(garde,13) log10(Features(garde,12))], T(garde), KKK, 1,2);

duree_vs_amplitude(cc, [Features(garde,29) Features(garde,24)], T(garde), KKK, 1,2);



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

