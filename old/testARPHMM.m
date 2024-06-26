clear all
% load y;
addpath PHMM/
addpath PHMM/utils

predc = 5;

if 0
   % ############################
    % charge des donnees ressemblant a des salves mises bout a bout
    % ############################
    
   load datatest/Descripteurs_voie1_allData.mat
   % figure, for uu=[45 50 93 187 202 220 257 263 301 299 294 288 277 275 309 314 317 320 321 322 323 327 330], x = savesalve{uu};, ,plot(x), uu, pause, clf, end
   s1=savesalve{314}; lnS1 = length(s1);
   s2=savesalve{317}; lnS2 = length(s2);
   s3=savesalve{320}; lnS3 = length(s3);
   n1=0; n2=0;
   %y=[normalise_0_1(s1) ; n1+normalise_0_1(s2) ; n1+normalise_0_1(s3)];
   y=[(s1) ; n1+(s2) ; n2+(s3)];
   gtruth = [ones(lnS1,1) ; 2*ones(lnS2,1); 3*ones(lnS3,1)];
   K=3;
   
elseif 1
   % ############################
   % charge des donnees ressemblant a des salves mises bout a bout et "deformées"
    % ############################
    
   load datatest/zz
   y=sig2(:);
   lnS1=480; lnS2=2430-480; lnS3=length(y)-2430;
   gtruth = [ones(lnS1,1) ; 2*ones(lnS2,1); 3*ones(lnS3,1)];
   K=3;
   
end
   
K
figure,subplot(211), plot(y), subplot(212), plot(gtruth)
SEUILHAUT=1.2;
SEUILBAS = -0.3;
D = 15;
T=length(y);


% ############################
% apprentissage d'un ARPHMM
% sans apriori sur les etats
% ############################

pl=ones(length(y),K);
[gamma2_vac, NDEI_vac, qualitymf_vac, signalt_vac, signalcomb_vac] =...
   arphmm(y, pl, K, predc);
mean((y(D:end)-signalcomb_vac(D:end)).^2)
[a b]=max(gamma2_vac,[],2); 
[AR,~,~,~, C]=valid_RandIndex(gtruth, b)

figure,subplot(211), plot(y), subplot(212), plot(gtruth), hold on, plot(b,'r')

% ici on impose un a priori aleatoire sur les etats caches
gtruth2=randi(3,[length(y),1]);

pl0=labels2matrix(gtruth2, K);
 [gamma2_sup, NDEI_sup, qualitymf_sup, signalt_sup, ...
         signalcomb_sup, B_sup, sigma_sup, A_sup, Pi_sup] =...
         arphmm(y, pl0,K,predc);      
mean((y(D:end)-signalcomb_sup(D:end)).^2)
[a b]=max(gamma2_sup,[],2); 
 [AR,~,~,~, C]=valid_RandIndex(gtruth2, b)

 
% ############################
% on test ici l'influence des a priori sur la reconnaissances des etats
% dans certains cas "non supervisés" les a priori ne servent a rien
% dans d'autres cas "plutot supervisé" les a priori peuvent aider a converger
E=[]; ARI=[];
lesnp=[0:0.1:1];
l=0; i=1; cont=true;
while cont
   try
      if 0
         pl=pl0;
         pl(find(pl==0))=lesnp(i);
      else
         % Bruitage des labels
         if lesnp(i)==0,
            perr=zeros(T,1);
         elseif lesnp(i)==1,
            perr=ones(T,1);
         else,
            [a,b]=param_beta(lesnp(i),(0.2).^2);
            perr=betarnd(a,b,T,1);
         end;
         [pluncertain,y1,plnoisy]=add_noise1(gtruth,perr,K);
      end
      plType = pluncertain; %plnoisy
      [gamma2_01, NDEI_01, qualitymf_01, signalt_01, ...
         signalcomb_01, B, sigma, A, Pi] =...
         arphmm(y, plType,K, predc);
      e=mean((y(D:end)-signalcomb_01(D:end)).^2);
      %r=sqrt(e)/std(y);
      E = [E ; [e lesnp(i)]];
      [a b]=max(gamma2_01,[],2);
      ARI = [ARI ; valid_RandIndex(gtruth, b)];
      [E, ARI]
      l=0; i=i+1; 
   catch
      disp('Retry...'), l=l+1;
      if l>20, error('PROBLEM DE CONVERGENCE ??'), end
   end
   if i>length(lesnp), cont=false; end
end
figure,plot(lesnp,E(:,1)), title('MSE prediction')
%figure,plot(lesnp,E(:,2)), title('NDEI prediction')
figure,plot(lesnp,ARI), title('ARI')
figure,plot(lesnp,(1-ARI).*(E(:,1))), title('(1-ARI)*MSE')
%figure,plot(lesnp,(1-ARI).*E(:,2)), title('(1-ARI)*NDEI')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
