function [gamma2, NDEIf, msemf, signalt, signalcomb, B, sigma, ...
    A, Pi, alpha2, gammaHMM, alphaHMM, NDEI, msem, LL]= ...
   AUTOREGRESSIVE(data,pl, Q, P)

% P = order of the mAR
% Q = number of states
% T = number of observations
% k = dimensions of signal

k = 1;
T = length(data);
USE_CRISP_GAMMA = false;

disp('Predata...'), tic;
% Zk = Vector with previous p data
% predata{i} = bloc de P data, décalage de 1
predata=cell(T,1);
%datmat=zeros(T+P,k)+0.01*(rand(T+P,k));
datmat=zeros(T+P,k);
datmat(1:P,:)=data(1:P,:);
datmat(P+1:T+P,:)=data;
for t=1:T
    predata{t} = zeros(P,size(data,2));
   for pp=1:P
      predata{t}(pp,:)=datmat(t+P-pp,:);
   end
end
disp(sprintf('OK en %f (%d pts, predc = %d)',toc,T,P));

%%% APLICATION MODEL CONTINUE AVEC PL=ones

%pl=ones(T,Q);
% pl=zeros(T,Q)+0.1;
% pl(1:125,1)=1;
% pl(126:250,2)=1;
% pl(251:375,3)=1;
% pl(376:500,4)=1;

CONSTEPS = 1e-5;

disp('PHMM...'), tic;

parametersAlgorithm = setHMMDefaultParameters;
parametersAlgorithm.nessai=5;
[parametersHMM, outputsInference] = ...
    phmm_gauss_mix_learn(data, pl, Q, 1, parametersAlgorithm);
Pi=parametersHMM.Pif;
A=parametersHMM.Af;
gamma=outputsInference.gamma;
alpha=outputsInference.alpha;
%[~,~,Pi,A,~,~,gamma,alpha]=phmm_gauss(data,pl,nessai,idiag,iltr,visu,iplot,init);
Pi=Pi+CONSTEPS;
A=A+CONSTEPS;
gammaHMM = gamma;
alphaHMM = alpha; 
disp(sprintf('OK en %f',toc));

disp('Mstep1...'), tic;
[B,sigma,obslik]=MstepARHMM(data,gamma,predata,k,P,Q,T);
disp(sprintf('OK en %f',toc));

disp('Autoreg...'), tic;
signal=autoregres(predata,B,sigma);
disp(sprintf('OK en %f',toc));
% eq 1, evalue yk dans chaque etat, figure,plot(signal{i})

ssmat=cell2mat(signal);
signalt=reshape(ssmat,T,Q);

% Signal with max gamma

% le signal de sortie est donné par signal{k} avec le max de gamma
if USE_CRISP_GAMMA
   [~,b]=max(gamma,[],2);
   signalcomb=zeros(T,1);
   for t=1:T
      signalcomb(t)=signalt(t,b(t));
   end
else
   signalcomb=zeros(T,1);
   for t=1:T
      for j=1:Q
         signalcomb(t)=signalcomb(t) + gamma(t,j)*signalt(t,j);
      end
   end
end
msem=mean(abs(data-signalcomb).^2);
NDEI=sqrt(msem)/std(data);
disp(sprintf('Erreurs initiales : E=%f, NDEI=%f',msem,NDEI))

%Plausibilities
% pl=zeros(T,Q)+0.1;
% pl(1:125,mu0)=1;
% pl(126:250,mu1)=1;
% pl(251:375,mu2)=1;
% pl(376:500,mu3)=1;

% Turning the algorithm n times
disp('Running optim...'), tic;
counter=1; converged=false; decrease=false;
LL = []; loglik=-inf;
while counter<100 && ~converged && ~decrease
   
    LL=[LL loglik];
   
    % Etape E
   [alpha2,~, gamma2, loglik, xi] = fwdback_phmm_mix(Pi, A, obslik, pl);
   A = xi;%sum(xi,1);
   A = squeeze(mk_stochastic(A+CONSTEPS));% AJOUT
   Pi = normalise(gamma2(1,:)+CONSTEPS);% AJOUT
   
   gamma3 = gamma2;
   if USE_CRISP_GAMMA
       gamma3(:,:)=0;
       for i=1:size(gamma3,1)
           [a b]=max(gamma2(i,:));
           gamma3(i,b)=1;
       end
   end
   
   % M step
   [B,sigma,obslik]=MstepARHMM(data,gamma3,predata,k,P,Q,T);
   signal=autoregres(predata,B,sigma);
   
   ssmat=cell2mat(signal);
   signalt=reshape(ssmat,T,Q);
   
   % le signal de sortie est donné par signal{k} avec le max de gamma
   if USE_CRISP_GAMMA
      [~,b]=max(gamma2,[],2);
      signalcomb=zeros(T,1);
      for t=1:T
         signalcomb(t)=signalt(t,b(t));
      end
   else
      signalcomb=zeros(T,1);
      for t=1:T
         for j=1:Q
            signalcomb(t)=signalcomb(t) + gamma2(t,j)*signalt(t,j);
         end
      end
   end
   
   msemf=mean(abs(data-signalcomb).^2);
   NDEIf=sqrt(msemf)/std(data);
   
   %figure(20),clf,plot(data-signalcomb,'ro')
   
   disp(sprintf('%d : mse=%f, ndei=%f, LL=%f', counter, msemf, NDEIf, loglik))
    
%    if NDEIf<NDEI
%       NDEI=NDEIf;
%       signaltf=signalt;
%       signalcombf=signalcomb;
%    end
   counter=counter+1;
   [converged, decrease] = em_converged(loglik, LL(end), 1e-4, 1);
   if counter<=3, converged=false; decrease=false; end
   

end
disp(sprintf('FIN en %f',toc));

%Ploting results

% figure, subplot(211),
% %plot(data), hold on
% plot(signalt);%,'r'), plot(signalt(:,2),'g'), plot(signalt(:,3),'y')
% title('Signal dans chaque état') 
% subplot(212)
% plot(data), hold on, plot(signalcomb,'r')
% title('Signal predit') 

end


