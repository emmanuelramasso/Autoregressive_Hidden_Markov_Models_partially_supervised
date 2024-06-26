function [lesreponsesLOGL_tmp, lesreponsesMSE_tmp, ...
    lesreponsesNDEI_tmp, lesreponsesMAXERR_tmp, ...
    lesreponsesPSNR_tmp, lesreponsesL2RAT_tmp] = ...
    appliquerModeles(models, z)


% ici on applique les modeles
lesreponsesLOGL_tmp = zeros(1,length(models));
lesreponsesMSE_tmp = zeros(1,length(models));
lesreponsesNDEI_tmp = zeros(1,length(models));
lesreponsesMAXERR_tmp = zeros(1,length(models));
lesreponsesPSNR_tmp = zeros(1,length(models));
lesreponsesL2RAT_tmp = zeros(1,length(models));


for lemodel=1:length(models)
    
    try % cas ou la taille de la salve est trop courte wrt apprentissage
        
        y = z;
        
        %%% INFERE
        y = y - models(lemodel).my;
        y = y / models(lemodel).stdy;
        Q = models(lemodel).Q;
        T =length(y);
        P = models(lemodel).predc;
        
        warning('off')
        predata = fliplr(hankel(y,y(1:P)));
        warning('on')
        predata(end-P+2:end,:) = [];
        predata = [zeros(P-1,P) ; predata];

        % Calcul of b
        obslik=zeros(T,Q); obslik(1,:)=1; % init predata
        %for i=1:Q
        %    for t=2:T
        %        obslik(t,i) = mvnpdf(y(t,:)+models(lemodel).B(i,:)*predata(t-1,:)',0,models(lemodel).S{i});
        %    end
        %end
        
        for i=1:Q
            X = y(2:end)+predata(1:end-1,:)*models(lemodel).B(i,:)';
            obslik(:,i)=[1 ; mvnpdf(X,0,models(lemodel).S{i})];
        end
        
        if isempty(find(sum(obslik,2)==0,1))
            
            if ~models(lemodel).usePriorGamma
                [~,~, gamma, loglik] = fwdback_phmm_mix(models(lemodel).Pi, ...
                    models(lemodel).A, obslik, ones(T,models(lemodel).Q));
            else
                [~,~, gamma, loglik] = fwdback_phmm_mix(models(lemodel).Pi, ...
                    models(lemodel).A, obslik, models(lemodel).gamma);
            end
            
            assert(isempty(find(isnan(gamma),1)))
            
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
            %figure,plot(signalcomb,'rx'), hold on, stem(y)
            
            [PSNR,MSE,MAXERR,L2RAT] = measerr(y,signalcomb);
            
            lesreponsesLOGL_tmp(1,lemodel) = loglik / T;
            lesreponsesMSE_tmp(1,lemodel) = MSE; %mean((signalcomb-y).^2);
            lesreponsesNDEI_tmp(1,lemodel) = sqrt(lesreponsesMSE_tmp(1,lemodel))/std(y);
            lesreponsesMAXERR_tmp(1,lemodel) = PSNR;
            lesreponsesPSNR_tmp(1,lemodel) = MAXERR;
            lesreponsesL2RAT_tmp(1,lemodel) = L2RAT;
            
        else
            lesreponsesLOGL_tmp(1,lemodel) = -inf;
            lesreponsesMSE_tmp(1,lemodel) = inf;
            lesreponsesNDEI_tmp(1,lemodel) = inf;
            lesreponsesMAXERR_tmp(1,lemodel) = inf;
            lesreponsesPSNR_tmp(1,lemodel) = inf;
            lesreponsesL2RAT_tmp(1,lemodel) = inf;
        end
    catch
        lesreponsesLOGL_tmp(1,lemodel) = -inf;
        lesreponsesMSE_tmp(1,lemodel) = inf;
        lesreponsesNDEI_tmp(1,lemodel) = inf;
        lesreponsesMAXERR_tmp(1,lemodel) = inf;
        lesreponsesPSNR_tmp(1,lemodel) = inf;
        lesreponsesL2RAT_tmp(1,lemodel) = inf;
    end
    % [lasalve, length(salvesOut), lemodel, length(models)]
    
end