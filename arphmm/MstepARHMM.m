function [B,sigma,obslik]=MstepARHMM(data,gamma,predata,k,P,Q,T)

if 1
    
    % B
    B=zeros(Q,P);    
    sigma=cell(Q,1);
    obslik=zeros(T,Q);
    
    for i=1:Q
        res1=zeros(1,P);
        res2=zeros(P,P);
        for t=1:T
            res1=res1+(gamma(t,i)*data(t,:)*predata{t}');
            res2=res2+(gamma(t,i)*predata{t}*predata{t}');
        end
        B(i,:)=-res1*pinv(res2);
        
        % SIGMA
        sigmastep = 0; if k>1, error('??'); end %zeros(k);
        for t=1:T
            sigmastep=sigmastep+(gamma(t,i)*((data(t,:)+(B(i,:)*predata{t}))'*(data(t,:)+(B(i,:)*predata{t}))))/sum(gamma(:,i));
        end
        sigma{i}=sigmastep;
        
        % Calcul of b    
        for t=1:T
            obslik(t,i)=mvnpdf(data(t,:)+(B(i,:)*predata{t}),0,sigma{i});
        end
        
    end
               
else
    
    % B
    B=zeros(Q,P);
    for i=1:Q
        res1=zeros(1,P);
        res2=zeros(P,P);
        for t=1:T
            res1=res1+(gamma(t,i)*data(t,:)*predata{t}');
            res2=res2+(gamma(t,i)*predata{t}*predata{t}');
        end
        B(i,:)=-res1*pinv(res2);
    end
    
    % SIGMA
    sigma=cell(Q,1);
    for i=1:Q
        [sigmastep]=zeros(k);
        for t=1:T
            sigmastep=sigmastep+(gamma(t,i)*((data(t,:)+(B(i,:)*predata{t}))'*(data(t,:)+(B(i,:)*predata{t}))))/sum(gamma(:,i));
        end
        sigma{i}=sigmastep;
    end
    
    % Calcul of b
    obslik=zeros(T,Q);
    for i=1:Q
        for t=1:T
            obslik(t,i)=mvnpdf(data(t,:)+(B(i,:)*predata{t}),0,sigma{i});
        end
    end
    % A VOIR
    %obslik=mk_stochastic(obslik);
end
end
