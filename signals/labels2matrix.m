function m=labels2matrix(labels, K)

n=length(labels);
m=zeros(K,n);
m(labels(:)' + [0:K:K*n-1])=1; 
m=m';