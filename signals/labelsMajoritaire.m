function [L1, L2, F1, F2, doute]=labelsMajoritaire(X)

L1=zeros(size(X,1),1); L2=L1; F1=L1; F2=L1; doute=L1;

for i=1:size(X,1)
    
    v=X(i,:);
    cc=tabulate(v(:));
    f=find(cc(:,2)==0);
    cc(f,:)=[];
    
    [a1 b1]=sort(cc(:,2),'descend');
    
    if length(b1)>1
        doute(i) = a1(1)==a1(2);
        L2(i)=cc(b1(2),1);
        F2(i)=cc(b1(2),2);%/100;
    end
    
    L1(i)=cc(b1(1),1);
    F1(i)=cc(b1(1),2);%/100;
    
end