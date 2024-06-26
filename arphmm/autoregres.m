function [signal, signalbruite]=autoregres(predata,B,Sigma)
%AUTOREGRESSIVE TEST
Q=length(Sigma);
signal=cell(Q,1);
signalbruite=signal;
w=cell(Q,1);
T=length(predata);
for i=1:Q
    signalstep=zeros(T,1);
    for t=1:T
        signalstep(t)=-(B(i,:)*predata{t});
    end
    w{i}=mvnrnd(0,Sigma{i},T);
    signal{i}=signalstep;
    signalbruite{i}=signal{i}+w{i};
end
end