function [c w]= visualiser_resultats_shm_clust_morceaux_2(b,T,nbc,W,tracerFig)

if nargin==4, tracerFig=1; end
if isempty(W), LB=b; 
else
    LB=[];
    d=1;
    f2=fix(10*length(b)/100);
    f=f2;
    P=fix(W/2);
    cont=true;
    while cont
        if f>length(b), f=length(b); cont=false; end
        [L1 N1] = decouvrir_clusters_sequences_shm(b(d:f), fix(mean(W(d:f))), max(1,fix(mean(P(d:f)))));
        if isempty(L1), L1=repmat(nan,f-d+1,1); end
        LB=[LB;L1];
        d=f+1;
        f=f+f2;
    end
end

[c w]=cherche_cumul_etat_shm(T,nbc,LB,tracerFig);
