d=dir('perf_*.mat'); nbNonclasses = [];
perfindic = [];
for i=1:length(d)
    clear perf
    load(d(i).name)
    perfindic = [perfindic ; [perf.a perf.b perf.c perf.c]];
    nbNonclasses = [nbNonclasses ; 100-sum(perf.r,1)];
end

figure,boxplot(perfindic), title('Indic. de perf')
figure,boxplot(nbNonclasses), title('Non classes par classe')

