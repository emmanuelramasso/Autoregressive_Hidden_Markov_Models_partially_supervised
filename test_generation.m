function [lessalves, lesentrees, labelInit] = test_generation

clear lessalves lesentrees entree salve energies
zeta = [0.02 0.05 0.2 0.35]
omega = [75 50 30 20]
Nsim = 1000
Fe = 200

for i=1:length(zeta)
    k=1;
    for j=1:Nsim 
        [salve{i}(k,:), entree{i}(k,:)] = genererSalves(omega(i), zeta(i), Fe);
        energies(k,i) = rms(salve{i}(k,:));
        k=k+1;
        if mod(j,round(10/100*Nsim))==0
            fprintf('\nz=%d/%d, sim=%d/%d',i,length(zeta),j,Nsim);
        end
    end
end

q=quantile(energies,0.10);
labelInit=[];
for i=1:length(zeta)
    f = find(energies(:,i)<=q(i));
    lessalves{i} = salve{i}(f,:);
    lesentrees{i} = entree{i}(f,:);
    labelInit = [labelInit ; ones(size(lessalves{i},1),1)*i];
end

return 






% figure,
% for i=1:length(zeta)
%     for j=1:size(lessalves{1},1)
%         plot(lessalves{i}(j,:)), pause(0.25), clf,
%     end
%     disp('attend')
%     pause
% end

% on utilise ensuite 10% des salves pour apprendre et 90% pour tester

system.MD=200e5;                 % Maximal duration [uS]
system.HDT=100*1000;                % Hit definition time [uS]
system.HLT=500*1000;                % Hit lockout time [uS]
system.PDT=20*1000;                % Peak definition time [uS]

system.TH=20;                 % threshold level in the same units of the waveform.

system.wavelet='db45';
system.lev=5;
system.SORH='s';
system.scal='mln';
system.TPTR='sqtwolog';
system.maxFreq = nan;
system.salves=1;
system.Vref=1e-6;


system.SR=Fe/1000

meanGroupDelay = 0;%round(meanGroupDelay);

NFFT=1024;
freq = system.SR/2*linspace(0,1,NFFT/2+1);

desc=[]; lesSalves=[]; label=[]; n=0;
for i=1:length(zeta)
    
    d2=[]; s2={};
    for k=1:size(lessalves{i},1)
        
        d=lessalves{i}(k,:);
        d=d/max(d);
        d=d*1000;
        
        deb = denoising(d, system.wavelet, system.lev, system.SORH, system.scal, system.TPTR, false);
        deb(find(abs(deb)<10^((system.TH)/20)))=0;
        
        deb(1:meanGroupDelay) = [];
        
        [H,salvesTmp] = hitDetection(deb, system, freq, ...
            0, 1, d(1:end-meanGroupDelay));
        
        d2=[d2 H];
        s2=[s2 ; salvesTmp];
        
    end
        
    disp(sprintf('A detecté %d salves vs %d salves reelles',size(d2,2),size(lessalves{i},1)));

    label=[label ; i*ones(size(d2,2),1)];

    desc = [desc d2]; 
    lesSalves = [lesSalves ; s2];    
    
end


[carac,T,ordre] = reordonneTri(desc);
    
lesSalves = lesSalves(ordre);
label = label(ordre);

plotmatrix_mine(carac(:,[13 12]), label, 1,0);


    