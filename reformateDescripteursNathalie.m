clear all
clear system

system.MD=200e5;                 % Maximal duration [uS]
system.HDT=40;                % Hit definition time [uS]
system.HLT=100;                % Hit lockout time [uS]
system.PDT=20;                % Peak definition time [uS]

system.TH=35;                 % threshold level in the same units of the waveform.

system.wavelet='db45';
system.lev=14;
system.SORH='s';
system.scal='mln';
system.TPTR='sqtwolog';
system.maxFreq = nan;
system.salves=1;
system.Vref=1e-6;


dirr='/home/emmanuel.ramasso/Documents/OndesStageFEMPTO/MinicompositeSiC/MC1/WF/';
system.SR=20e6/1000

% dirr='/home/emmanuel.ramasso/Documents/OndesStageFEMPTO/MinicompositeSiC/MC1/WF_locs/'
% system.SR=20e6/1000

% dirr='/home/emmanuel.ramasso/Documents/OndesStageFEMPTO/MinicompositeSiC/MC2/WF_FLT8/';
% system.SR=2e6/1000

% dirr='/home/emmanuel.ramasso/Documents/OndesStageFEMPTO/MinicompositeSiC/MC2/WF_FLT8_locs/';
% system.SR=2e6/1000

lesfichiers=dir([dirr '*.csv'])


% question : les amplitudes saturent à 100dB chez toi, 120 chez moi ??


preampParVoie=40;
fact = 10^(preampParVoie / 20);

[~,meanGroupDelay] = denoising(rand(1000,1) , system.wavelet, system.lev, system.SORH, system.scal, system.TPTR, true);
meanGroupDelay = round(meanGroupDelay);

NFFT=1024;
freq = system.SR/2*linspace(0,1,NFFT/2+1);


desc=[]; lesSalves=[];
for i=1:length(lesfichiers)
    d=csvread([dirr lesfichiers(i).name],12,0);
    d=d(:,2);
    %diff(d(:,1))
    f=fopen([dirr lesfichiers(i).name],'rt');
    for j=1:9, g=fgets(f); end
    l=length('CHANNEL NUMBER: ');
    chan=str2double(g(l:end));
    g=fgets(f);
    g=fgets(f);
    l=length('TIME OF TEST: ');
    tps=str2double(g(l:end));
    
    d=d/fact;% preamp cancellation
    d=d/system.Vref;% in mV
    
    deb = denoising(d, system.wavelet, system.lev, system.SORH, system.scal, system.TPTR, false);
    deb(find(abs(deb)<10^((system.TH)/20)))=0;
    
    deb(1:meanGroupDelay) = [];
    
    [H,salvesTmp] = hitDetection(deb, system, freq, ...
        tps, chan, d(1:end-meanGroupDelay));
    
    if isempty(H),
        figure(1),clf,plot(d), if input('Valider H=[] ??'), else, error('??'), end%H=zeros(23,1);%error('??')
    end
    
    %     figure,plot(salvesTmp{1})
    %
    %     %d=d/1000; % en mv
    %     try
    %         [D nomsDesColonnes] = features_extraction_salves(d, tps, 2E6, 1e-3, chan, 1e-6, 40);
    %     catch
    %         D=zeros(1,30); D(1)=tps;
    %     end
    
    desc=[desc H];
    lesSalves=[lesSalves ; salvesTmp];
    if mod(i,round(10/100*length(lesfichiers)))==0, fprintf('%d/%d\n',i,length(lesfichiers)); end
    
end

disp(sprintf('A detecté %d salves vs %d files',size(desc,2),length(lesfichiers)))

% % retri a cause des fichiers lus dans le desordre
% [a b]=sort(desc(:,1),'ascend');
% desc=desc(b,:);

[carac,T,ordre] = reordonneTri(desc);
salvesOut = lesSalves(ordre);

% save caracMC1_NG carac system lesSalves
% save caracMC1_loc_NG carac system lesSalves
% save caracMC2_NG carac system lesSalves
% save caracMC2_loc_NG carac system lesSalves

% save 20110622_MiniComp2_descripteursNew D nomsDesColonnes
% save 20110622_MiniComp2_descripteursNew D nomsDesColonnes


figure,plot(carac(:,5),carac(:,13),'*')
