function [carac,T,ordre] = reordonneTri(Descripteurs)

% Descripteurs: 23*nbHits
% 1 - chanNum
% 2 - hit start [s]
% 3 - hit max [s]
% 4 - hit stop [s]
% 5 - Rise Time [s]
% 6 - Counts
% 7 - PAC-Energy [uVs]
% 8 - Duration [s]
% 9 - Amplitude [dB]
% 10 - Average Frequency [kHz]
% 11 - RMS [V]
% 12 - ASL [dB] ???
% 13 - Counts to Peak
% 14 - Reverberation Freq. [kHz]
% 15 - Initiation Freq. [kHz]
% 16 - Signal Strength [pV.s]
% 17 - Absolute Energy [aJ]
% 18-21 - Partial Power [%]
% 22 - Frequency Centroid [kHz]
% 23 - Peak Freq. [kHz]


if isempty(Descripteurs)
    carac = [];
    T = [];
    ordre = [];
    return
end

Descripteurs = Descripteurs';
thelength = size(Descripteurs,1);

carac = [Descripteurs zeros(thelength,6)];

carac(:,19:28) = Descripteurs(:,14:23);
carac(:,18) = 40*ones(thelength,1); % THR
carac(:,9:17) = Descripteurs(:,5:13);
carac(:,8) = Descripteurs(:,1); % chanNum
carac(:,7) = Descripteurs(:,4); % hit stop
carac(:,6) = Descripteurs(:,3); % hit max
carac(:,5) = Descripteurs(:,2); % hit start
carac(:,1:4) = zeros(thelength,4);

carac(:,29) = sqrt(carac(:,27).*carac(:,28)); % produit racine Pfrq.*Cfrq


T = Descripteurs(:,2);

[T,ordre]=sort(T,'ascend');
carac = carac(ordre,:);

T = T - T(1);
carac(:,7) = carac(:,7) - carac(1,5);
carac(:,6) = carac(:,6) - carac(1,5);
carac(:,5) = carac(:,5) - carac(1,5);




% Carac nbHits*29
% 1 - 0
% 2 - 0
% 3 - 0
% 4 - 0
% 5 - hit start (T dans l'ordre avec T(1) = 0)
% 6 - hit max
% 7 - hit stop
% 8 - chanNum
% 9 - Rise Time [us]
% 10 - Counts
% 11 - PAC-Energy [uVs]
% 12 - Duration [us]
% 13 - Amplitude [dB]
% 14 - Average Frequency [kHz]
% 15 - RMS [V]
% 16 - ASL [dB] ???
% 17 - Counts to Peak
% 18 - THR
% 19 - Reverberation Freq. [kHz]
% 20 - Initiation Freq. [kHz]
% 21 - Signal Strength [pV.s]
% 22 - Absolute Energy [aJ]
% 23-26 - Partial Power [%]
% 27 - Frequency Centroid [kHz]
% 28 - Peak Freq. [kHz]
% 29 - sqrt(Pfreq*CFreq)