function [salve, entree] = genererSalves(omega, zeta, Fe)

sys = tf(1,[1/omega^2 2*zeta/omega 1]);
x0  = [0,0]';

% Subject to a white noise excitation with:
S  = 1;                  % White noise spectral intensity
T  = 30;                 % Duration of the excitation, s
dt = 1/Fe; %0.02;        % Time increment, s
t  = (0:dt:T)';          % discrete time instants
n  = length(t);          % n points ~ number of random variables

mu = 0;
sig = T/7;
y = normpdf(t, mu, sig);
y = y / max(y);
W1similaire = 5;
% W2lissage = 20;

% figure,plot(t, y)
amplitude = 1/10*fliplr(logspace(0,1,n))';

W = zeros(n,1);

sigLog = 1;

% PHASE 1 : CALME PLAT
idt=10; 

% PHASE 2 : ON EXCITE ALEATOIREMENT SANS SE SOUCIER DU PASSE
while idt<=round(10/100*n)
        
    if rand > y(idt)
        
        r = 0 + exp(-0.01*t(idt))*lognrnd(0,sigLog);%randn;
         
        for i=idt:min(idt+W1similaire,n)
            
            W(idt)  = sign(randn)*sqrt(2*pi*S/dt)*r;
                     
            idt=idt+1;
        end
        
    else
        idt=idt+1;
    end
end

% PAHSE 3 : ON S'ASSURE QUE LES DEPLACEMENTS DECROISSENT SUFFISAMMENT
% (RETOUR AU CALME PLAT)
while idt<=n
        
    if rand > y(idt)
        
        r = 0 + exp(-0.1*t(idt))*lognrnd(0,sigLog);%randn;
         
        for i=idt:min(idt+W1similaire,n)
            
            tmp = amplitude(idt)*sign(randn)*sqrt(2*pi*S/dt)*r;
           
            W(idt) = tmp;%sign(tmp)*min(max(abs(W(max(idt-W2lissage,1):idt-1)))*0.99, abs(tmp));
            
            idt=idt+1;
        end
        
    else
        idt=idt+1;
    end
end

% W=abs(W);

salve = lsim(sys, W, t, x0, 'zoh')';
entree = W; 

% figure,subplot(211),plot(salve),title('Salve')
% subplot(212),plot(W),title('Excitation')
      
    
