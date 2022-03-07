

for i=1:30
    if encoursdetraitement(sprintf('perf_%d.mat',i))
        [lessalves, ~, labelInit] = test_generation;
        perf = testARPHMM_salves_silmulees(lessalves, labelInit);
        sauverresultatstestARPHMMsimu(sprintf('perf_%d',i),perf);
    else
        disp('Fichier existant, passe au suivant...')
    end
end

