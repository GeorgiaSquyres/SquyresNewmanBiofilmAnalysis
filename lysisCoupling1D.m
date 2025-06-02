% 1D Model: Coupling of lysis to biofilm volume

% Georgia Squyres, Newman Lab, Caltech

function lysisCoupling1D

height = 1:100; % specify linearly increasing biofilm height
time = 1:length(height); % time steps to model

% Patterned model:
lysisPatterned = ones(length(height),1); 
lysisPatterned = cumsum(lysisPatterned);
lysisPatterned = lysisPatterned./max(lysisPatterned);

% Uniform model:
lysisUniform = height;
lysisUniform = cumsum(lysisUniform);
lysisUniform = lysisUniform./max(lysisUniform);

figure; hold on;
plot(height,lysisPatterned,'LineWidth',2);
plot(height,lysisUniform,'LineWidth',2);

xlabel('Biofilm height')
ylabel('Cumulative lysed cells (normalized)')
set(gca,'LineWidth',2,'FontSize',18,'TickDir','out')

temp = legend({'Patterned model','Uniform model'}); temp.Location = 'northwest'; temp.Box = 'off';

end