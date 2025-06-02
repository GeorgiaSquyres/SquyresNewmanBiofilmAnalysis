% 1D Model: eDNA distribution

% Georgia Squyres, Newman Lab, Caltech

function matrixDistribution1D

depth = 0:100;

matrixPatterned = (1/2)*(1+erf((depth-5)/(3*sqrt(2))));
matrixUniform = depth./max(depth);

figure; hold on;
plot(depth,matrixPatterned,'LineWidth',2);
plot(depth,matrixUniform,'LineWidth',2);

xlabel('Depth')
ylabel('[eDNA] (normalized)')
set(gca,'LineWidth',2,'FontSize',18,'TickDir','out')

temp = legend({'Patterned model','Uniform model'}); temp.Location = 'southeast'; temp.Box = 'off';