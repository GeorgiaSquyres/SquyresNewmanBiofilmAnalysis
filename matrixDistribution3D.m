% 3D Model: Coupling of lysis to biofilm volume
% This turns out to be identical to the 1D case

% Georgia Squyres, Newman Lab, Caltech

function matrixDistribution3D

radius = 1:400; % specify linearly increasing biofilm radius
% Here, radius units are 0.1 um
time = 1:length(radius); % time steps to model

% Patterned model:
lysisPDF = normpdf(0:max(radius),50,30); % create lysis PDF, gaussian with mu=5, sigma=3
matrixPatterned = zeros(length(radius),1);
layerVolume = zeros(length(radius),1);
for t = time
    r = radius(t):-1:0;
    % Weight volume of each layer by lysis probability in that layer
    for i = 1:length(r)-1 
        currVol = ((4/3)*pi.*(r(i).^3))./2-((4/3)*pi.*(r(i+1).^3))./2;
        currHeight = r(i);
        matrixPatterned(currHeight) = matrixPatterned(currHeight)+currVol.*lysisPDF(i);
        layerVolume(currHeight) = currVol;
     end
end
% Finally, divide by volume at each layer to get a concentration. Notice
% that we multipled by layer volume above, and divide here: hence, this
% is identical to the 1D case and we need only consider the cumulative
% probability distribution. 
matrixPatterned = matrixPatterned./layerVolume;
matrixPatterned = matrixPatterned./mean(matrixPatterned); % normalize 

% Uniform model:
lysisPDF = ones(length([0:1:max(radius)]),1); % create uniform lysis PDF
matrixUniform = zeros(length(radius),1);
layerVolume = zeros(length(radius),1);
for t = time
    r = radius(t):-1:0;
    % Weight volume of each layer by lysis probability in that layer
    for i = 1:length(r)-1 
        currVol = ((4/3)*pi.*(r(i).^3))./2-((4/3)*pi.*(r(i+1).^3))./2;
        currHeight = r(i);
        matrixUniform(currHeight) = matrixUniform(currHeight)+currVol.*lysisPDF(i);
        layerVolume(currHeight) = currVol;
     end
end
% Finally, divide by volume at each layer to get a concentration
matrixUniform = matrixUniform./layerVolume;
matrixUniform = matrixUniform./mean(matrixUniform); % normalize

figure; hold on;
plot(linspace(0,40,length(matrixPatterned)),flipud(matrixPatterned),'LineWidth',2)
plot(linspace(0,40,length(matrixUniform)),flipud(matrixUniform),'LineWidth',2)

xlabel('Biofilm depth')
ylabel('[eDNA] (normalized)')
set(gca,'LineWidth',2,'FontSize',18,'TickDir','out')

temp = legend({'Patterned model','Uniform model'}); temp.Location = 'northwest'; temp.Box = 'off';

%}