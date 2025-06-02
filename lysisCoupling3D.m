% 3D Model: Coupling of lysis to biofilm volume

% Georgia Squyres, Newman Lab, Caltech

radius = 1:1000; % specify linearly increasing biofilm radius
% Here, radius units are 0.1 um for better precision
time = 1:length(radius); % time steps to model
vol = ((4/3)*pi.*(radius.^3))./2; % hemisphere biofilm volume over time
vol = vol./max(vol); % normalize

% Patterned model:
lysisPDF = normpdf([0:1:max(radius)],5,3); % create lysis PDF, gaussian with mu=5, sigma=3
lysisPatterned = zeros(length(time),1);
for t = time
    r = radius(t):-1:0;
    % Weight volume of each layer by lysis probability in that layer
    for i = 1:length(r)-1
        currVol = ((4/3)*pi.*(r(i).^3))./2-((4/3)*pi.*(r(i+1).^3))./2;
        lysisPatterned(t) = lysisPatterned(t)+currVol.*lysisPDF(i);
    end
end
lysisPatterned = cumsum(lysisPatterned);
lysisPatterned = lysisPatterned./max(lysisPatterned); % normalize

% Uniform model:
lysisPDF = ones(length([0:1:max(radius)]),1); % create uniform lysis PDF
lysisUniform = zeros(length(time),1);
for t = time
    r = radius(t):-1:0;
    % Weight volume of each layer by lysis probability in that layer
    for i = 1:length(r)-1
        currVol = ((4/3)*pi.*(r(i).^3))./2-((4/3)*pi.*(r(i+1).^3))./2;
        lysisUniform(t) = lysisUniform(t)+currVol.*lysisPDF(i);
    end
end
lysisUniform = cumsum(lysisUniform);
lysisUniform = lysisUniform./max(lysisUniform);

figure; hold on;
plot(vol,lysisPatterned,'LineWidth',2);
plot(vol,lysisUniform,'LineWidth',2);

xlabel('Biofilm volume (normalized)')
ylabel('Cumulative lysed cells (normalized)')
set(gca,'LineWidth',2,'FontSize',18,'TickDir','out')

temp = legend({'Patterned model','Uniform model'}); temp.Location = 'northwest'; temp.Box = 'off';
