% 2D models of cell lysis, for visualization purposes
% These are stocastic simulations according to the patterned and uniform 
% lysis probability distributions

% Georgia Squyres, Newman Lab, Caltech

function viewModels

%% Patterned cell lysis

currHeight = 0:0.5:40;
temp = linspace(0.5,1,length(currHeight))';
temp2 = linspace(0,0.7,length(currHeight))';
colors = [temp2, temp, temp2];
simR = []; simTheta = []; simT = [];
for i = 1:length(currHeight)
    if currHeight(i) < 3
        areaDeathRegion = 0; 
    else
        areaDeathRegion = 2*pi*(currHeight(i)-3);
    end
    numCells = round(areaDeathRegion/30);
    currDepths = randn(numCells,1)+4;
    currZ = currHeight(i)-currDepths;
    currZ(currZ<0) = [];
    currZ(currZ>currHeight(i)) = [];
    simR = [simR;currZ];
    currTheta = rand(length(currZ),1)*pi;
    simTheta = [simTheta;currTheta];
    simT = [simT;repmat(colors(i,:),[length(currZ) 1])];
end

figure('Name','Patterned model'); hold on;
theta = linspace(0,pi,100);
r = currHeight(end);
x = r.*cos(theta); x(end+1) = x(1);
y = r.*sin(theta); y(end+1) = y(1);
fill(x,y,[0.756 0.756 0.756],'LineWidth',2)
colors = parula(length(currHeight));
scatter(simR.*cos(simTheta),simR.*sin(simTheta),100,simT,'o','MarkerFaceColor','flat','MarkerEdgeColor','k');
axis equal;

%% Uniform cell lysis 

currHeight = 0:0.5:40;
temp = linspace(0.5,1,length(currHeight))';
temp2 = linspace(0,0.7,length(currHeight))';
colors = [temp, temp2, temp2];
simR = []; simTheta = []; simT = [];
for i = 1:length(currHeight)
    areaDeathRegion = pi*(currHeight(i)+3).^2;
    numCells = round(areaDeathRegion/400);
    currRadii = sqrt(rand(numCells,1)).*currHeight(i);
    simR = [simR;currRadii];
    currTheta = rand(length(currRadii),1)*pi;
    simTheta = [simTheta;currTheta];
    simT = [simT;repmat(colors(i,:),[length(currRadii) 1])];
end

figure('Name','Uniform model'); hold on;
theta = linspace(0,pi,100);
r = currHeight(end);
x = r.*cos(theta); x(end+1) = x(1);
y = r.*sin(theta); y(end+1) = y(1);
fill(x,y,[0.756 0.756 0.756],'LineWidth',2)
colors = parula(length(currHeight));
scatter(simR.*cos(simTheta),simR.*sin(simTheta),100,simT,'o','MarkerFaceColor','flat','MarkerEdgeColor','k');
axis equal;
