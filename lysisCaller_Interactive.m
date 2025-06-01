% User interface for manual annotation of lysis events

% Georgia Squyres, Newman Lab, Caltech

function pointsList = lysisCaller_Interactive(img)

pointsList = [];

sizeZ = size(img,3); sizeT = size(img,4);

%pointsList = [];

disp('Initializing...');

maxScale = zeros(sizeT);
for i = 1:size(img,4)
    currImg = img(:,:,:,i);
    maxScale(i) = prctile(currImg(:),99);
end
%}

f = figure('Position',[342.5000 51 755.5000 740],'KeyPressFcn',@keyPress_Callback); 
showAllPts = 0;
hold on;

% Z scroll GUI elements
sldZ = uicontrol('Style', 'slider','Min',1,'Max',sizeZ,'Value',1,'SliderStep',[1/sizeZ,10/sizeZ], ...
    'Units','pixels','Position', [77.5 70 600 20], 'Callback',@sldZ_Callback);
textZ = uicontrol('Style', 'text','String',['Z 1/',num2str(sizeZ)],'FontSize',12,'HorizontalAlignment','left', ...
    'Units','pixels','Position', [20 72 60 20]);
addlistener(sldZ,'ContinuousValueChange',@sldZ_Callback);

% T control GUI elements
sldT = uicontrol('Style', 'slider','Min',1,'Max',sizeT,'Value',1,'SliderStep',[1/sizeT,10/sizeT], ...
    'Units','pixels','Position', [77.5 30 600 20], 'Callback',@sldT_Callback);
textT = uicontrol('Style', 'text','String',['T 1/',num2str(sizeT)],'FontSize',15,'HorizontalAlignment','left', ...
    'Units','pixels','Position', [20 32 60 20]);
addlistener(sldT,'ContinuousValueChange',@sldT_Callback);

allPtsButton = uicontrol('Style','radiobutton','String','Show All','FontSize',15,'HorizontalAlignment','left', ...
    'Units','pixels','Position', [530 705 100 20],'Callback',@allPts_Callback);

% display first figure
currZ = 1; currT = 1;
currImg = img(:,:,currZ,currT);
h = imshow(currImg,[0 maxScale(currT)],'InitialMagnification',250);
set(h,'ButtonDownFcn',@click_Callback);
h2 = plot(0,0,'.','Visible','off');

uiwait(f); % delay return until figure is closed 

%% GUI callbacks

% Z slider
function sldZ_Callback(varargin)
    % show current image
	currZ = round(get(sldZ,'Value'));
    set(textZ,'String',['Z ',num2str(currZ),'/',num2str(sizeZ)])
    updateFigure;
end

% T slider
function sldT_Callback(varargin)
    % show current image
	currT = round(get(sldT,'Value'));
    set(textT,'String',['T ',num2str(currT),'/',num2str(sizeT)])
    updateFigure;
end

% Arrow key callback 
function keyPress_Callback(~,event)
    if strcmp(event.Key,'rightarrow')
        if currT < sizeT
            currT = currT + 1;     
            set(sldT,'Value',currT)
            set(textT,'String',['T ',num2str(currT),'/',num2str(sizeT)])
            updateFigure;
        end
    elseif strcmp(event.Key,'leftarrow')
        if currT > 1
            currT = currT - 1;
            set(sldT,'Value',currT)
            set(textT,'String',['T ',num2str(currT),'/',num2str(sizeT)])
            updateFigure;
        end
    elseif strcmp(event.Key,'uparrow')
        if currZ < sizeZ
            currZ = currZ + 1;
            set(sldZ,'Value',currZ)
            set(textZ,'String',['Z ',num2str(currZ),'/',num2str(sizeZ)])
            updateFigure;
        end
    elseif strcmp(event.Key,'downarrow')
        if currZ > 1
            currZ = currZ - 1;
            set(sldZ,'Value',currZ)
            set(textZ,'String',['Z ',num2str(currZ),'/',num2str(sizeZ)])
            updateFigure;
        end
    end
end

% Click to add points
function click_Callback(~,event)
    coords = event.IntersectionPoint; % pixel coordinates of intersection
    pointsList = [pointsList;coords(1) coords(2) currZ currT];
    updateFigure;
end

% Show/hide all points
function allPts_Callback(varargin)
    if get(allPtsButton,'Value') == 0
        showAllPts = 0;
    else
        showAllPts = 1;
    end
    updateFigure
end

% Redraw figure and points after any GUI update
function updateFigure
    currImg = img(:,:,currZ,currT);
    delete(h);
    h = imshow(currImg,[0 maxScale(currT)],'InitialMagnification',250);
    set(h,'ButtonDownFcn',@click_Callback);
    hold on;
    if ~isempty(pointsList)
        delete(h2);
        if showAllPts
            nearbyPoints = and(pointsList(:,3) > currZ-5,...
                pointsList(:,3)<currZ+5);
            h2 = plot(pointsList(nearbyPoints,1),pointsList(nearbyPoints,2),'.r','MarkerSize',20);
        else
            nearbyPoints = and(...
                and((pointsList(:,3) > currZ-5),...
                (pointsList(:,3)<currZ+5)),...
                and((pointsList(:,4)>currT-2),...
                (pointsList(:,4)<currT+2)));
            h2 = plot(pointsList(nearbyPoints,1),pointsList(nearbyPoints,2),'.r','MarkerSize',20);
        end
    end
end

end