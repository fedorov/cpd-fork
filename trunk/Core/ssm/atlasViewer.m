function atlasViewer(atlas)

% if(~exist('atlas.tri'))
%     atlas.tri = 1;
% end

%% Use system background color for GUI components
panelColor = get(0,'DefaultUicontrolBackgroundColor');

%% ------------ Callback Functions ---------------
% [FileName,PathName,FilterIndex] = uigetfile();
% res = load([PathName '\' FileName]);
% atlas = res.atlas;
% atlas = res.nonrigidAtlas;
% atlas = res.affineAtlas;
% res.mean = res.mean*res.scale;
% res.mods = res.mods*res.scale;
% res.meanVar = res.meanVar*res.scale;
% res.latent = res.latent*res.scale;
range = 3*(atlas.latent.^0.5);

% m_t = myCrustOpen(atlas.mean);
% Figure resize function
function figResize(src,evt)
    fpos = get(f,'Position');
    set(botPanel,'Position',...
        [1/20 1/20 fpos(3)-.1 fpos(4)*8/35])
    set(rightPanel,'Position',...
        [fpos(3)*85/120 fpos(4)*8/35 fpos(3)*35/120 fpos(4)*27/35])
    set(centerPanel,'Position',...
        [1/20 fpos(4)*8/35 fpos(3)*85/120 fpos(4)*27/35]);
end
% Bottom panel resize function
function botPanelResize(src, evt)
    bpos = get(botPanel,'Position');
    set(plotButton,'Position',...
        [bpos(3)*10/120 bpos(4)*2/8 bpos(3)*24/120 2])
    set(varPanel,'Position',...
        [bpos(3)*80/120 bpos(4)*1/8 bpos(3)*20/120 bpos(4)*6/8])
end

% Right panel resize function
function rightPanelResize(src,evt)
    rpos = get(rightPanel,'Position');
    set(firstModeSlider,'Position',...
        [rpos(3)*4/32 rpos(4)*16/27 rpos(3)*24/32 rpos(4)*2/27]);
    set(secondModeSlider,'Position',...
        [rpos(3)*4/32 rpos(4)*12/27 rpos(3)*24/32 rpos(4)*2/27]);
    set(thirdModeSlider,'Position',...
        [rpos(3)*4/32 rpos(4)*8/27 rpos(3)*24/32 rpos(4)*2/27]);
    set(smoothRadioButton,'Position',...
        [rpos(3)*4/32 rpos(4)*24/27 150 10]);
    set(varMag1RadioButton,'Position',...
        [rpos(3)*4/32 rpos(4)*8/27 150 10]);
    set(varMag2RadioButton,'Position',...
        [rpos(3)*40/32 rpos(4)*8/27 150 10]);
end
%% Callback for plot button
function plotButtonCallback(src,evt)
    plotShape();
end
function firstModeSliderCallback(src, evt)
    plotShape();
end
function secondModeSliderCallback(src, evt)
    plotShape();
end
function thirdModeSliderCallback(src, evt)
    plotShape();
end
function smoothRadioButtonCallback(src, evt)
    plotShape();
end       
function varMag1RadioButtonCallback(src, evt)
    plotShape();
end       
function varMag2RadioButtonCallback(src, evt)
    plotShape();
end       
function openMenuCallback(src, evt)
    [FileName,PathName,FilterIndex] = uigetfile();
    res = load([PathName '\' FileName]);
end       

%% ------------ GUI layout ---------------
%% Set up the figure and defaults
f = figure('Units','characters',...  
        ... %'Position',[30 30 120 35],... EDIT: Sanchez,(puts offscreen)
        'Color',panelColor,...
        'HandleVisibility','callback',...
        'IntegerHandle','off',...
        'Renderer','painters',...
        'Toolbar','figure',...
        'NumberTitle','off',...
        'Name','Workspace Plotter',...
        'ResizeFcn',@figResize);
%% Create the main menu   
men = uimenu('Label','Workspace', 'parent', f);
openMenu = uimenu(men,'Label','Open Atlas','Callback',@openMenuCallback);    
%% Create the main uipanel    
centerPanel = uipanel('bordertype','etchedin',...
    'Units','characters',...
    'Position', [17/20 8 88 27],...
    'Parent',f);
%% Create the right side panel
rightPanel = uipanel('bordertype','etchedin',...
    'BackgroundColor',panelColor,...
    'Units','characters',...
    'Position',[88 8 32 27],...
    'Parent',f,...
    'ResizeFcn',@rightPanelResize);
%% Create the bottom uipanel
botPanel = uipanel('BorderType','etchedin',...
    'BackgroundColor',panelColor,...
    'Units','characters',...
    'Position',[1/20 1/20 119.9 8],...
    'Parent',f,...
    'ResizeFcn',@botPanelResize);
%% Create the bottom uipanel
varPanel = uipanel('BorderType','etchedin',...
    'BackgroundColor',panelColor,...
    'Units','characters',...
    'Position',[80 1/20 20 20],...
    'Parent',botPanel);
%% Add an axes to the center panel
a = axes('parent',centerPanel);
b = axes('parent',varPanel);
%% Add plot buttons      
plotButton = uicontrol(f,'Style','pushbutton','Units','characters',...
        'Position',[10 2 24 2],...
        'String','Create Plot',...
        'Parent',botPanel,...
        'Callback',@plotButtonCallback);
%% Add first mode slider    
firstModeSlider = uicontrol(f,'Style','slider','Units','characters',...
        'Position',[2 40 24 2],...
        'String','Create Plot',...
        'Parent',rightPanel,...
        'Callback',@firstModeSliderCallback,...
        'Min', -range(1), 'Max', range(1));
    
%% Add second mode slider    
secondModeSlider = uicontrol(f,'Style','slider','Units','characters',...
        'Position',[2 40 24 2],...
        'String','Create Plot',...
        'Parent',rightPanel,...
        'Callback',@secondModeSliderCallback,...
        'Min', -range(2), 'Max', range(2));   
%% Add third mode slider    
thirdModeSlider = uicontrol(f,'Style','slider','Units','characters',...
        'Position',[2 40 24 2],...
        'String','Create Plot',...
        'Parent',rightPanel,...
        'Callback',@thirdModeSliderCallback,...
        'Min', -range(3), 'Max', range(3)); 
%% Add smoothing radio button
smoothRadioButton = uicontrol('Style','checkbox','String','Surface Smoothing',...
    'Position',[2 40 20 2],'parent',rightPanel,'HandleVisibility','off',...
    'Callback',@smoothRadioButtonCallback);
%% Add button to show variation magnitude
varMag1RadioButton = uicontrol('Style','checkbox','String','1st mode Magnitude',...
    'Position',[2 40 24 2],'parent',rightPanel,'HandleVisibility','off',...
    'Callback',@varMag1RadioButtonCallback);
varMag2RadioButton = uicontrol('Style','checkbox','String','2nd mode Magnitude',...
    'Position',[2 40 24 2],'parent',rightPanel,'HandleVisibility','off',...
    'Callback',@varMag2RadioButtonCallback);

%% A function for plotting data    
function plotShape()
    axes(a);
    s = zeros(length(atlas.latent), 1);
    s(1) = get(firstModeSlider,'Value');
    s(2) = get(secondModeSlider,'Value');
    s(3) = get(thirdModeSlider,'Value');
    ch = atlas.mods*s;
    ad = reshape(ch, 3, length(ch)/3)';
    
    newP = atlas.mean+ad;
    if(get(smoothRadioButton, 'Value'))
        FV.vertices = newP;
        FV.faces = m_t;
        FV2=smoothpatch(FV,0,1);
        m_t = FV2.faces;
        newP = FV2.vertices;
    end
    if(get(varMag1RadioButton, 'Value'))
        mag = reshape(atlas.mods(:,1), 3, length(atlas.mods)/3)';
        mag = sum(mag.^2, 2).^0.5;
        cla
        p1 = patch('Faces',m_t,'Vertices',newP);
        set(p1,'FaceColor','flat','FaceVertexCData',mag);
    elseif(get(varMag2RadioButton, 'Value'))
        mag = reshape(atlas.mods(:,2), 3, length(atlas.mods)/3)';
        mag = sum(mag.^2, 2).^0.5;
        cla
        p1 = patch('Ffaces',m_t,'Vertices',newP);
        set(p1,'FaceColor','flat','FaceVertexCData',mag);
    else
%         if(atlas.tri == 1)
%             trisurf(m_t, newP(:,1), newP(:,2), newP(:,3), 'facecolor',[0.9 0.9 0.9]);
%             trisurf(m_t, newP(:,1), newP(:,2), newP(:,3), 'facecolor',[1 1 1]);

             scatter3(newP(:,1), newP(:,2), newP(:,3), 'r');
%             view([-40 90]);
%             view([-180 80]);

        view([180 -90]);
%         view([180 -88]);
%             view([-180 0]);
%  view([90 60]);
%         else
%             scatter3(newP(:,1), newP(:,2), newP(:,3), 'r');
% %             view([-92 16]);
%         end
            
        grid off;
    end
    scale = max(atlas.mean(:));
    axis(scale*[-1.4 1.4 -1.4 1.4 -1.4 1.4]);
    daspect([1,1,1])
    
%     
%     title(FileName);
    axes(b);
    
    for i=1:length(atlas.latent)
        varSum(i) = sum(atlas.latent(1:i))/sum(atlas.latent);
    end
    plot(varSum);
end
end