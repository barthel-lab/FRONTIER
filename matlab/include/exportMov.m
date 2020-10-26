function exportMov(filename, title)
    
    % Axis equalize
    axis equal
    
    % Change Axis Font for clarity; add labels
    set(gca,'linewidth',2)
    set(gca,'fontsize',14)
    set(get(gca,'XLabel'),'string','X (mm)')
    set(get(gca,'YLabel'),'string','Y (mm)')
    set(get(gca,'ZLabel'),'string','Z (mm)')
    
    % Set title
    set(get(gca,'title'),'string', title)
    
    % Set background ('w' = white, 'none' = transparent)
    set(gcf,'color','w')
    
    vidfp = strcat('../results/mov/',filename)
    
    % Save video (3d lineage + tumor volume + brain volume)
    OptionZ.FrameRate=5;OptionZ.Duration=15;OptionZ.Periodic=true; % Rendering settings
    CaptureFigVid([-20,10;-110,35;-190,80;-290,10;-380,10],vidfp,OptionZ)

end