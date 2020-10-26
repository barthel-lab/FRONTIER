function plotPhylo3D(nodes,seg)
    
    % Plot segments connecting nodes
    for i = 1:size(seg,1)
        p1 = nodes(nodes.Node==seg.Node1(i),:)
        p2 = nodes(nodes.Node==seg.Node2(i),:)
        plot3([p1.X; p2.X], [p1.Y; p2.Y], [p1.Z; p2.Z], '--black', 'MarkerFaceColor', 'black', 'MarkerSize', 15, 'LineWidth', 1)
        %arrow([p1.X p1.Y p1.Z], [p2.X p2.Y p2.Z], 12, 'BaseAngle', 60, 'LineStyle', '--', 'MarkerSize', 1, 'LineWidth', 1) %, 'LineStyle', '--', 'facealpha', 0.5, 'color', 'black', 'stemWidth', 0.4) %'--black', 'MarkerFaceColor', 'black', 'MarkerSize', 15, 'LineWidth', 1)
    end
    
    % Plot nodes in 3d space
    for i = 1:size(nodes,1)
        x = nodes.X(i)
        y = nodes.Y(i)
        z = nodes.Z(i)
        label = nodes.Biopsy(i)
        plot3(x, y, z, 'o', 'MarkerFaceColor', findColor(nodes.Subtype(i)), 'MarkerEdgeColor', 'black', 'MarkerSize', findSize(string(nodes.Tip(i))), 'LineWidth', 2.5)
        if(strcmp(label,'NA') == 0)
            text(x+0.3, y+0.3, z+0.3, label)
        end
    end
end