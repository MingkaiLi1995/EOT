function [handle_extent] = plot_extent(ellipse,line_style, color, line_width)
% plot a elliptical extension
    centroid_position = [ellipse(1) ; ellipse(2)];
    orientation = ellipse(3);
    L = [ellipse(4) ; ellipse(5)];
    R_matrix = [cos(orientation) -sin(orientation); sin(orientation) cos(orientation)]; 
    
    alpha = 0:pi/100:2*pi;
    x_norotated = L(1)*cos(alpha);
    y_norotated = L(2)*sin(alpha);
    rotated = R_matrix* [x_norotated; y_norotated];
    
    xpoints = rotated(1,:) + centroid_position(1);
    ypoints = rotated(2,:) + centroid_position(2);
    handle_extent = plot(xpoints,ypoints,'LineStyle',line_style,'color',color,'LineWidth',line_width);
end