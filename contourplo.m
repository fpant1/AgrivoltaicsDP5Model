

function contourplo(x,y,z)
%     tri = delaunay(x,y);
%     trisurf(tri,x,y,z);
%     xlabel('X-axis');
%     ylabel('Y-axis');
%     

    x_unique = unique(x);
    y_unique = unique(y);
    heat = tomato_sunny_SD.sum_outputs
    counter = 1;
    for i = 1:length(x_unique)
        for j=1:length(y_unique)
            heat(i,j) = z(counter)
            counter = counter+1

        end
    end

    heatmap(heat)


%     [X, Y] = meshgrid(x_unique, y_unique);
%     Z = griddata(x, y, z, X, Y);
%     figure;
%     contourf(X, Y, Z);
%     title('Filled Contour Plot');
%     xlabel('X-axis');
%     ylabel('Y-axis');


%     stem3(x,y,z)



end
