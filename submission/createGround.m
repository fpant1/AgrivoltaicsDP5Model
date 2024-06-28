classdef createGround < handle

%you pass:
%     % Meshes an area of ground into squares
%     % INPUTS:
%     %   xmin: minimum x-coordinate of the ground
%     %   xmax: maximum x-coordinate of the ground
%     %   ymin: minimum y-coordinate of the ground
%     %   ymax: maximum y-coordinate of the ground
%     %   nx: number of squares in the x-direction
%     %   ny: number of squares in the y-direction
%     % OUTPUTS:
%     %   vertices: an (nx+1)*(ny+1) by 2 matrix of vertex coordinates
%     %   squares: an nx*ny by 4 matrix of square indices
%     

properties
    xmin;
    xmax;
    ymin;
    ymax;
    nx;
    ny;
    vertices
    X
    Y
    points

end
  

methods 
    
    % constructor
    function obj = createGround(xmin, xmax, ymin, ymax, nx, ny)

        obj.xmin = xmin;
        obj.xmax = xmax;
        obj.ymin = ymin;
        obj.ymax = ymax;
        obj.nx = nx;
        obj.ny = ny;
        
           % Create a vector of x and y values using the "linspace" function
        x_values = linspace(xmin, xmax, nx);
        y_values = linspace(ymin, ymax, ny);
        
        % Use the "meshgrid" function to create a grid of all possible (x,y) pairs
        [x_grid, y_grid] = meshgrid(x_values, y_values);
        
        % Reshape the x and y grids into column vectors
        x_points = reshape(x_grid, [], 1);
        y_points = reshape(y_grid, [], 1);
        
        % Combine the x and y points into a single array of 2D points
        points = [x_points, y_points];
        
        obj.points = points;
        % Plot the array of 2D points
%         scatter(points(:,1), points(:,2));
%         xlabel('X');
%         ylabel('Y');
%         title('2D points between X and Y max/min values');
        
    end
    
    function draw(obj)
        % Plot the mesh
       
%         patch('Faces', obj.squares, 'Vertices', obj.vertices, 'FaceColor', 'none', 'EdgeColor', 'k');
        plot3(obj.points(:,1), obj.points(:,2), zeros(length(obj.points)),'MarkerSize',10,'Marker','x','LineStyle','none','Color',[1 0 0] );
        
        ax = gca;
        ax.FontSize = 28;

        xlabel('x - North - South','FontSize',28);
        ylabel('y - East - West','FontSize',28);

        axis equal;
    end

 

end

end


    
   
