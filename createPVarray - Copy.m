classdef createPVarray < handle
    
    properties
        modules = {};
        rows = 0;
        columns = 0;
        rowPitch = 0;
        azimuth = 0;
        modulePitch = 0;
    end
    
    methods
        function obj = createPVarray(module,numberOfrows,numberOfcolumns,rowPitch,modulePitch,azimuth, displacement)
            obj.rows = numberOfrows;
            obj.columns = numberOfcolumns;
            obj.rowPitch = rowPitch;
            obj.modulePitch = modulePitch;

            obj.modules = cell(numberOfrows,numberOfcolumns);
            
            row_distance = 2;%offset between rows

            for i=1:numberOfrows
                column_distance = 2;  % offset of panel along row (y vector)
                for j=1:numberOfcolumns
                    panel = createPanel(module,displacement); %replicates module
                    panel.move([row_distance,column_distance,0]); % placing the new module
                    obj.modules(i,j) = {panel}; %saving in cell
                    column_distance=column_distance+module.width+modulePitch;
                end
                
                row_distance = row_distance-rowPitch;
            end
           
        %rotates array (in case there is an azimuth)
        cellfun(@(x) x.rotateFromOrigin(azimuth),obj.modules);
        obj.azimuth = azimuth;
        obj.modules{1,1};
        end
        
     function drawArray(obj)
            %3D representation of PV array
            for i=1:obj.rows
                for j=1:obj.columns
                   obj.modules{i,j}.draw;
                end
            end
            grid on
            axis equal

        end
        
%         function draw(obj,obstacle, vertices,xdat,ydat,sumout)
%             %3D representation of PV array
%             init_time = datetime('24-Jul-2023 07:00:00');
%             counter = 0;
%             w=4069
% %             for w = 4063:4074
% 
%             figure('Position', get(0, 'Screensize'))
%             hold on
%             for i=1:obj.rows
%                 for j=1:obj.columns
%                    obj.modules{i,j}.draw;
%                 end
%             end
% 
% %                 for k = 1:9
% %                     createShadow.drawshadow(k,vertices,w);
% %                     
% %     
% %                 end
% %                   Add a title to the figure with the current time
%             sumout = sumout./max(sumout)
% %             sumout = sumout.*1.3
% 
%             % Ensure your data is in column vectors: xdat, ydat, and sumout
% 
%             % Create a grid from the x and y data
%             [xgrid, ygrid] = meshgrid(unique(xdat), unique(ydat));
%             
%             % Create a matrix for the z data
%             zgrid = griddata(xdat, ydat, sumout, xgrid, ygrid, 'cubic');
%             
%             inverse_summer = flipud(colormap('summer'));
%             num_colors = 64;  % Number of colors in the colormap
%             green_colormap = [linspace(0, 1, num_colors)' linspace(1, 0, num_colors)' linspace(0, 1, num_colors)'];
%             
%             % Create the contour plot on the XY plane at z = 0
%             contourf(xgrid, ygrid, zgrid);
%             colormap(inverse_summer);
%             cb = colorbar;
%             title(cb, 'Relative Growth', 'FontSize', 28);
% 
%             stem3(xdat,ydat,sumout,'Color',[0 0.5 0])
% 
% 
% %             colorbar;
%             % Create the contour plot on the XY plane at z = 0
% %             [C, h] = contour(xgrid, ygrid, zgrid, 'LineColor', 'green');
% %             h.ZData = zeros(size(h.ZData));  % Set ZData to zero for all points
%             % Create the contour plot with a green colormap
% % %             figure;
% %             contourf(xgrid, ygrid, zgrid);
% %             colormap('green');
% %             colorbar;
% 
%             view([-37.5 70])
% %             view(90,0)
% %             view(0,0)
%             curr_time = init_time + hours(counter);
%             titlestr = datestr(curr_time);
%             title(titlestr);
% 
% 
% %             obstacle.draw;
% %                 view(0,90);
%                % Save a PNG image of the figure with a unique filename
%             filename = sprintf('goat2growthfigurefigure_%d.png', w);
%             print(filename, '-dpng');
% 
% %                 FigH = figure('Position', get(0, 'Screensize'));
% %                 saveas(FigH, 'Foos.png','png');     
% 
% %                 print(filename, '-dpng');
%             % Resize the saved PNG image to match the size of the figure
% %                 fig = gcf;
% %                 pos = fig.Position;
% %                 set(fig, 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 pos(3:4)]);
% %                
% %                % Save the resized figure as a PNG image
% %                 im = imread(filename);
% %                 imwrite(im, filename, 'png');
% 
%             hold off
%             clf;
%             counter = counter+1;
% %             end
        %         end
        function obj = plotcrops(obj, all_panels, all_beams, all_pylons, x_points, y_points)
            % Plot the original panel and the modified panels
            figure;
            hold on;
            for i = 1:length(all_panels)
                patch('Vertices', all_panels{i}.vertices, ...
                      'Faces', all_panels{i}.faces, ...
                      'FaceColor', [1.0, 0.7, 0.0], ...
                      'EdgeColor', [0.1, 0.1, 0.1], ...
                      'FaceAlpha', 0.7);
            end
            % Plot the beams
            for i = 1:length(all_beams)
                patch('Vertices', all_beams{i}.vertices, ...
                      'Faces', all_beams{i}.faces, ...
                      'FaceColor', [0.7, 0.0, 1.0], ...
                      'EdgeColor', [0.1, 0.1, 0.1], ...
                      'FaceAlpha', 0.7);
            end
        
            % Plot the pylons
            for i = 1:length(all_pylons)
                patch('Vertices', all_pylons{i}.vertices, ...
                      'Faces', all_pylons{i}.faces, ...
                      'FaceColor', [0.0, 0.7, 1.0], ...
                      'EdgeColor', [0.1, 0.1, 0.1], ...
                      'FaceAlpha', 0.7);
            end
        
            % Add labels and title
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('Original Panel, Modified Panels, and Modified Beams');
            axis equal;
            grid on;
            legend('Original Panel', 'Modified Panels', 'Modified Beam 1', 'Modified Beam 2');
        
            % Now plot the extra bit
            init_time = datetime('24-Jul-2023 07:00:00');
            counter = 0;
            for w = 4063:4074
                figure('Position', get(0, 'Screensize'))
                hold on
                for i = 1:obj.rows
                    for j = 1:obj.columns
                        obj.modules{i,j}.draw;
                    end
                end
        
                % Plot the sum output
                sumout = obj.sum_outputs(:,w) ./ max(obj.sum_outputs(:,w));
                stem3(x_points, y_points, sumout, 'Color', [0 0.5 0]);
        
                % Create a grid from the x and y data
                [xgrid, ygrid] = meshgrid(unique(x_points), unique(y_points));
                
                % Create a matrix for the z data
                zgrid = griddata(x_points, y_points, sumout, xgrid, ygrid, 'cubic');
                
                % Create the contour plot
                contourf(xgrid, ygrid, zgrid);
                colormap(flipud(colormap('summer')));
                cb = colorbar;
                title(cb, 'Relative Growth', 'FontSize', 28);
        
                view([-37.5 70]);
                curr_time = init_time + hours(counter);
                titlestr = datestr(curr_time);
                title(titlestr);
        
                % Save the figure
                filename = sprintf('goat2growthfigurefigure_%d.png', w);
                print(filename, '-dpng');
        
                hold off
                clf;
                counter = counter + 1;
            end

end
        
     
    end
    
end

%calculate diffuse B = cellfun(@(x) x.diffuse,A);


