classdef createTable
    %creates a table with shading coefficient results
    
    properties
        coefficients
        azRange
        elRange
        module
        obstacle
        timetables
        vertices
        x_min
        x_max
        y_min
        y_max
       


    end
    
    methods
        function obj = createTable(module,obstacle,resolution)
            azRange = -180:resolution:180;
            elRange = 0:resolution:90;
            all_tables = struct();  % create structure to save tables
            table = zeros(length(elRange),length(azRange));
            [az,el] = meshgrid(azRange,elRange);
            obstacle
            for a = 1:2
                for i=1:size(el,1)
                    for j=1:size(az,2)
                        [sunX, sunY, sunZ] = sph2cart(-az(i,j)*pi/180,el(i,j)*pi/180,1);%SPH2CART - Transform spherical to Cartesian coordinates.
                        -az(i,j);
                        el(i,j);

                        sun = [sunX sunY sunZ];
                        square = module(4*a-3:4*a, :);
%                         for i = obstacles
%                          
                        sombra = createShadow(square,obstacle,sun); % in create shadow extract the projected points for every combination of az/el at given resolution
                        table(i,j) = sombra.coef(square); % in create shadow, extract the shading coefficient for that setup
%                         all_tables.(['table_',a]) = table;
                        a_str = num2str(a,'%02d');
                        all_tables.(['table_',a_str]) = table;
                         
                           
%                         end
%                         all_tables.(['table_',a_str]) = table;
                    end
                end
%                 AZrange = -180:1:180;
%                 ELrange = 0:1:90;
%                 if resolution~=1
%                             
%                     tables = zeros(length(ELrange),length(AZrange))
%                     % interpolates values
%                     [AZ,EL] = meshgrid(AZrange,ELrange);
%                     tables = interp2(az,el,table,AZ,EL);
%                   
% 
%                 end
               
            end
            AZrange = -180:1:180;
            ELrange = 0:1:90;
            azRange = -180:resolution:180;
            elRange = 0:resolution:90;
            [az,el] = meshgrid(azRange,elRange);
%             sort this so that interpolation doesnt change it from a
%             structure to just a table
            names = fieldnames(all_tables);
            
            if resolution~=1
                for i = 1:length(names)
                    
                    table = all_tables.(names{i});
                    table(isnan(table)) = 0;
            % interpolates values
                    [AZ,EL] = meshgrid(AZrange,ELrange);
                    all_tables.(names{i}) = interp2(az,el,table,AZ,EL);
                    obj.coefficients = all_tables;
                end
            else
                obj.coefficients = all_tables;
            end
%             obj.coefficients = all_tables
    
            
            % fixes negative values
            for i=1:length(AZrange)
                if AZrange(i)<0
                    AZrange(i) = AZrange(i) + 360;
                end
            end
            obj.azRange = AZrange;
            obj.elRange = ELrange;
            obj.module = module;
            obj.obstacle = obstacle;
            
        end
        
        function coeff = coeff(obj,elevation,azimut,a)
            % returns the coefficient for a given elevation and azimut
            % solar angles
            elevation = round(elevation);
            azimut = round(azimut);
            
            el = obj.elRange;
            az = obj.azRange;

            if azimut==360
                azimut=0;
            elseif azimut<0
                azimut = azimut+360;
            end
            
            table_Name = sprintf('table_%02d', a);     
            if elevation>=1
                elevIndex = find(el==elevation);
                azIndex = find(az==azimut);
                coeff = obj.coefficients.(table_Name)(elevIndex,azIndex(1));
            else
                coeff = nan;
            end

        end
        
        function heatPlot(obj)
            %plots the shade coefficient table
            az = -180:1:180;
            el = 0:1:90;

            surf(az,el,obj.coefficients,'EdgeColor','none')
            xlabel('Azimuth');
            ylabel('Elevation');
        end
        
        function show(obj)
            % shows the shading situation
            hold on
            obj.module.draw;
            obj.obstacle.draw;
            grid on
            axis equal
            hold off
            title('ct show')
        end
        
       
%         
    end

    methods(Static)

        function new_modules = createfullModules(num_modules_x, num_modules_y, module_spacing, row_spacing, tracked_geom)
            % Inputs:
            % - num_modules_x: Number of modules in a row
            % - num_modules_y: Number of modules in a column
            % - module_spacing: Spacing between modules in the x-direction
            % - row_spacing: Spacing between rows in the y-direction
            % - tracked_geom: Structure containing fields panels, beams, pylons
        
            % Create cell array to store new modules
            new_modules = cell(num_modules_y, num_modules_x);
        
            % Apply column (x-axis) offsets for panels, beams, and pylons
            new_modules{1,1} = tracked_geom;
            for j = 1:num_modules_x
                % Offset for this module
                x_offset = (j - 1) * module_spacing;
    
                % Create new module structure
                new_module = struct();
    
                % Reference the correct field name for panels, beams, and pylons

                for pan = 1 : length(fieldnames(tracked_geom.panels))
                    panel_field = sprintf('panel_%d', pan);
                
                    original_vertices = tracked_geom.panels.(panel_field);
                    new_vertices = original_vertices;
                    new_vertices(:, 1) = new_vertices(:, 1) + x_offset;
                    new_module.panels.(panel_field) = new_vertices;
                end

                for beam = 1 : length(fieldnames(tracked_geom.beams))
                    beam_field = sprintf('beam_%d', beam);

                    original_vertices = tracked_geom.beams.(beam_field);
                    new_vertices = original_vertices;
                    new_vertices(:, 1) = new_vertices(:, 1) + x_offset;
                    new_module.beams.(beam_field) = new_vertices;
                end
                
                for pylon = 1: length(fieldnames(tracked_geom.pylons))
                    pylon_field = sprintf('pylon_%d', pylon);
          
                    original_vertices = tracked_geom.pylons.(pylon_field);
                    new_vertices = original_vertices;
                    new_vertices(:, 1) = new_vertices(:, 1) + x_offset;
                    new_module.pylons.(pylon_field) = new_vertices;
                end
    
                % Store the new module in the cell array
                new_modules{1, j} = new_module;
            end
           
        
            % Apply row (y-axis) offsets for panels, beams, and pylons
            for i = 2:num_modules_y
                for j = 1:num_modules_x
                    % Offset for this module
                    y_offset = (i - 1) * row_spacing;

                    modulebaseline = new_modules{1,j};

                    for pan = 1 : length(fieldnames(modulebaseline.panels))
                        panel_field = sprintf('panel_%d', pan);
                    
                        original_vertices = modulebaseline.panels.(panel_field);
                        new_vertices = original_vertices;
                        new_vertices(:, 2) = new_vertices(:, 2) + y_offset;
                        new_modules{i,j}.panels.(panel_field) = new_vertices;
                    end

                    for beam = 1 : length(fieldnames(modulebaseline.beams))
                        beam_field = sprintf('beam_%d', beam);
                    
                        original_vertices = modulebaseline.beams.(beam_field);
                        new_vertices = original_vertices;
                        new_vertices(:, 2) = new_vertices(:, 2) + y_offset;
                        new_modules{i,j}.beams.(beam_field) = new_vertices;
                    end

                    for pylon = 1 : length(fieldnames(modulebaseline.pylons))
                        pylon_field = sprintf('pylon_%d', pylon);
                    
                        original_vertices = modulebaseline.pylons.(pylon_field);
                        new_vertices = original_vertices;
                        new_vertices(:, 2) = new_vertices(:, 2) + y_offset;
                        new_modules{i,j}.pylons.(pylon_field) = new_vertices;
                    end
                 
                end

            end

        end

        
            % Store new modules in the same format as tracked_geom
%             new_vertices_array = struct('panels', cell(num_modules_x, num_modules_y), ...
%                                         'beams', cell(num_modules_x, num_modules_y), ...
%                                         'pylons', cell(num_modules_x, num_modules_y));
%         
%             for i = 1:num_modules_y
%                 for j = 1:num_modules_x
%                     new_vertices_array(i, j).panels = new_modules{i, j}.panels;
%                     new_vertices_array(i, j).beams = new_modules{i, j}.beams;
%                     new_vertices_array(i, j).pylons = new_modules{i, j}.pylons;
%                 end
%             end
%         end

        function plotarray(tracked_geom, panel_faces, beam_faces,pylon_faces, vertices, g )
                     % Function to plot an array of modules and their shadows for a single hour
            % Arguments:
            % - tracked_geom: Structure containing rotated geometry data for each module
            % - panel_faces: Faces for panels
            % - beam_faces: Faces for beams
            % - pylon_faces: Faces for pylons
            % - vertices: Cell array containing shadow vertices for each module
        
            % Initialize figure for plotting

            no_panels = size(fieldnames(tracked_geom{1,1}.panels),1);
            no_beams = size(fieldnames(tracked_geom{1,1}.beams),1);
            no_pylons = size(fieldnames(tracked_geom{1,1}.pylons),1);
            no_vertices = no_pylons+no_beams+no_panels;
            
%             figure();
            figure('WindowState', 'maximized');
            hold on;
            grid on;
        
            % Get the number of rows and columns in the tracked geometry
            [num_rows, num_cols] = size(tracked_geom);
        
            % Loop through each module in the array
            for j = 1:num_cols
                for i = 1:num_rows
                    % Get the rotated geometry for the current module
                    module_geom = tracked_geom{i, j};
        
                    % Check if the module is empty
                    if isempty(module_geom)
                        continue;
                    end
        
                    % Get the number of panels, beams, and pylons for this module
                    no_panels = size(fieldnames(module_geom.panels), 1);
                    no_beams = size(fieldnames(module_geom.beams), 1);
                    no_pylons = size(fieldnames(module_geom.pylons), 1);
        
                    % Loop through each panel
                    for p = 1:no_panels
                        panel_name = ['panel_' num2str(p)];
                        panel_vertices = module_geom.panels.(panel_name);
        
                        % Plot panel vertices using patch
                        patch('Vertices', panel_vertices, 'Faces', panel_faces, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.8);
                    end
        
                    % Loop through each beam
                    for m = 1:no_beams
                        beam_name = ['beam_' num2str(m)];
                        beam_vertices = module_geom.beams.(beam_name);
        
                        % Plot beam vertices using patch
                        patch('Vertices', beam_vertices, 'Faces', beam_faces, 'FaceColor', [0.1, 0.1, 0.1], 'FaceAlpha', 0.5);
                    end
        
                    % Loop through each pylon
                    for c = 1:no_pylons
                        pylon_name = ['pylon_' num2str(c)];
                        pylon_vertices = module_geom.pylons.(pylon_name);
        
                        % Plot pylon vertices using patch
                        patch('Vertices', pylon_vertices, 'Faces', pylon_faces, 'FaceColor', [0.1, 0.1, 0.1], 'FaceAlpha', 0.5);
                    end
        
%                     for k = 1:no_vertices
%                         
%                         % Plot the shadow for this module
%                         shadow_vertices = vertices{i, j}{1,k};
%                         if ~isempty(shadow_vertices) && ~any(isnan(shadow_vertices(:)))
%                             % Check if shadow_vertices is in the correct format
%                             if size(shadow_vertices, 2) == 2
%                                 % Add Z=0 for 2D plot
%                                 shadow_vertices = [shadow_vertices, zeros(size(shadow_vertices, 1), 1)];
%                             end
%                             
%                             % Calculate convex hull for the shadow
%                             [K,~] = convhull(shadow_vertices(:, 1), shadow_vertices(:, 2));
%                             
%                             % Extract X, Y, Z coordinates for plotting
%                             x = shadow_vertices(K, 1);
%                             y = shadow_vertices(K, 2);
%                             z = shadow_vertices(K, 3);
%                             
%                              % Plot the shadow using plot3 with grey color
%                             patch('XData', x, 'YData', y, 'ZData', z, 'FaceColor', [0.7, 0.7, 0.7], ...
%                                  'EdgeColor', [0.7, 0.7, 0.7], 'FaceAlpha', 0.8);
%                             hold on;
%                         end
%                     end

                end
            end
            % Subtract 3600 from x
            x_modified = g - 3600;
            
            % Convert x_modified to a string
            time_str = sprintf('%02d:00 01/06/2023 2x2 Array Shading', x_modified);
            
%             % Replace spaces and special characters with underscores in the title
%             title_str = strrep(time_str, ' ', '_');
%             title_str = strrep(title_str, '/', '_');
%             title_str = strrep(title_str, ':', '_');

            % Customize plot
            title('x Array', 'FontSize', 30);
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            axis equal;
            xlim([-5, 35]);
            ylim([-5, 35]);
            view(-44.9287671232877, 32.7706849315068);
            % Uncomment the line below to set view to a single camera angle
% view(azimuth_angle, elevation_angle);
        % Save the plot with filename as x.png
        % Maximize the figure window
            
            filename = sprintf('%d.png', g);
            saveas(gcf, filename);
        end

     


        function plotmodule_shading(tracked_geometry, panel_faces, beam_faces, pylon_faces, vertices)
        % Function to plot 24 hours worth of geometry and shading data
            % Arguments:
            % - tracked_geometry: Structure containing rotated geometry data for each hour
            % - panel_faces: Faces for panels
            % - beam_faces: Faces for beams
            % - pylon_faces: Faces for pylons
        
            % Initialize figure for plotting
            figure;
            hold on;
            grid on;
            
            % Loop through each hour of the day
            for hour = 3600:3624
                % Get the rotated geometry for the current hour
                hour_name = ['field_' num2str(hour)];
                rotated_hour = tracked_geometry.(hour_name);
                
                no_panels = size(fieldnames(rotated_hour.panels),1);
                no_beams = size(fieldnames(rotated_hour.beams),1);
                no_pylons = size(fieldnames(rotated_hour.pylons),1);
                no_vertices = no_pylons+no_beams+no_panels;
                size(rotated_hour.panels,2) ;

                % Loop through each panel
                for p = 1:no_panels
                    panel_name = ['panel_' num2str(p)];
                    panel_vertices = rotated_hour.panels.(panel_name);
                    
                    % Plot panel vertices using patch
                    patch('Vertices', panel_vertices, 'Faces', panel_faces, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.8);
                end
                
                % Loop through each module
                for m = 1:no_beams
                    beam_name = ['beam_' num2str(m)];
                    beam_vertices = rotated_hour.beams.(beam_name);
                    
                    % Plot module vertices using patch
                    patch('Vertices', beam_vertices, 'Faces', beam_faces, 'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.5);
                end
                
                % Loop through each pylon
                for c = 1:no_pylons
                    pylon_name = ['pylon_' num2str(c)];
                    pylon_vertices = rotated_hour.pylons.(pylon_name);
                    
                    % Plot pylon vertices using patch
                    patch('Vertices', pylon_vertices, 'Faces', pylon_faces, 'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.5);
                end
% %                   Plot the shadow for each panel
%                 for k = 1:no_vertices
%                     shadow_vertices = vertices{hour, k};
%                     x = shadow_vertices(:,1);
%                     y = shadow_vertices(:,2);
%                     [K,A] = convhull(x,y);
%                     if ~isempty(shadow_vertices)
% %                         Check if panel_vertices is in the correct format
%                         if size(shadow_vertices, 2) == 2
%                             shadow_vertices = [shadow_vertices, zeros(size(shadow_vertices, 1), 1)]; % Add Z=0 for 2D plot
%                         end
%                         patch('Vertices', shadow_vertices(K), 'Faces', panel_faces, 'FaceColor', 'r', 'FaceAlpha', 0.5);
%                     end
%                 end
               for k = 1:no_vertices
                    shadow_vertices = vertices{hour, k};
                    
                    if ~isempty(shadow_vertices) && ~any(isnan(shadow_vertices(:)))
                        % Check if shadow_vertices is in the correct format
                        if size(shadow_vertices, 2) == 2
                            % Add Z=0 for 2D plot
                            shadow_vertices = [shadow_vertices, zeros(size(shadow_vertices, 1), 1)];
                        end
                        
                        % Calculate convex hull for the shadow
                        [K,~] = convhull(shadow_vertices(:, 1), shadow_vertices(:, 2));
                        
                        % Extract X, Y, Z coordinates for plotting
                        x = shadow_vertices(K, 1);
                        y = shadow_vertices(K, 2);
                        z = shadow_vertices(K, 3);
                        
                         % Plot the shadow using plot3 with grey color
                        patch('XData', x, 'YData', y, 'ZData', z, 'FaceColor', [0.7, 0.7, 0.7], ...
                             'EdgeColor', [0.7, 0.7, 0.7], 'FaceAlpha', 0.5);
                        hold on;
                    end
                end
                
                % Customize plot
                title(['Geometry and Shading Data - Hour ' num2str(hour)]);
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
                axis equal;
                
              % Save the plot with hour identifier
                saveas(gcf, ['Hour_' num2str(hour) '_plot.png']);
                
                % Clear the plot for the next hour
                clf;
            end
        
        % End of loop
    end

        function [new_vertices_array, max_x_offset, max_y_offset] = createNewModules(num_modules_x, num_modules_y, module_spacing, row_spacing, vertices_array)
            % Inputs:
            % - num_modules_x: Number of modules in a row
            % - num_modules_y: Number of modules in a column
            % - module_spacing: Spacing between modules in the x-direction
            % - row_spacing: Spacing between rows in the y-direction
            % - vertices_array: Cell array containing vertices of the original modules
            
            array_size = num_modules_x * num_modules_y;
    
            % Initialize variables
            new_vertices_array = cell(length(vertices_array), 1);
    
             % Initialize max offsets for x and y
            max_x_offset = -inf;
            max_y_offset = -inf;
    
            % Loop through each row of the original vertices array
            for time_step = 1:length(vertices_array)
                % Create cell array to store new modules for this row
                new_row_modules = cell(num_modules_x, num_modules_y);
    
                new_row_modules{1,1} = vertices_array(time_step,:);
        
                for col = 2:num_modules_x
                    
                    % Calculate the offset for this module along the x-direction
                    x_offset = (col-1) * module_spacing;
    
                    for module_part = 1:length(new_row_modules{1,1})
                        module_part_baseline = new_row_modules{1, 1}{1, module_part};
    
                        if any(any(isnan(module_part_baseline)))
                            new_row_modules{1, col}{1, module_part} = new_row_modules{1, 1}{1, module_part};
                        else
                            new_row_modules{1, col}{1, module_part} = new_row_modules{1, 1}{1, module_part};
                            new_row_modules{1, col}{1, module_part}(:, 1) = new_row_modules{1, 1}{1, module_part}(:, 1) + x_offset;
    
                             % Update max offset for x
                            max_x_offset = max(max_x_offset, x_offset);
    
                        end
    
                        new_row_modules{1, col}{1, module_part}(:, 2:end) = new_row_modules{1, 1}{1, module_part}(:, 2:end);
    
                    end
    
                end
    
                for col = 1:num_modules_x
    
                    for row = 2:num_modules_y
    
                        y_offset = (row-1) * row_spacing;
    
                        for module_part = 1:length(new_row_modules{1,1})
                            module_part_baseline = new_row_modules{1, col}{1, module_part};
        
                            if any(any(isnan(module_part_baseline)))
                                new_row_modules{row, col}{1, module_part} = new_row_modules{1, col}{1, module_part};
                            else
                                new_row_modules{row, col}{1, module_part} = new_row_modules{1, col}{1, module_part};
                                new_row_modules{row, col}{1, module_part}(:, 2) = new_row_modules{1, 1}{1, module_part}(:, 2) + y_offset;
                                 
                               % Update max offset for y
                                max_y_offset = max(max_y_offset, y_offset);
                            end
                        end
                    end
                end
    
                new_vertices_array{time_step,1} = new_row_modules;
            end
    %       % Return the overall max offsets for x and y
            max_x_offset = max_x_offset + module_spacing;
            max_y_offset;  
        end

        function turn_angle = tracking_angle(azapparent, elapparent, latitude, axis_angle)

                    % Inputs:
        %   SunZen - a scalar or vector of apparent (refraction-corrected) zenith
        %     angles in decimal degrees. If SunZen is a vector it must be of the
        %     same size as all other vector inputs. SunZen must be >=0 and <=180.
        %   SunAz - a scalar or vector of sun azimuth angles in decimal degrees.
        %     If SunAz is a vector it must be of the same size as all other vector
        %     inputs. SunAz must be >=0 and <=360. The Azimuth convention is defined
        %     as degrees East of North (e.g. North = 0, East = 90, West = 270).
        %   Latitude - a scalar value denoting which hempisphere the tracker is
        %     in. The exact latitude is NOT required, any positive number denotes
        %     the northern hemisphere, any negative number denotes the southern
        %     hemisphere, a value of 0 is assumed to be northern hemisphere.
        %   AxisTilt - a scalar value denoting the tilt of the axis of rotation
        %     (i.e, the y-axis defined by AxisAzimuth) with respect to horizontal, 
        %     in decimal degrees. AxisTilt must be >=0 and <=180.
        %   AxisAzimuth - a scalar value denoting the compass direction along which
        %     the axis of rotation lies, in decimal degrees. Again, the convention 
        %     is defined as degrees East of North (e.g. North = 0, East = 90, 
        %     West = 270). AxisAzimuth must be >=0 and <=360.
        %   MaxAngle - a scalar value denoting the maximum rotation angle, in
        %     decimal degrees, of the one-axis tracker from its horizontal position
        %     (horizontal if AxisTilt = 0). MaxAngle must be <=180 and >=0. A
        %     MaxAngle of 90 degrees allows the tracker to rotate to a vertical
        %     position to point the panel towards a horizon.  
        %     MaxAngle of 180 degrees allows for full rotation.
       
        axis_tilt = 0;
        max_angle = 50;
      
        zenith = 90 - elapparent;

        [TrkrTheta, AOI, SurfTilt, SurfAz]= pvl_singleaxis(zenith, azapparent, latitude, axis_tilt, axis_angle, max_angle);

        turn_angle= TrkrTheta;


        end

        function [is_shaded, vertices, max_x_offset,max_y_offset] = shadingFactor(year,resolution,UTC,latitude,...
                longitude,altitude,blocking_object,setup_ground, num_modules_y, num_modules_x, module_spacing, row_spacing, panel_faces, beam_faces, pylon_faces, x)
            % returns the shading factors for a given time and location
            all_TT = struct();  % create structure to save tables
            t1 = datetime(year, 1, 1, 0, 0, 0);                             % start time
            t2 = datetime(year, 12, 31, 24, 0, 0);                          % end time
            t = t1:minutes(resolution):t2;                                  % time resolution 
            time_num = datenum(t);
            
            time = pvl_maketimestruct(time_num,UTC);  % time structure
            location = pvl_makelocationstruct(latitude,longitude,altitude); % Create a structure to define a site location
            
            square = [0, 0, 0;          0, 4, 0;         -2, 4, 0;         -2, 0, 0];
%             names = fieldnames(obj.coefficients)
            % sun position angles
            [azimut, elevation, elAparent]=pvl_spa(time, location);%Calculates the position of the sun given time, location, and optionally pressure and temperature
            
            sun = zeros(numel(azimut), 3);
            
            for i = 1:numel(azimut)
                [sunX, sunY, sunZ] = sph2cart(-azimut(i)*pi/180,elAparent(i)*pi/180,1);%SPH2CART - Transform spherical to Cartesian coordinates.
                

                sun(i,1:3) = [sunX, sunY, sunZ];
                x1= createTable.tracking_angle(azimut(i),elAparent(i),latitude, 0);
                zenith = 90 - elAparent(i);
                sun(i,4) = x1;
                sun(i,5) = azimut(i);
                sun(i,6) = zenith;
                sun(i,7) = 90;


                
            end
                        

            tracked_geometry = struct();
            for i = 1:length(sun(:,4))
                tracking_angle = round(sun(i, 4));
                
                % Construct angle name based on convention
                if tracking_angle < 0
                    anglestring = ['minus_' num2str(abs(tracking_angle))];
                else
                    anglestring = num2str(tracking_angle);
                end
                angle_name = ['angle_' anglestring];
                
                % Use i as a string for the field name
                field_name = ['field_' num2str(i)];
                
                tracked_geometry.(field_name) = blocking_object.(angle_name);
            end

            field_name = ['field_', num2str(x)];

            vertices = createShadow.createShadow2(square, tracked_geometry, sun);
            
%             createTable.plotmodule_shading(tracked_geometry, panel_faces, beam_faces, pylon_faces, vertices)
            
            [new_vertices_array, max_x_offset, max_y_offset] = createTable.createNewModules(num_modules_x, num_modules_y, module_spacing, row_spacing, vertices);

            new_modules = createTable.createfullModules(num_modules_x, num_modules_y, module_spacing, row_spacing, tracked_geometry.(field_name));%, new_vertices_array{3613,1}

            createTable.plotarray(new_modules,panel_faces, beam_faces, pylon_faces,  new_vertices_array{x,1}, x);
            
%             is_shaded = 1;
            is_shaded = createShadow.check_shaded(new_vertices_array,setup_ground);



        end

        function [diffuse_factor, gcr] = diffuseappr(max_x_offset, max_y_offset, panel_area,total_modules, max_height, num_panels)
            
            top_area = (max_x_offset*max_y_offset);
            side_area_x =  2*max_x_offset*max_height;
            side_area_y = 2*max_y_offset*max_height;
            total_area = top_area + side_area_y + side_area_x;
            

            array_area = panel_area * num_panels * total_modules;
            
            diffuse_factor = array_area/total_area;
            gcr = (array_area/ top_area);


    end
    end

end


