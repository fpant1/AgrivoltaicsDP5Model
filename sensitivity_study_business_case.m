clear
clc

% Fixed system inputs
latitude = 37.193;
longitude = -5.853;
startYear = 2000;
endYear = 2020;
UTC = 1;
altitude = 36;
resolution = 60; % Resolution in hours
growthtosen = 300;

planting = 77;
harvesting = 217;



% Define the array parameters
num_modules_x = 4;
num_modules_y = 5;
module_spacing = 10.5;
total_modules = num_modules_y * num_modules_x;

field_length = 400;
field_width = 250;
field_area = (field_length * field_width / 10000);


% Load the weather and data files
WEATHERFILE = importdata('WEATHERFILE.mat');
DATA = importdata('DATA.mat');

% Initialize table to store the results
results = table();

% Sensitivity study ranges
num_panels_range = [3, 5, 6, 7];
row_spacing_range = [5, 7.5,10,12.5, 15];
crop_protection_range = [0, 1];





% Tracking algorithm inputs -----------

% parameters dont change

Tracking_limit = 60; % Tracking rotation angle limit (clockwise)
Position_A = 30; % Wind stow position


%  for electrical yield put electrical put eto high and position b to 60
%  critical times set them to 0, 0
Electrical_focussed_tracking = [];

% for standard tracking put critical values to zero



Position_B = 0; % heavy precipitation position
Position_C = 60; % low precipitation position
heavy_rain = 0.004; 
CT_s = 3648; % Start hour for critical growth period
CT_e = 4368; % End hour for critical growth period
ETo_average = 0.06976680980;

%%%%%%%%%%% SUN ANGLES %%%%%%%%%%%
year = 2019; % Integer representing the year
resolution = 60; % Resolution in hours
UTC = 1;

altitude = 36;

%%%%%%%%%%% Load the weather and data files %%%%%%%%%%%
WEATHERFILE = importdata('WEATHERFILE.mat');
DATA = importdata('DATA.mat');

% Assuming WEATHERFILE is already imported and contains the hourly rainfall data in its second column.
rainfall_hourly_table = WEATHERFILE(:,2); % Extract hourly rainfall table
rainfall_hourly = table2array(rainfall_hourly_table); % Convert table to array

% Multiply each value by 100
rainfall_hourly = rainfall_hourly * 1000;
% Initialize an empty array to store daily rainfall totals.
rainfall_daily = [];

% Calculate the number of days in the dataset (assuming complete days).
numberOfDays = floor(height(rainfall_hourly_table) / 24);

% Loop over the hourly data to sum every 24 hours into one day.
for day = 1:numberOfDays
    dailyTotal = sum(rainfall_hourly((day-1)*24+1 : day*24));
    rainfall_daily = [rainfall_daily; dailyTotal]; % Append the daily total to the array.
end

% rainfall_daily now contains the daily rainfall totals.

%%%%%%%%%%% Call the function %%%%%%%%%%%
[total_power_yearly, tracking_angles] = Tracking_Algorithm(Tracking_limit, Position_A, Position_B, Position_C, heavy_rain, CT_s, CT_e, ETo_average, year, resolution, UTC, latitude, longitude, altitude, WEATHERFILE, DATA);

%%%%%%%%%%% Display the result %%%%%%%%%%%
disp(['Total yearly power output: ', num2str(total_power_yearly), ' kWh']);





%  ---------------- API call data  ---------------------------

% Create an instance of the PVGISData class
pvgis = PVGISData();
tmy_data = pvgis.getTMYData(latitude, longitude, startYear, endYear);
myData = Irradiance(tmy_data);

% %  API call for temperature data
temp_data = pvgis.gettempData(latitude,longitude,  1);


% ------------- Create ground points to be analysed -------------------
% set an area of ground to be analysed, mesh a numver of points within this
% area depending on level of fidelity


setup_ground = createGround(-20,-10,-20,-10,4,4);
points_x = setup_ground.points(:,1);
points_y = setup_ground.points(:,2);
% setup_ground.draw

% Loop through each combination of parameters
counter = 1; % Counter for indexing into the results table
num_panels = 4;
row_spacing = 6;
% for num_panels = num_panels_range
%     for row_spacing = row_spacing_range
%         for crop_protection = crop_protection_range
            

        % --------------------------------Setup Structure -----------
        
        % ---------------------BEAMS
        
        
        
        all_beams{1}.vertices = [1, 2, 3];
        all_beams{1}.Faces = [1, 2, 3];
        
        % Create thermal models
        model = createpde("structural");
        model1 = createpde("structural");
        model2 = createpde("structural");
        model3 = createpde("structural");
        
        
        
        
        
        % Import a new geometry (boxsectiontest.stl) and include it in the second model
        beam = importGeometry(model1, "boxsectiontest.stl");
        
        % Extract the vertices and faces of the imported beam geometry
        [F_beam, V_beam] = beam.allDisplayFaces();
        
        % Scale the vertices from mm to meters
        V_beam = V_beam / 1000; % Convert from mm to meters
        
        min_x = inf;
        max_x = -inf;
        % Loop through each vertex to find min and max x values
        for i = 1:size(V_beam, 1)
            x = V_beam(i, 1); % x coordinate of the vertex
            
            % Update min and max x values
            if x < min_x
                min_x = x;
            end
            
            if x > max_x
                max_x = x;
            end
        end
        
        % Calculate the difference between max and min x values
        x_difference = max_x - min_x;
        
        % Beam offset values in meters
        beam_offset1 = [0.3803, 0.45975, 4.262]; % Converted from mm to meters
        beam_offset2 = [0.3803, 1.85925, 4.262]; % Converted from mm to meters
        
        
        
        % Store the modified beam vertices and faces
        all_beams{1}.vertices = V_beam + beam_offset1;
        all_beams{1}.faces = F_beam;
        
        all_beams{2}.vertices = V_beam + beam_offset2;
        all_beams{2}.faces = F_beam;
        
        beam2 = importGeometry(model2,"boxsectiontestcross.stl");
        
        % Extract the vertices and faces of the imported beam geometry
        [F_beam, V_beam] = beam2.allDisplayFaces();
        
        % Scale the vertices from mm to meters
        V_beam = V_beam / 1000; % Convert from mm to meters
        
        
        % Beam offset values in meters
        crossbeam_offset1 = [0.3778, 0.52725, 4.262]; % Converted from mm to meters
        crossbeam_offset2 = [0.3778 + x_difference+ 0.0675, 0.52725, 4.262];
        
        % Store the modified beam vertices and faces
        all_beams{3}.vertices = V_beam + crossbeam_offset1;
        all_beams{3}.faces = F_beam;
        
        all_beams{4}.vertices = V_beam + crossbeam_offset2;
        all_beams{4}.faces = F_beam;
        
        
        %  ----------------PYLONS 
        
        % Import a new geometry (pylons.stl) and include it in the fourth model (model3)
        pylons = importGeometry(model3, "pylon.stl");
        
        % Extract the vertices and faces of the imported pylons geometry
        [F_pylons, V_pylons] = pylons.allDisplayFaces();
        
        % Scale the vertices from mm to meters
        V_pylons = V_pylons / 1000; % Convert from mm to meters
        
        % Create space in all_pylons for storing pylon data
        all_pylons = cell(1, 2);
        
        % Pylon offset values in meters (to be edited)
        pylon_offset1 = [0, 1.102,4]; % Placeholder
        pylon_offset2 = [x_difference + 0.5291,1.102, 4]; % Placeholder
        
        % Store the modified pylon vertices and faces (placeholders for now)
        all_pylons{1}.vertices = V_pylons + pylon_offset1;
        all_pylons{1}.faces = F_pylons;
        
        max_height = max(all_pylons{1}.vertices(:,3));
        
        all_pylons{2}.vertices = V_pylons + pylon_offset2;
        all_pylons{2}.faces = F_pylons;
        
        
        % --------------------PANELS 
        
        % Import the original geometry and include it in the model
        original_panel = importGeometry(model, "testpanel.stl");
        
        % Extract the vertices and faces of the imported original geometry
        [F, V] = original_panel.allDisplayFaces();
        
        % Scale the vertices from mm to meters
        V = V / 1000; % Convert from mm to meters
        
        % Panel offset values in meters
        panel_offset = [1.6833, 0, 4.262]; % Converted from mm to meters
        
        % Create a new structure to store the panel geometry
        % panel_geometry = struct('vertices', V, 'faces', F);
        
        % Offset the panel geometry into position
        all_panels{1}.vertices(:, 1) = V(:, 1) + panel_offset(1); % Offset in x
        all_panels{1}.vertices(:, 2) = V(:, 2) + panel_offset(2); % Offset in y
        all_panels{1}.vertices(:, 3) = V(:, 3) + panel_offset(3); % Offset in z
        
        % Calculate the area for the first panel
        % Area in x-direction
        panel_width_x = max(all_panels{1}.vertices(:, 1)) - min(all_panels{1}.vertices(:, 1));
        % Area in y-direction (assuming same y-values for all vertices)
        panel_width_y = max(all_panels{1}.vertices(:, 2)) - min(all_panels{1}.vertices(:, 2));
        
        % Calculate the total area of the first panel
        panel_area = panel_width_x * panel_width_y;
        
        
        
        all_panels{1}.faces = F;
        
        
        
        
        % Loop through each panel to create and append
        for k = 1:num_panels
            % Calculate the offset for this panel
            if k == 1
                % First panel has no offset
                panel_offset_x = 0;
        %         all_panels{k}.vertices(:, 1) = all_panels{1}.vertices(:, 1) + panel_offset_x;
        %         all_panels{k}.vertices(:, 2) = all_panels{1}.vertices(:, 2) + panel_offset(2); % Same y-offset
        %         all_panels{k}.vertices(:, 3) = all_panels{1}.vertices(:, 3)
            else
                % Calculate even distribution of offsets up to (x_difference - 1.303)
                panel_offset_x = (k - 1) * ((x_difference - 1.303) / (num_panels - 1));
                % Offset the panel geometry in the x-direction
                all_panels{k}.vertices(:, 1) = all_panels{1}.vertices(:, 1) + panel_offset_x;
                all_panels{k}.vertices(:, 2) = all_panels{1}.vertices(:, 2);  % Same y-offset
                all_panels{k}.vertices(:, 3) = all_panels{1}.vertices(:, 3);  % Same z-offset
            end
            
            
            all_panels{k}.faces = F; % Same faces for all panels
            
        
        end
        
        
        
        
        % -----------------------------  tracking ------------------------
        
        % Extract vertices of beam 1
        beam1_vertices = all_beams{3}.vertices;
        
        % Extract vertices of beam 2
        beam2_vertices = all_beams{4}.vertices;
        
        % Find the minimum and maximum x and z coordinates of both beams
        min_y = min([min(beam1_vertices(:, 2)), min(beam2_vertices(:, 2))]);
        max_y = max([max(beam1_vertices(:, 2)), max(beam2_vertices(:, 2))]);
        min_z = min([min(beam1_vertices(:, 3)), min(beam2_vertices(:, 3))]);
        max_z = max([max(beam1_vertices(:, 3)), max(beam2_vertices(:, 3))]);
        
        % Calculate the midpoint
        midpoint_x = 1 ;
        midpoint_y = (min_y + max_y) / 2; % Y-coordinate is 1
        midpoint_z = (min_z + max_z) / 2;
        
        % Create the midpoint as a 1x3 vector
        midpoint = [midpoint_x, midpoint_y, midpoint_z];
        
        % Direction vector along the y-axis
        v = [1, 0, 0];
        
        % Choose a value to shift along the y-axis
        shift_distance = 10;  % Adjust this value as needed
        
        
        
        
        
        % Choose rotation angles from -50 to 50 degrees (in radians)
        angles_degrees = -60:1:60;
        angles_radians = deg2rad(angles_degrees);
        
        % Initialize structure to store rotated vertices for each angle
        rotated_struct = struct();
        
        % Loop through each angle
        for i = 1:length(angles_radians)
            % Calculate rotated vertices for the current angle for all panels
            angle = angles_degrees(i);
            if angle < 0
                anglestring = ['minus_' num2str(abs(angles_degrees(i)))];
                angle_name = ['angle_' anglestring];
            else
                angle_name = ['angle_' num2str(angles_degrees(i))];
            end
            
            % Initialize structure for current angle
            angle_rotated = struct();
            
            % Initialize structure for panels and beams at current angle
            panels = struct();
            beams = struct();
            pylons = struct();
            
            % Loop through each panel
            for p = 1:length(all_panels)
                % Panel vertices
                panel_vertices = all_panels{p}.vertices.';
                
                % Calculate rotated vertices for the current panel
                [rotated_panel, R, t] = AxelRot(panel_vertices, angle, v, midpoint);
                
                % Store rotated vertices for current panel in panels structure
                panel_name = ['panel_' num2str(p)];
                panels.(panel_name) = rotated_panel.';
            end
            
            % Loop through each beam
            for b = 1:length(all_beams)
                % Beam vertices
                beam_vertices = all_beams{b}.vertices.';
                
                % Calculate rotated vertices for the current beam
                [rotated_beam, R, t] = AxelRot(beam_vertices, angle, v, midpoint);
                
                % Store rotated vertices for current beam in beams structure
                beam_name = ['beam_' num2str(b)];
                beams.(beam_name) = rotated_beam.';
            end
            
            for c = 1:length(all_pylons)
                pylon_vertices = all_pylons{c}.vertices;
                pylon_name = ['pylon_' num2str(c)];
                pylons.(pylon_name) = pylon_vertices;
            end
        
            % Store panels and beams structure within current angle structure
            angle_rotated.panels = panels;
            angle_rotated.beams = beams;
            angle_rotated.pylons = pylons;
            
            % Store rotated panel and beam structure in main structure for current angle
            rotated_struct.(angle_name) = angle_rotated;
        end


        

            % ------------------- Calculate shading of ground --------------------

            % % year,resolution,UTC,latitude,longitude,altitude,blocking_object,setup_ground
            [shading_factor, vertices, max_x_offset, max_y_offset] = createTable.shadingFactor(2019,60,1,37.193,-5.853,36,rotated_struct, ...
                    setup_ground.points, num_modules_y, num_modules_x, module_spacing, row_spacing, all_panels{1}.faces, all_beams{1}.faces, all_pylons{1}.faces, num_modules_y, tracking_angles);
       
            [diffuse_factor,gcr] = createTable.diffuseappr(max_x_offset, max_y_offset, panel_area,total_modules, max_height, num_panels);
            
            myData = myData.calcdata(shading_factor,diffuse_factor);
            myData = myData.dailyenergy();

            tomato_sunny_SD = CropModel(2800, 0.68, 520, 400, 6, 26, 1.00, 100, 5, 32, 45, 0.07, 2.5, 'case_studies_temp_modelling.csv', 1.0923, rainfall_daily, 0, gcr);
            tomato_sunny_SD = tomato_sunny_SD.calculate_temp_diff(0);
            tomato_sunny_SD = tomato_sunny_SD.fTemp();
            tomato_sunny_SD = tomato_sunny_SD.fHeat();
            tomato_sunny_SD = tomato_sunny_SD.calculate_fco2_list();
            tomato_sunny_SD = tomato_sunny_SD.fWater();
            tomato_sunny_SD = tomato_sunny_SD.calculate_i50_list();
            tomato_sunny_SD = tomato_sunny_SD.fSolar(growthtosen);
            tomato_sunny_SD = tomato_sunny_SD.runModel(myData.energyvalues,planting, harvesting);
            

            total_module_power_yearly = total_power_yearly * num_panels;
            num_modules = floor((field_length/module_spacing));
            num_rows = floor((field_width/row_spacing));
            total_panels = num_panels * num_modules * num_rows;
            final_site_kWh_output = total_module_power_yearly * num_modules * num_rows;
            installed_capacity_kW = total_panels * 0.715;
            
            final_agricultural_output = tomato_sunny_SD.mean_outputs*field_area;

            results(counter,:) = table(total_panels,num_modules, num_rows, num_panels, row_spacing, gcr, installed_capacity_kW, final_site_kWh_output, final_agricultural_output);
                        counter = counter + 1;
%         end
%     end
% end

% Save the results table to a file
save('no_protect_plotensitivity_study_results.mat', 'results');

% Display the first few rows of the results for verification
disp(results(1:5,:));