clear
clc




% ---------------------------- SYSTEM INPUTS --------------------

% Call the getTMYData method
latitude = 40;
longitude = -5;
startYear = 2000;
endYear = 2020;

% Define the number of panels
num_panels = 4;  % Adjust as needed


% Define the array parameters
num_modules_x = 3;
num_modules_y = 4;
module_spacing = 10.5;
row_spacing = 15;
total_modules = num_modules_y * num_modules_x;


%  field paramteters

field_length = 35;
field_width = 30;
field_area = (field_length * field_width / 10000);


% % ------------------------- Setup Agrivoltaic geometry-----------------
% 
% % -------------------theory ----------------------------------
% % create each part individually in inventor with [0,0,0] as reference point
% % offset the vertices of each part by the required distance in [x,y,z]
% % results in a fully formed assembly in MATLAB
% %-----------------------------------------------------------
% 


% ---------------------------- get solar data



% Create an instance of the PVGISData class
pvgis = PVGISData();

% pv_data = pvgis.getPVData(latitude, longitude, 25.2);
% disp(pv_data);
% pv_data = pvgis.getPVData(latitude, longitude, 0.715);
% % disp(pv_data);

tmy_data = pvgis.getTMYData(latitude, longitude, startYear, endYear);
% G_h_ - Global irradiance on the horizontal plane - 'W/m2'
% Gb_n_ - Beam/direct irradiance on a plane always normal to sun rays - 'W/m2'
% Gd_h_  - Diffuse irradiance on the horizontal plane - 'W/m2'


% Display or process the TMY data as needed
disp(tmy_data);

% % Call the getPVData method
% peakPower = 4000; % kW
% 
% pv_data = pvgis.getPVData(latitude, longitude, peakPower);

% E_d - Average daily energy production from the given system - 'kWh/d'
% E_m - Average monthly energy production from the given system - 'kWh/mo'
% Display or process the PV data as needed
% disp(pv_data);

temp_data = pvgis.gettempData(latitude,longitude,  1);

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



% % Plot the original panel and the modified panels
% figure;
% hold on;
% for i = 1:(num_panels)
%     i;
%     patch('Vertices', all_panels{i}.vertices, ...
%           'Faces', all_panels{i}.faces, ...
%           'FaceColor', [1.0, 0.7, 0.0], ...
%           'EdgeColor', [0.1, 0.1, 0.1], ...
%           'FaceAlpha', 0.7);
% end
% % Plot the beams
% for i = 1:length(all_beams)
%     patch('Vertices', all_beams{i}.vertices, ...
%           'Faces', all_beams{i}.faces, ...
%           'FaceColor', [0.7, 0.0, 1.0], ...
%           'EdgeColor', [0.1, 0.1, 0.1], ...
%           'FaceAlpha', 0.7);
% end
% 
% % Plot the pylons
% for i = 1:length(all_pylons)
%     patch('Vertices', all_pylons{i}.vertices, ...
%           'Faces', all_pylons{i}.faces, ...
%           'FaceColor', [0.0, 0.7, 1.0], ...
%           'EdgeColor', [0.1, 0.1, 0.1], ...
%           'FaceAlpha', 0.7);
% end
% 
% % Add labels and title
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('Original Panel, Modified Panels, and Modified Beams');
% axis equal;
% grid on;
% legend('Original Panel', 'Modified Panels', 'Modified Beam 1', 'Modified Beam 2');
% 
% hold off;





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
angles_degrees = -50:1:50;
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




% ---------------------------- plot rotated structure ---------------

% % Define colors for the panels
% panel_colors = {'r', 'g', 'b', 'c', 'm', 'y', 'k'};
% 
% % Loop through each angle in the rotated_struct
% angles = fieldnames(rotated_struct);
% for i = 1:numel(angles)
%     current_angle = angles{i};
%     current_angle_struct = rotated_struct.(current_angle);
%     
%     % Create a new figure for each angle
%     figure;
%     hold on;
%     title(['Panels and Beams at Rotation Angle ' current_angle]);
%     
%     % Loop through each panel in the current angle struct
%     if isfield(current_angle_struct, 'panels')
%         panel_names = fieldnames(current_angle_struct.panels);
%         for j = 1:numel(panel_names)
%             current_panel = current_angle_struct.panels.(panel_names{j});
%             
%             % Plot the panel using vertices and faces
%             patch('Vertices', current_panel, ...
%                   'Faces', all_panels{j}.faces, ...
%                   'FaceColor', panel_colors{j}, ...
%                   'EdgeColor', 'k', ...
%                   'FaceAlpha', 0.7);
%         end
%     end
%     
%     % Loop through each beam in the current angle struct
%     if isfield(current_angle_struct, 'beams')
%         beam_names = fieldnames(current_angle_struct.beams);
%         for k = 1:numel(beam_names)
%             current_beam = current_angle_struct.beams.(beam_names{k});
%             
%             % Plot the beam using vertices and faces
%             patch('Vertices', current_beam, ...
%                   'Faces', all_beams{k}.faces, ...
%                   'FaceColor', 'y', ...
%                   'EdgeColor', 'k', ...
%                   'FaceAlpha', 0.7);
%         end
%     end
%     
%     % Set axis labels and properties
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');
%     axis equal;
%     grid on;
%     
%     hold off;
% end
% 




% ------------- Create ground points to be analysed -------------------
% set an area of ground to be analysed, mesh a numver of points within this
% area depending on level of fidelity


setup_ground = createGround(10,20,5,20,10,10);
points_x = setup_ground.points(:,1);
points_y = setup_ground.points(:,2);
% setup_ground.draw
     

% 
% sol_vector = [-0.31357059  0.71199223  0.62828381];
% find_shade = createShadow(setup_ground.squares(1:4,:),original_panel.Vertices,sol_vector);
% find_shade.coef(setup_ground.squares(1:4,:))
% % find_shade.vertices;
%find_shade.draw


% ------------------- Calculate shading of ground --------------------



% for x = 3606:3622

% % year,resolution,UTC,latitude,longitude,altitude,blocking_object,setup_ground
    [shading_factor, vertices, max_x_offset, max_y_offset] = createTable.shadingFactor(2019,60,1,37.193,-5.853,36,rotated_struct, ...
        setup_ground.points, num_modules_y, num_modules_x, module_spacing, row_spacing, all_panels{1}.faces, all_beams{1}.faces, all_pylons{1}.faces, num_modules_y);

% end

[diffuse_factor,gcr] = createTable.diffuseappr(max_x_offset, max_y_offset, panel_area,total_modules, max_height, num_panels);
% % 
% array.drawArray;
% 
% % ------------------ Turn Shading into Irradiance ----------------------
% diffuse_factor = 0.15;

myData = Irradiance(tmy_data);
% % 

% % Specify the filename
% filename = '2x2shading_factor_data.mat';
% 
% % Save the variable to a .mat file
% save(filename, 'shading_factor');

% shading_factor = importdata('mod3_4row7space_5pan.mat');
% shading_factor = shading_factor.shade_array;

myData = myData.calcdata(shading_factor,diffuse_factor);
% 
myData = myData.dailyenergy();

% myData.shadingfactorplot();




% % ------------------ Calculate Crop yield ------------------------------

% Instantiate TomatoModel with parameters for SunnySD cultivar
tomato_sunny_SD = CropModel(2800, 0.68, 520, 400, 6, 26, 1.00, 100, 5, 32, 45, 0.07, 2.5, 'tmyiradata.csv', 1.0923);
tomato_sunny_SD = tomato_sunny_SD.calculate_temp_diff(0);
tomato_sunny_SD = tomato_sunny_SD.fSolar();
tomato_sunny_SD = tomato_sunny_SD.fTemp();
tomato_sunny_SD = tomato_sunny_SD.fHeat() ;


tomato_sunny_SD = tomato_sunny_SD.calculate_i50_list();
tomato_sunny_SD = tomato_sunny_SD.calculate_fco2_list();
tomato_sunny_SD = tomato_sunny_SD.runModel(myData.energyvalues,45,210);

% contourplo(setup_ground.points(:,1),setup_ground.points(:,2),myData.energyvalues(:,166));
% ---------------------- OUTPUTS ETC. --------------------------------

peakmodulePower = (num_panels*0.7);
% Call the getPVData method
% peakPower = 4000; % kW

pv_data = pvgis.getPVData(latitude, longitude, peakmodulePower);
disp(pv_data);

tomato_sunny_SD = tomato_sunny_SD.plotcrops(all_panels, all_beams,all_pylons,points_x,points_y);


%  outputs

yearly_module_output = pv_data.outputs.totals.vertical_axis.E_y;

final_site_kWh_output = yearly_module_output * floor((field_length/module_spacing)) * floor((field_width/row_spacing));

final_agricultural_output = tomato_sunny_SD.mean_outputs*field_area;




% % Instantiate TomatoModel with parameters for Agriset761 cultivar
% tomato_agriset761 = CropModel(2300, 0.50, 550, 300, 6, 26, 1.00, 100, 5, 32, 45, 0.07, 2.5);
% 
% 
% wheat_model = CropModel(2200, 0.36, 480, 200, 0, 15, 1.24, 100, 25, 34, 45, 0.08, 0.4, '2005iradata.csv',1.0923);
% wheat_model = wheat_model.calculate_temp_diff(0);
% wheat_model = wheat_model.fSolar()
% wheat_model = wheat_model.fTemp()
% wheat_model = wheat_model.fHeat() 
% 
% 
% wheat_model = wheat_model.calculate_i50_list()
% wheat_model = wheat_model.calculate_fco2_list()
% wheat_model = wheat_model.runModel(myData.energyvalues,45,210)
% displacement
% array = createPVarray(setup_panel,1,1,5,5,0, displacement);
% 
% array.draw(setup_ground, vertices,points_x,points_y,wheat_model.sum_outputs)