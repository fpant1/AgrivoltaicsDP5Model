
function [total_power_yearly, panel_angle] = Tracking_Algorithm(Tracking_limit, Position_A, Position_B, Position_C, heavy_rain, CT_s, CT_e, ...
    ETo_average, year, resolution, UTC, latitude, longitude, altitude, weather_file_path, data_file_path)

        t1 = datetime(year, 1, 1, 0, 0, 0);                             % start time
        t2 = datetime(year, 12, 31, 23, 0, 0);                          % end time
        t = t1:minutes(resolution):t2;                                  % time resolution 
        time_num = datenum(t);
                    
        time = pvl_maketimestruct(time_num,UTC);  % time structure
        location = pvl_makelocationstruct(latitude,longitude,altitude); % Create a structure to define a site location
                    
        [azimut, elevation, elAparent]=pvl_spa(time,location);
        
        sun = zeros(numel(azimut),3);
        
        azaparent = 0;
        axis_angle = 0;
        axis_tilt = 0;
        max_angle = 90;
        
        for i = 1:numel(azimut)
                        [sunX, sunY, sunZ] = sph2cart(-azimut(i)*pi/180,elAparent(i)*pi/180,1);%SPH2CART - Transform spherical to Cartesian coordinates.
                        zenith = 90 - elAparent(i);
                        name = azimut(i);
                        sun(i,1:3) = [sunX, sunY, sunZ];
                        [TrkrTheta, AOI, SurfTilt, SurfAz]= pvl_singleaxis(zenith, name, latitude, axis_tilt, axis_angle,max_angle);
                        sun(i,4) = TrkrTheta;  
                        sun(i,5) = AOI;
        end
        
        %%%%%%%%%%% EXTRACTING WEATHER FILE DATA %%%%%%%%%%%
        WEATHERFILE = importdata('WEATHERFILE.mat');
        weather = table2array(WEATHERFILE);
        
        Sun_angle = sun(:,4);
        Wind_speed = weather(:,1);
        Water_level = weather(:,2);
        ETo = weather(:,3);
        Wind_dir = weather(:,4);
        
        %%%%%%%%%%% PANEL ANGLE ALGORITHM %%%%%%%%%%%
        
        % Combine wind speed and direction into a single variable (+- speed depending on
        % direction)
        Wind_speed_dir = Wind_speed;
        index = Wind_dir > 180;
        Wind_speed_dir(index) = -Wind_speed_dir(index);
        
        % Integrate ETo with water_level to limit tracking due to rain 
        average_ETo = mean(ETo);
        index_below_average = ETo < average_ETo;
        index_low_precipitation = Water_level < heavy_rain;
        combined_index = index_below_average & index_low_precipitation;
        modified_water_level = Water_level;
        modified_water_level(combined_index) = 0;
        
        % Initialize panel_angle with Sun_angles
        panel_angle = Sun_angle;
        
        % Define night_time variable
        night_time = (Sun_angle == 0) & (([0; Sun_angle(1:end-1)] == 0) | ([Sun_angle(2:end); 0] == 0));
        
        % Adjust panel_angle during critical crop growth period if Sun_angle is not 0 before moving
        for i = CT_s:CT_e
            if ~night_time(i)
                if Sun_angle(i) > 0
                    panel_angle(i) = Sun_angle(i) - 90;
                elseif Sun_angle(i) < 0
                    panel_angle(i) = Sun_angle(i) + 90;
                end
            end
        end
        
        % Adjust panel_angle based on Water_level
        indices = modified_water_level > 0 & modified_water_level <= 0.004;
        panel_angle(indices) = Position_C;
        panel_angle(modified_water_level > heavy_rain) = Position_B;
        
        % Panel hail protection (60 deg)
        indices = Water_level > heavy_rain;
        panel_angle(indices) = Position_C;
        
        % Adjust panel_angle based on Wind_speed_dir
        panel_angle(Wind_speed_dir > 16) = Position_A;
        panel_angle(Wind_speed_dir < -16) = -Position_A;
        
        % Limit panel_angle between tracking angle limits
        panel_angle(panel_angle > Tracking_limit) = Tracking_limit;
        panel_angle(panel_angle < -Tracking_limit) = -Tracking_limit;
        
        %%%%%%%%%%%%%%%%%%%% POWER CALCS %%%%%%%%%%%%%%%%
        
        % Convert the table to an array
        DATA = importdata('DATA.mat');
        data = table2array(DATA);
        % Extract the header information (angles) from the first row of the table
        angles = data(1, :);
        
        % Remove the header row from the data
        data = data(2:end, :);
        
        % Round the angles in the panel_angles vector
        rounded_panel_angles = round(panel_angle);
        
        % Initialize total_power variable
        total_power = NaN(8760, 1);
        
        % Iterate over each hour in a year
        for  hour = 1:8760
            % Extract the power output values for the current hour from the Excel data
            power_output_hour = data(hour, :);
            
            % Find the index of the column corresponding to the rounded angle in the header
            column_index = find(rounded_panel_angles(hour) == angles, 1);
            
            % If the rounded angle is found, assign the power output value for the corresponding hour and angle
            if ~isempty(column_index)
                total_power(hour) = power_output_hour(column_index);
            else
                % Handle cases where the rounded angle is not in the range -60 to 60
                total_power(hour) = NaN; % Or any other appropriate handling
            end
        end
        
        % Calculate total power output daily
        total_power_daily = reshape(total_power, 24, [])'; % Reshape total_power into a 2D array with 24 columns representing each hour of the day
        total_power_daily = sum(total_power_daily, 2); % Sum the power output for each day
        
        % Calculate total power output yearly
        total_power_yearly = sum(total_power, 'omitnan');
        
        %%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%
        
        % time_vector = 1:numel(panel_angle); % Assuming time vector ranges from 1 to the number of hours
        % 
        % % Plot sun angles
        % plot(time_vector, Sun_angle, 'r', 'LineWidth', 1);
        % 
        % hold on;
        % 
        % % Plot panel angles
        % plot(time_vector, panel_angle, 'b', 'LineWidth', 1);
        % 
        % % Set labels and title
        % xlabel('Time (hours)');
        % ylabel('Angle (degrees)');
        % title('Panel Angles vs Sun Angles');
        % legend('Sun angles', 'Panel angles', 'Location', 'best');
        % xlim([min(time_vector), max(time_vector)]);
        % 
        % % Show grid
        % grid on;
        % 
        % % Show plot
        % hold off;
        % 
        % figure
        % hold on
        % plot(time_vector, Wind_speed_dir, 'b','LineWidth', 1);
        % yline(16, 'r', 'LineWidth', 1);
        % yline(-16, 'r', 'LineWidth', 1);
        % hold off
        % 
        % figure
        % hold on
        % plot(time_vector, Water_level, 'b','LineWidth', 1);
        % yline(0.004, 'r', 'LineWidth', 1);
        % hold off
        % 
        % figure
        % plot((1:length(total_power)),total_power)
end
