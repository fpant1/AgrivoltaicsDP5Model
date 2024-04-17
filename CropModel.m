classdef CropModel
    properties
        data
        temps
        Tsum    % Cumulative temperature sum (degree days)
        HI      % Harvest index (dimensionless)
        I50A    % Parameter for calculating thermal time for vegetative growth (degree days)
        I50B    % Parameter for calculating thermal time for reproductive growth (degree days)
        Tbase   % Base temperature (°C)
        Topt    % Optimum temperature (°C)
        RUE     % Radiation use efficiency (g/MJ)
        I50maxH % Maximum value for I50A during vegetative growth (degree days)
        I50maxW % Maximum value for I50B during vegetative growth (degree days)
        Theat    % Maximum temperature (°C)
        Text    % Temperature at which thermal time is calculated (°C)
        SCO2    % CO2 concentration scaling factor (dimensionless)
        Swater  % Water stress factor (dimensionless)
        co2

        % ------------------- Calculated -----------
        TT % cumulative temp diff
        mean_temp
        maxtemp
        mintemp
        sum_outputs
        
        I50B_list

        fsolar
        ftemp
        fheat
        fco2
        fco22
        fwater
        bio_cum_array
        bio_rate_array
        awc

        rain
        eto
        irrigation = 0 % mm/day (constant value for each day)
        RZD = 800 % Root Zone Depth (mm)
        DDC = 0.1 % Deep Drainage Coefficient
        RCN = 85 % Runoff Curve Number

    end
    
    methods
        function obj = CropModel(Tsum, HI, I50A, I50B, Tbase, Topt, RUE, I50maxH, I50maxW, Tmax, Text, SCO2, Swater,filename,co2)
            obj.data = readtable(filename);
            obj.temps =obj.data.('T2m');
            obj.rain = obj.data.('rain (mm)');
            obj.eto = obj.data.('Eto');
            obj.Tsum = Tsum;
            obj.HI = HI;
            obj.I50A = I50A;
            obj.I50B = I50B;
            obj.Tbase = Tbase;
            obj.Topt = Topt;
            obj.RUE = RUE;
            obj.I50maxH = I50maxH;
            obj.I50maxW = I50maxW;
            obj.Theat = Tmax;
            obj.Text = Text;
            obj.SCO2 = SCO2;
            obj.Swater = Swater;
            obj.TT = 1;
            obj.co2 = 421; % ppm
            obj.fco22 = co2;
    
        end
        
    
        function obj = calculate_temp_diff(obj,base_temp)
            % Calculates the mean temperature, temperature difference, and cumulative
            % temperature difference over 24-hour intervals, compared to a base temperature.
            % Inputs:
            %   temps: a vector containing the temperature at each hour for a long time
            %   base_temp: the base temperature to compare against
            % Outputs:
            %   mean_temp: a vector containing the mean temperature over 24-hour intervals
            %   temp_diff: a vector containing the temperature difference over 24-hour intervals
            %   cum_diff: a vector containing the cumulative temperature difference over time
            
            % Set the interval length to 24 hours
            interval_length = 24;
            
            % Calculate the number of intervals
            num_intervals = floor(length(obj.temps) / interval_length);
            
            % Initialize the output vectors
            mean_temp = zeros(num_intervals, 1);
            temp_diff = zeros(num_intervals, 1);
            cum_diff = zeros(num_intervals, 1);
            
            % Loop over each interval
            for i = 1:num_intervals
                % Get the indices for the current interval
                start_idx = (i - 1) * interval_length + 1;
                end_idx = i * interval_length;
                
                % Calculate the mean temperature for the current interval
                obj.mean_temp(i) = mean(obj.temps(start_idx:end_idx));
                obj.maxtemp(i) = max(obj.temps(start_idx:end_idx));
                obj.mintemp(i) = min(obj.temps(start_idx:end_idx));
                % Calculate the temperature difference for the current interval
                temp_diff(i) = obj.mean_temp(i) - base_temp;
                
                % If the temperature difference is negative, set it to zero
                if temp_diff(i) < 0
                    temp_diff(i) = 0;
                end
    
                % Update the cumulative temperature difference
                if i == 1
                    cum_diff(i) = temp_diff(i);
                else
                    cum_diff(i) = cum_diff(i-1) + temp_diff(i);
                end
                
                
            end
            obj.TT = cum_diff;

        end






        function obj = fTemp(obj)

%             fTemp = zeros(1, n); % initialize output array to zeros

            for i = 1:length(obj.mean_temp)
                Tmean(i) = obj.mean_temp(i);


                if Tmean(i) < obj.Tbase;
                    obj.ftemp(i) = 0;
                elseif Tmean(i) >= obj.Tbase && Tmean(i) < obj.Topt
                    obj.ftemp(i) = (Tmean(i) - obj.Tbase)/(obj.Topt - obj.Tbase);
                else  Tmean(i) >= obj.Topt;
                    obj.ftemp(i) = 1;
                end
                
            end

        end


        function obj = fHeat(obj)

            for i = 1:length(obj.maxtemp)
                Tmax = obj.maxtemp(i);


                if obj.maxtemp(i) < obj.Theat
                    obj.fheat(i) = 1;
                elseif obj.maxtemp(i) > obj.Text
                    obj.fheat(i) = 0;
                else
                    obj.fheat(i) = 1 - ((obj.maxtemp(i) - obj.Theat) / (obj.Text - obj.Theat));
                end

            end
        end
        
        function obj = fWater(obj)
           
            numberOfDays = length(obj.rain);
            S = (25400 / obj.RCN) - 254; % Calculate potential maximum retention
            Ia = 0.2 * S; % Initial abstraction (initially 20% of S)
            
            % Initialize variables
            obj.awc = zeros(numberOfDays, 1); % Initialize available water content array
            Wi = 0; % Initialize available water in the root zone for the first day
            
            for i = 1:numberOfDays
                % Calculate deep drainage Di using modified DDC (0.1)
                theta_ai_bd = Wi / obj.RZD; % Available soil water content before drainage
                % Assuming theta_m (soil water holding capacity) needs to be calculated or set here
                theta_m = 0.19; % This is a placeholder value and should be calculated or set accordingly
                Di = obj.DDC * (obj.RZD * (theta_ai_bd - theta_m));
                
                % Calculate surface runoff Ri
                Pi = obj.rain(i); % Precipitation for the ith day
                Ti = obj.eto(i); % Transpiration for the ith day (ETO)
                initial_abstraction = Ti + Di; % Initial abstraction = Transpiration + Deep drainage
                if (Pi - initial_abstraction) > 0
                    Ri = ((Pi - initial_abstraction)^2) / (Pi - initial_abstraction + S);
                else
                    Ri = 0; % No runoff if precipitation is less than initial abstraction
                end
                
                % Update available water Wi for the ith day using the water balance equation
                Wi = Wi + Pi + obj.rain - Ti - Di - Ri;
                
                % Make sure Wi does not go negative
                Wi = max(Wi, 0);
                
                % Calculate and store the available water content (awc) for the ith day
                obj.awc(i) = Wi / obj.RZD;

                % Calculate ARID index for the ith day
                PAW = obj.awc(i); % Plant Available Water
                ET_o = obj.eto(i); % Reference evapotranspiration
                ARID = 1 - min(ET_o, 0.096 * PAW) / ET_o;
                
                % Calculate the water stress factor for the ith day
                obj.fwater(i) = 1 - obj.Swater * ARID;


            end
        end

        function obj = calculate_i50_list(obj)
            obj.I50B_list = zeros(size(obj.fheat)); % initialize list of i50 values
            obj.I50B_list(1) = obj.I50B; % set initial i50 value to 0
            
            for i = 2:length(obj.fheat)
                obj.I50B_list(i) = obj.I50B_list(i-1) + (obj.I50maxH*(1-obj.fheat(i-1)));
            end
            for i = 2:length(obj.fwater)
                obj.I50B_list(i) = obj.I50B_list(i-1) + (obj.I50maxW*(1-obj.fwater(i-1)));
            end
        end


        function obj = calculate_fco2_list(obj)
            value = obj.fco22;
            for i = 1:length(obj.mean_temp)
                obj.co2(i) = value ;
               
            end

        end

        function obj = fSolar(obj,growthtosen)
            solar_max = 0.95;
            length(obj.TT);
            for i = 1:length(obj.TT)

                obj.TT(i)
                if i < growthtosen
                    fsolar = solar_max / (1 + exp(-0.01 * (obj.TT(i) - obj.I50A))); % calculate fsolar using equation 1
                else
                    fsolar = solar_max / (1 + exp(0.01 * (obj.TT(i) - (obj.Tsum * obj.I50B_list(i))))); % calculate fsolar using equation 2
                end
                obj.fsolar(i) = fsolar; % save the value of fsolar to obj.fsolar at the correct position
             
            end
            % Create an array of time values (assuming one day intervals)
            time = 1:365;
            
            % Plot the solar values against time
            plot(time, obj.fsolar);
            
            % Set the axis labels and title
            xlabel('Time (days)');
            ylabel('Solar Values');
            title('Solar Values vs. Time');
        end

        function obj = runModel(obj, energyvalues, start_day, end_day)
           
    
            % turn to megajoules/m^2       
            energyvalues = energyvalues/1000000;


            % turns irradiance in to photsynthetically active radiation
            energyvalues = energyvalues*0.48;
            % Loop through each situation
            sum_outputs = zeros(length(energyvalues(:,1)),1);
            bio_cum_array = zeros(length(energyvalues(:,1)),length(start_day:end_day));
            bio_rate_array = zeros(length(energyvalues(:,1)),length(start_day:end_day));
%             output_array_cum = zeros(length(start_day:end_day),length(energyvalues(:,1)));
            for i = 1:length(energyvalues(:,1))
              
                % Initialize the output for this situation
                output = 0;
                
%                 output_array = (length(start_day:end_day));
                % Loop through each day of radiation
                for j = start_day:end_day
                    % Calculate the output for this day using the equation and constants
                    bio_rate = energyvalues(i,j)*obj.fsolar(j)*obj.RUE*obj.co2(j)*obj.ftemp(j)*obj.fheat(j);
                    
                    output = output + bio_rate;
                    bio_cum_array(i,j) = output;
                    bio_rate_array(i,j) = bio_rate;
                end
                output = output*obj.HI*0.01;
                
                                % Store the sum of the outputs for this situation
                sum_outputs(i) = output;
            end
            obj.sum_outputs = sum_outputs;
            obj.bio_cum_array = bio_cum_array;
            obj.bio_rate_array = bio_rate_array;

        end



    end



end