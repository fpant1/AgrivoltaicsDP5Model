classdef PVGISData
    properties
        base_url = 'https://re.jrc.ec.europa.eu/api/v5_1/';
    end
    
    methods
        function obj = PVGISData()
            % Constructor
        end
        
        function tmy_data = getTMYData(obj, latitude, longitude, startYear, endYear)
            % Define the tool name as 'tmy'
            tool_name = 'tmy';

            % Define input parameters
            params = struct();
            params.lat = latitude;
            params.lon = longitude;
            params.usehorizon = 0;
%             params.userhorizon = [0, 10, 20, 30, 40, 15, 25, 5];
%             params.startyear = startYear;
%             params.endyear = endYear;
            params.outputformat = 'json';

            % Construct the query string
            query_string = obj.constructQueryString(params);

            % Construct the full URL
            full_url = [obj.base_url tool_name '?' query_string]

            % Make the GET request
            try
                tmy_data = webread(full_url);
            catch ME
                error('Error retrieving TMY data from PVGIS API: %s', ME.message);
            end
        end
        
        function pv_data = getPVData(obj, latitude, longitude, peakPower)
            % Define the tool name as 'PVcalc'
            tool_name = 'PVcalc';

            % Define input parameters
            params = struct();
            params.lat = latitude;
            params.lon = longitude;
            params.usehorizon = 0;
            params.raddatabase = 'PVGIS-SARAH';
            params.peakpower = peakPower;
            params.pvtechchoice = 'crystSi';
            params.mountingplace = 'free';
            params.loss = 14;
            params.fixed = 0;
            params.angle = 0;
            params.aspect = 0;
            params.optimalinclination = 0;
            params.optimalangles = 0;
            params.inclined_axis = 0;
            params.inclined_optimum = 0;
            params.inclinedaxisangle = 0;
            params.vertical_axis = 0;
            params.vertical_optimum = 0;
            params.verticalaxisangle = 0;
            params.twoaxis = 1;
            params.pvprice = 0;
            params.systemcost = 0;
            params.interest = 0;
            params.lifetime = 25;
            params.outputformat = 'json';
            params.browser = 0;

            % Construct the query string
            query_string = obj.constructQueryString(params);

            % Construct the full URL
            full_url = [obj.base_url tool_name '?' query_string];

            % Make the GET request
            try
                pv_data = webread(full_url);
            catch ME
                error('Error retrieving PV data from PVGIS API: %s', ME.message);
            end
        end

       function temp_data = gettempData(obj, latitude, longitude, showtemperatures)
            % Define the tool name as 'DRcal'
            tool_name = 'DRcalc';

            % Define input parameters
            params = struct();
            params.lat = latitude;
            params.lon = longitude;
            params.month = 0;
            params.showtemperatures = showtemperatures;
            params.outputformat = 'json';

            % Construct the query string
            query_string = obj.constructQueryString(params);

            % Construct the full URL
            full_url = [obj.base_url tool_name '?' query_string];

            % Make the GET request
            try
                temp_data = webread(full_url);
            catch ME
                error('Error retrieving temperature data from PVGIS API: %s', ME.message);
            end
        end
        
        function query_string = constructQueryString(~, params)
            % Construct the query string from input parameters
            query_string = '';
            fields = fieldnames(params);
            for i = 1:length(fields)
                query_string = [query_string fields{i} '=' num2str(params.(fields{i})) '&'];
            end

            % Remove the last '&' character
            query_string = query_string(1:end-1);
        end
    end
end

% % % Define the base URL for PVGIS API
% base_url = 'https://re.jrc.ec.europa.eu/api/v5_1/';
% 
% % Specify the tool name (e.g., 'PVcalc')
% tool_name = 'PVcalc';
% 
% % Define input parameters
% params = struct();
% params.lat = 40.0;  % Example latitude in decimal degrees
% params.lon = -3.0;  % Example longitude in decimal degrees
% params.usehorizon = 0;  % Calculate taking into account shadows from high horizon (1 for "yes")
% % params.userhorizon = '0,10,20,30,40,15,25,5';  % Example horizon heights
% params.raddatabase = 'PVGIS-SARAH';  % Radiation database name
% params.peakpower = 4000;  % Nominal power of the PV system in kW
% params.pvtechchoice = 'crystSi';  % PV technology choice
% params.mountingplace = 'free';  % Type of mounting of PV modules
% params.loss = 20;  % System losses in percent
% params.fixed = 0;  % Calculate a fixed mounted system (1 for "yes")
% params.angle = 0;  % Inclination angle from horizontal plane (if fixed)
% params.aspect = 0;  % Orientation (azimuth) angle (if fixed)
% params.optimalinclination = 0;  % Calculate optimum inclination angle (0 for "no")
% params.optimalangles = 0;  % Calculate optimum inclination and orientation angles (0 for "no")
% params.inclined_axis = 0;  % Calculate a single inclined axis system (0 for "no")
% params.inclined_optimum = 0;  % Calculate optimum angle for inclined axis system (0 for "no")
% params.inclinedaxisangle = 0;  % Inclination angle for inclined axis system (if not calculating optimum)
% params.vertical_axis = 1;  % Calculate a single vertical axis system (0 for "no")
% params.vertical_optimum = 0;  % Calculate optimum angle for vertical axis system (0 for "no")
% params.verticalaxisangle = 0;  % Inclination angle for vertical axis system (if not calculating optimum)
% params.twoaxis = 0;  % Calculate a two axis tracking system (0 for "no")
% params.pvprice = 0;  % Calculate PV electricity price in currency introduced by user for system cost
% params.systemcost = 0;  % Total cost of installing PV system in user's currency (if calculating PV price)
% params.interest = 0;  % Interest rate in %/year (if calculating PV price)
% params.lifetime = 25;  % Expected lifetime of PV system in years
% params.outputformat = 'json';  % Output format as CSV
% params.browser = 0;  % Output as stream (0 for "no")
% 
% % Construct the query string
% query_string = '';
% fields = fieldnames(params);
% for i = 1:length(fields)
%     query_string = [query_string fields{i} '=' num2str(params.(fields{i})) '&'];
% end
% 
% % Remove the last '&' character
% query_string = query_string(1:end-1);
% 
% % Construct the full URL
% full_url = [base_url tool_name '?' query_string];
% 
% % Make the GET request
% try
%     data = webread(full_url);
%     
%     % Display or process the data as needed
%     disp(data);
%     
%     % Write the CSV data to a file (optional)
%     csvwrite('output_data.csv', data);
%     
% catch ME
% %     error('Error retrieving data from PVGIS API: %s', ME.message);
% end