classdef Irradiance
    properties
        Gd
        Gb
        globalir
        data
        time
        irradvalues
        powervalues
        energyvalues
        fullvalue
        globalhoriz
        multiplied_array
        sum_globalhoriz
        maxenergyvalues
        dailyshadingfactor
        monthlyshadingfactor
        totalshadingfactorpoints
        totaltotalshadingfactor 

  
    end
    
    methods
        function obj = Irradiance( filename)
            obj.data = filename.outputs.tmy_hourly;
            obj.data = struct2table(obj.data);
            
            obj.Gb = obj.data.('Gb_n_');
            obj.Gd = obj.data.('Gd_h_');
            obj.time = obj.data.('time_UTC_');
            obj.globalhoriz = obj.data.('G_h_');
            obj.Gb = obj.globalhoriz - obj.Gd;
            obj.irradvalues = 1;
            obj.multiplied_array = obj.globalhoriz *3600;
            obj.sum_globalhoriz = sum(obj.multiplied_array);
            

        end


        function obj = calcdata(obj,shade_array,diffuse)
            
%             shade_array = shade_array.shade_array;%.direct_shaded;
            direct_irad = obj.Gb;
            diffuse_irad = obj.Gd;

            
            for i = 1 : length(direct_irad)
                shade = num2cell(shade_array{i});

               
                shade(:, 4) = num2cell([shade{:, 3}] .* direct_irad(i))';
                shade(:, 5) = num2cell([diffuse_irad(i)].* (1-diffuse));
                time1 = obj.time{i};
                X = convertCharsToStrings(time1);
                shade{1,6} = X;

%                 numRows = size(shade(:,2));
%                 newColumn = repmat({obj.time{i}}, numRows, 1);
%                 shade = horzcat(newColumn, shade,1);
                shade_array{i} = shade;
                


                
            end
            obj.irradvalues = shade_array;
        end

        function obj = dailyenergy(obj)
            irradvalues = obj.irradvalues;
            powervalues = zeros(length(obj.irradvalues),length(obj.irradvalues{1}(:,1)));
            globalirad = obj.globalir;

            globaldiffuse = obj.Gd;
            globaldirect = obj.Gb;
            globalval = zeros(length(globalirad));
            for i = 1: 8760 
                
                for j = 1: length(obj.irradvalues{i}(:,1))

                    powervalue = irradvalues{i}{j,4}+ irradvalues{i}{j,5};

                    powervalues(i,j) = powervalue*3600;
                    
                    globalval(i) = (globaldiffuse(i)+ globaldirect(i))*3600;
                end
                
            end
           
            obj.powervalues = powervalues;
            energy = obj.powervalues;
            runtime = round(length(powervalues(:,1)) ./ 24);
         
            energyvalues = zeros(length(powervalues(1,:)),runtime);

            
            for k = 1:length(powervalues(1,:)) 
                
                
                for l = 1:runtime
                    
                    start_index = (l-1)*24 + 1;
                    end_index = l*24;
                    energyvalues(k,l) = sum(energy(start_index:end_index,k));
                    globalenergy(l) = sum(globalval(start_index:end_index));

                end
            end
                % Calculate maxenergyvalues
            maxenergyvalues = zeros(1, 365);
            for day = 1:365
                start_index = (day-1)*24 + 1;
                end_index = day*24;
                maxenergyvalues(day) = sum(obj.multiplied_array(start_index:end_index));
            end

                % Calculate daily shading factor
            dailyshadingfactor = zeros(size(obj.energyvalues));
            for i = 1:365
                dailyshadingfactor(:,i) = energyvalues(:,i) / maxenergyvalues(i);
            end

               % Calculate monthly shading factor
            monthlyshadingfactor = zeros(size(dailyshadingfactor, 1), 12);
            days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31];
            day_count = 1;
            for month = 1:12
                start_day = day_count;
                end_day = start_day + days_in_month(month) - 1;
                monthlyshadingfactor(:, month) = mean(dailyshadingfactor(:, start_day:end_day), 2);
                day_count = end_day + 1;
            end

            totalshadingfactor1 = mean(dailyshadingfactor, 2);
            totalshadingfactor2 = mean(totalshadingfactor1,1)


            obj.globalir = globalenergy;
            obj.energyvalues = energyvalues;
            obj.maxenergyvalues = maxenergyvalues;
            obj.dailyshadingfactor = dailyshadingfactor;
            obj.monthlyshadingfactor = monthlyshadingfactor;
            obj.totalshadingfactorpoints = totalshadingfactor1;
            obj.totaltotalshadingfactor = totalshadingfactor2;
          


        end

        function obj = shadingfactorplot(obj)

%             X = obj.irradvalues{1, 1}(:, 1);
%             Y = obj.irradvalues{1, 1}(:, 2);
            X = cellfun(@double, obj.irradvalues{1, 1}(:, 1));
            Y = cellfun(@double, obj.irradvalues{1, 1}(:, 2));
            Z = obj.totalshadingfactorpoints(:);
            % Create a table with cell array variables
            tbl = table(X, Y, Z, 'VariableNames', {'X', 'Y', 'Z'});

        
            % Create a new figure for the heatmap
            figure;
            h = heatmap(tbl, 'X', 'Y', 'ColorVariable', 'Z');
        
            % Customize the heatmap as needed
            % For example:
            % h.Title = 'Shading Factor Heatmap';
            % h.XLabel = 'X';
            % h.YLabel = 'Y';

        end
    end
end
