classdef createTable
    %creates a table with shading coefficient results
    
    properties
        coefficients
        azRange
        elRange
        module
        obstacle
    end
    
    methods
        function obj = createTable(module,obstacle,resolution)
            azRange = -180:resolution:180;
            elRange = 0:resolution:90;
            all_tables = struct()  % create structure to save tables
            table = zeros(length(elRange),length(azRange));
            [az,el] = meshgrid(azRange,elRange);
            for a = 1:2
                for i=1:size(el,1)
                    for j=1:size(az,2)
                        [sunX, sunY, sunZ] = sph2cart(-az(i,j)*pi/180,el(i,j)*pi/180,1);%SPH2CART - Transform spherical to Cartesian coordinates.
    %                     z = r .* sin(elev);
    %                     rcoselev = r .* cos(elev);
    %                     x = rcoselev .* cos(az);
    %                     y = rcoselev .* sin(az);
                        sun = [sunX sunY sunZ];
                        square = module(4*a-3:4*a, :);
%                         for i = obstacles
%                          
                        sombra = createShadow(square,obstacle,sun); % in create shadow extract the projected points for every combination of az/el at given resolution
                        table(i,j) = sombra.coef(square); % in create shadow, extract the shading coefficient for that setup
                        all_tables.(['table_',a]) = table;
                    end
                end
            end
            AZrange = -180:1:180;
            ELrange = 0:1:90;
            if resolution~=1
            % interpolates values
                [AZ,EL] = meshgrid(AZrange,ELrange);
                obj.coefficients = interp2(az,el,table,AZ,EL);
            else
                obj.coefficients = all_tables
            end
            
            all_tables
            
            % fixes negative values
            for i=1:length(AZrange)
                if AZrange(i)<0
                    AZrange(i) = AZrange(i) + 360;
                end
            end
            obj.azRange = AZrange
            obj.elRange = ELrange;
            obj.module = module;
            obj.obstacle = obstacle;
            
        end
        
        function coeff = coeff(obj,elevation,azimut)
            % returns the coefficient for a given elevation and azimut
            % solar angles
            elevation = round(elevation);
            azimut = round(azimut);
            
            el = obj.elRange;
            az = obj.azRange;
            
            %fixes bug
            if azimut==360
                azimut=0;
            elseif azimut<0
                azimut = azimut+360;
            end
            
            if elevation>=1
                elevIndex = find(el==elevation);
                azIndex = find(az==azimut);
                coeff = obj.coefficients(elevIndex,azIndex(1));
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
        
        function dif = diffuse(obj)
            %integra la tabla para calcular la componente difusa
            den = 0;
            nom = 0;

            for i=1:length(obj.elRange)
                for j=1:length(obj.azRange)
                    AOI = pvl_getaoi(obj.module.tilt,obj.module.azimuth,90-obj.elRange(i),obj.azRange(j));
                    if ~isnan(obj.coefficients(i,j))
                        nom = nom+obj.coefficients(i,j)*cosd(AOI)*cosd(obj.elRange(i));
                        den = den+cosd(AOI)*cosd(obj.elRange(i));
                    end
                end
            end
            dif = nom/den;
        end
        
        function TT = directFactor(obj,year,resolution,UTC,latitude,longitude,altitude)
            % returns the shading factors for a given time and location
            
            t1 = datetime(year, 1, 1, 0, 0, 0);                             % start time
            t2 = datetime(year, 12, 31, 24, 0, 0);                          % end time
            t = t1:minutes(resolution):t2;                                  % time resolution 
            time_num = datenum(t);

            time = pvl_maketimestruct(time_num,UTC);                        % time structure
            location = pvl_makelocationstruct(latitude,longitude,altitude); % Create a structure to define a site location

            % sun position angles
            [azimut, elevation, elAparent]=pvl_spa(time, location);%Calculates the position of the sun given time, location, and optionally pressure and temperature
            fB = zeros(length(elevation),1);
            
            % pass azimuth and apparent elevation into coeff function
            for i=1:length(elevation)
                fB(i) = obj.coeff(elAparent(i),azimut(i));
            end
            TT = timetable(fB,'RowTimes',t);
        end
        
    end
end

