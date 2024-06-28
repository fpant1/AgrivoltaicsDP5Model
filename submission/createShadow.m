classdef createShadow

properties
    vertices = {};
end

methods
    
    function obj = createShadow(module,blocking_object,vectorSol)
        
        obj.vertices = cell(length(blocking_object),1);
        %This function returns a 6x3 matrix with the 6 shadow vertices
        %corresponding to the 8 vertices of the obstacle
        blocking_object;
        %Inputs: 
        % - vertices of ground square (list of 4 3d squares)
        % - vertices of blocking panel
        % - Unit vector of solar location
        
        % for the length/ amount of blocking objects
        for j=1:length(blocking_object)
            %extract 3 vectors to calculate a normal vector to the panel
%             p1 = module.vertices(1,:);
%             p2 = module.vertices(2,:);
%             p3 = module.vertices(3,:);
            p1 = module(1,:);
            p2 = module(2,:);
            p3 = module(3,:);
            
            %does the vector product, calculates the vector "a" normal to the panel
            a = cross((p2-p1),(p3-p1));
    
            s = vectorSol;
            vertice = blocking_object(j).vertices;
            
            % find the ordinate of the plane "d"
            d0 = -(a(1)*p1(1)+a(2)*p1(2)+a(3)*p1(3));
      
            % assign vertices
            for i = 1:size(vertice,1)
                
                % cut the obstacle with the plane of the panel
                % I keep the portion above the panel
                % checks the height of the z coord of a point
                % if it is below the panel, replace it witht the panel z
                % coord
                inter=-(a(1)*vertice(i,1)+a(2)*vertice(i,2)+d0)/a(3);
                if vertice(i,3)<inter
                    vertice(i,3)=inter;
                end
                
                u(:,i) = vertice(i,:) - p1;    % distance between vertex and point on plane
                g = dot(a,u(:,i));             % dot of plane vector distance between plane point and vertex
                h = dot(a,s);                   % must always be positive
                f = g/h;  
                b = f*s;                        % distance * direction of sun ray
                ps(:,i) = vertice(i,:)-b;  % position on plane by combining vertex with its project in s direction
                % if h<0 the vertex is invalid as the ray isntfacing the
                % panel
                if  h<0         % if negative the sun is behind the plane
                    ps(:,i) = NaN;
                end
            end
            vertices = ps.';
            vertices(isnan(vertices)) = [];
            vertices(isinf(vertices)) = [];
            obj.vertices{j} = vertices; % assign vertices to obj.vertices
            
            % convex hull
%             if size(vertices,1)>3 && ~sum(sum(isinf(vertices)))
%                 extIndex = convhull(vertices(:,1),vertices(:,2));
%                 vertices = vertices(extIndex,:);
%                 obj.vertices{j} = vertices;
%             else
%                 obj.vertices{j} = [];
%             end
        end
        
    end


       
    
%     function draw(obj)
%         for i = 1:length(obj.vertices)
%             vertex = obj.vertices{i};
%             if ~isempty(vertex)
%                 hold on
%                 patch(vertex(:,1),vertex(:,2),vertex(:,3),[0.5 0.5 0.5]);
%                 title('cs patch')
%             end
%         end
%         
%     end
    
    
    %obj = vertices projected onto plane
    function coef = coef(obj,panel)
        % computes the area of intersection between the shadow and solar
        % module

        obj.vertices{1};
        
        if ~isempty(obj.vertices{1}) 
            
            %creates shadow polygon (sum of shadows of each obstacle) if
            %there are multiple obstacles
            if length(obj.vertices)>1
                for i=1:length(obj.vertices)-1
                    shadow1 = polyshape(obj.vertices{i}(:,1),obj.vertices{i}(:,2));
                    shadow2 = polyshape(obj.vertices{i+1}(:,1),obj.vertices{i+1}(:,2));
                    poly2 = union(shadow1,shadow2);
                end
            else
                % if there is only 1 obstacle, only 1 shadow
                poly2 = polyshape(obj.vertices{1}(:,1),obj.vertices{1}(:,2));
            end
            % poly1 = ground square
%             poly1 = polyshape(panel.vertices(:,1),panel.vertices(:,2));
            poly1 = polyshape(panel(:,1),panel(:,2));
            polyout = intersect(poly1,poly2);  % area the two regions intersect
            shadowArea = polyout.area;
            panelArea = poly1.area;
            coef = shadowArea/panelArea;  % ratio of shading

            % option = save if the region is shaded or not
        else
            coef = nan;  %
        end
    end

  
end %methods

methods(Static)
 

   function vertices_array = createShadow2(module, blocking_object, vectorSol)
                
        % Inputs: 
        % - module: vertices of ground square (list of 4 3d squares)
        % - blocking_object: rotated_struct with panels and beams at different angles
        % - vectorSol: Unit vector of solar location with tracking angle
    
        p1 = module(1,:);
        p2 = module(2,:);
        p3 = module(3,:);
                
        a = cross((p2-p1), (p3-p1));
        % find the ordinate of the plane "d"
        d0 = -(a(1)*p1(1) + a(2)*p1(2) + a(3)*p1(3));
        
        num_cells = numel(blocking_object);
        vertices_array = cell(length(vectorSol), num_cells);
        
   
        % for every solar position
        for m = 1:length(vectorSol(:,1))
            s = vectorSol(m, 1:3);
            
            x =0;
            % Get the tracking angle from vectorSol and round to nearest integer
            tracking_angle = round(vectorSol(m, 4));
            
    
           
            angle_name = ['field_' num2str(m)];
           
    
            % Assign current_setup based on the found angle_name in blocking_object
            current_setup = blocking_object.(angle_name);
            
                      % Get the number of panels and beams within current_setup
        % Get the number of panels and beams within current_setup
            num_panels = 0;
            if isfield(current_setup, 'panels')
                num_panels = numel(fieldnames(current_setup.panels));
            end
            
            num_beams = 0;
            if isfield(current_setup, 'beams')
                num_beams = numel(fieldnames(current_setup.beams));
            end

            num_pylons = 0;
            if isfield(current_setup, 'pylons')
                num_pylons = numel(fieldnames(current_setup.pylons));
            end


            % Loop through all panels first
            for k = 1:num_panels
                panel_name = ['panel_' num2str(k)];
                vertice = current_setup.panels.(panel_name);
                
                % for every vertex in panel
                for i = 1:size(vertice,1)
                    inter = -(a(1)*vertice(i,1) + a(2)*vertice(i,2) + d0) / a(3);
                    if vertice(i,3) < inter
                        vertice(i,3) = inter;
                    end
    
                    u(:,i) = vertice(i,:) - p1;    
                    g = dot(a, u(:,i));             
                    h = dot(a, s);                   
                    f = g / h;  
                    b = f * s;                        
                    ps(:,i) = vertice(i,:) - b;  
    
                    if h < 0         
                        ps(:,i) = NaN;
                    end
                end
    
                x = x + 1;
                vertices = ps.';
                vertices_array{m, x} = vertices;
            end
            
            % Loop through all beams
            for k = 1:num_beams
                beam_name = ['beam_' num2str(k)];
                vertice = current_setup.beams.(beam_name);
                
                % for every vertex in beam
                for i = 1:size(vertice,1)
                    inter = -(a(1)*vertice(i,1) + a(2)*vertice(i,2) + d0) / a(3);
                    if vertice(i,3) < inter
                        vertice(i,3) = inter;
                    end
    
                    u(:,i) = vertice(i,:) - p1;    
                    g = dot(a, u(:,i));             
                    h = dot(a, s);                   
                    f = g / h;  
                    b = f * s;                        
                    ps(:,i) = vertice(i,:) - b;  
    
                    if h < 0         
                        ps(:,i) = NaN;
                    end
                end
    
                x = x + 1;
                vertices = ps.';
                vertices_array{m, x} = vertices;
            end
             % Loop through all pylons
            for k = 1:num_pylons
                pylon_name = ['pylon_' num2str(k)];
                vertice = current_setup.pylons.(pylon_name);
                
                % for every vertex in pylon
                for i = 1:size(vertice,1)
                    inter = -(a(1)*vertice(i,1) + a(2)*vertice(i,2) + d0) / a(3);
                    if vertice(i,3) < inter
                        vertice(i,3) = inter;
                    end
    
                    u(:,i) = vertice(i,:) - p1;    
                    g = dot(a, u(:,i));             
                    h = dot(a, s);                   
                    f = g / h;  
                    b = f * s;                        
                    ps(:,i) = vertice(i,:) - b;  
    
                    if h < 0         
                        ps(:,i) = NaN;
                    end
                end
    
                x = x + 1;
                vertices = ps.';
                vertices_array{m, x} = vertices;
            end
            
        end
        vertices_array;
        obj.vertices = vertices_array;
    end

    function drawshadow(i,vertices,k)
      
       v = vertices{k,i};
       x = [v(:,1); v(1,1)];
       y = [v(:,2); v(1,2)];
       z = [v(:,3); v(1,3)];

       p = fill3(x,y,z,[0.5 0.5 0.5]);

       p.FaceAlpha = 0.7;

    end

    function direct_shaded = check_shaded(vertices, ground_points)
        num_cells = numel(vertices);
        binary_size = size(vertices);
%         binary_array = cell(binary_size);
        binary_array = zeros(length(ground_points(:,1)),length(vertices(1,:)));
        
        

        shaded_binary_array = ones(length(ground_points(:,1)), 3);
        shaded_binary_array(:,1) = ground_points(:,1);
        shaded_binary_array(:,2) = ground_points(:,2);
     

        direct_shaded = cell(length(vertices(:,1)),1);
        
        

        for j = 1:length(vertices(:,1))
            
            current_array = vertices{j,1};
%             for number of modules rows and columns
            index = 0;
            for col = 1:size(current_array,2)

                for row = 1:size(current_array,1)

                    current_part = current_array{col,row};


                    for module_part = 1:length(current_part)
                        
                        
                        vertexes = current_array{col,row}{1,module_part};
                        
                        
                        for k = 1: length(ground_points(:,1))
                            x = vertexes(:,1);
                            y = vertexes(:,2);
        
                            if isnan(x(1))
                                ground_points(k,3) = 1;
        
                            else 
        
                                [K,A] = convhull(x,y);
                                
                                test_point = ground_points(k,:);
            
                                is_inside = inpolygon(test_point(1), test_point(2), x(K), y(K));
                    
                            
                                if is_inside
        %                         disp(['The test point is inside the shape.' test_point])
                                    ground_points(k,3) = 0;
                            
        %                     else
        % %                         disp('The test point is outside the shape.')
        %                         ground_points(k,3) = 1;
                                end
                            end
                            
        
        
                        end
                            
                    
                    end
                    index = index+1;
                    binary_array(:,index) = ground_points(:,3);
                    randomArrayNew = ones(length(ground_points(:,1)), 1);
                    randomArrayNew(any(binary_array == 0, 2)) = 0;
    
                    shaded_binary_array(:,3) = randomArrayNew;
    
                    direct_shaded{j} = shaded_binary_array;

                end

            end



        end



end


end %classdef
end

