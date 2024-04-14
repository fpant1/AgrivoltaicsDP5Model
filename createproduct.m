classdef createproduct
    properties
        panels
        beams
    end
    
    methods
        function obj = createproduct()
            % Constructor
            obj.panels = {};
            obj.beams = {};
        end
        
        function createPanels(obj, num_panels, panel_offset, y_offset)
            % Create thermal models
            model = createpde("structural");

            % Import the original geometry and include it in the model
            original_panel = importGeometry(model, "testpanel.stl");

            % Extract the vertices and faces of the imported original geometry
            [F, V] = original_panel.allDisplayFaces();

            % Create a new structure array to store the modified panels
            obj.panels = cell(1, num_panels + 1);

            % Create the original panel
            panel_geometry = struct('vertices', V, 'faces', F);
            panel_geometry.vertices(:, 1) = panel_geometry.vertices(:, 1) + panel_offset(1); % Offset in x
            panel_geometry.vertices(:, 2) = panel_geometry.vertices(:, 2) + panel_offset(2); % Offset in y
            panel_geometry.vertices(:, 3) = panel_geometry.vertices(:, 3) + panel_offset(3); % Offset in z
            obj.panels{1} = panel_geometry;

            % Create additional modified panels
            for i = 2:(num_panels + 1)
                obj.panels{i} = panel_geometry;
                obj.panels{i}.vertices(:, 2) = obj.panels{i}.vertices(:, 2) + (i - 1) * y_offset;
            end
        end
        
        function createBeams(obj, num_beams, x_offsets, z_offsets)
            % Create thermal models
            model1 = createpde("structural");

            % Import a new geometry (boxsectiontest.stl) and include it in the model
            beam = importGeometry(model1, "boxsectiontest.stl");

            % Extract the vertices and faces of the imported beam geometry
            [F_beam, V_beam] = beam.allDisplayFaces();

            % Create a new structure array to store the modified beams
            obj.beams = cell(1, num_beams);

            % Create modified beams
            for i = 1:num_beams
                obj.beams{i} = struct('vertices', V_beam, 'faces', F_beam);
                obj.beams{i}.vertices(:, 1) = obj.beams{i}.vertices(:, 1) + x_offsets(i);
                obj.beams{i}.vertices(:, 3) = obj.beams{i}.vertices(:, 3) + z_offsets(i);
            end
        end
        
        function draw(obj)
            % Plot the panels and beams
            figure;
            hold on;
            for i = 1:(length(obj.panels) + length(obj.beams))
                if i <= length(obj.panels)
                    patch('Vertices', obj.panels{i}.vertices, ...
                          'Faces', obj.panels{i}.faces, ...
                          'FaceColor', [1.0, 0.7, 0.0], ...
                          'EdgeColor', [0.1, 0.1, 0.1], ...
                          'FaceAlpha', 0.7);
                else
                    j = i - length(obj.panels);
                    patch('Vertices', obj.beams{j}.vertices, ...
                          'Faces', obj.beams{j}.faces, ...
                          'FaceColor', [0.7, 0.0, 1.0], ...
                          'EdgeColor', [0.1, 0.1, 0.1], ...
                          'FaceAlpha', 0.7);
                end
            end

            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('Original Panel, Modified Panels, and Modified Beams');
            axis equal;
            grid on;
            legend('Original Panel', 'Modified Panels', 'Modified Beam 1', 'Modified Beam 2');

            hold off;
        end
    end
end
