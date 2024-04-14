classdef createObstacle < handle

    
% This class creates an obstacle object.
% The object is a rectangular prism.

%INPUTS:   - Size x,y,z
%          - Position x,y,z
%          - Azimuth angle
%          - Tilt angle

%OUTPUT:   - Obstacle object

properties
        %dimentions
        dx = 0
        dy = 0
        dz = 0
        %azimuth angle
        azimuth = 0
        %tilt angle
        tilt = 0
        %displacement vector
        position = zeros(3,1);
        %vertices matrix
        vertices = zeros(3,8)
    end
    
methods
    %constructor
    function obj = createObstacle(obstacleInput)

                %Asigna valores a las propiedades
                obj.dx = obstacleInput.dimX;
                obj.dy = obstacleInput.dimY;
                obj.dz = obstacleInput.dimZ;
                obj.position = obstacleInput.position;
                obj.azimuth = obstacleInput.azimuth;
                obj.tilt = obstacleInput.tilt;
                
                % assign short names
                dx = obj.dx;
                dy = obj.dy;
                dz = obj.dz;

                % definition of vertices for the rectangular prism
                obj.vertices(:,1) = [0 0 0];
                obj.vertices(:,2) = [dx 0 0];
                obj.vertices(:,3) = [dx dy 0];
                obj.vertices(:,4) = [0 dy 0];
                obj.vertices(:,5) = [0 0 dz];
                obj.vertices(:,6) = [dx 0 dz];
                obj.vertices(:,7) = [dx dy dz];
                obj.vertices(:,8) = [0 dy dz];

                az = -obj.azimuth*pi/180;
                tilt = -obj.tilt*pi/180;

                %create azimuth rotation matrix
                matrizAzimut = [cos(az) -sin(az) 0; sin(az) cos(az) 0; 0 0 1];
                %crea matriz de rotacion inclinacion
                matrizInclinar = [cos(tilt) 0 sin(tilt); 0 1 0; -sin(tilt) 0 cos(tilt)];

                %rota y desplaza la matriz con los vertices
                for i=1:8

                    % rotation
                    obj.vertices(:,i) = matrizInclinar*obj.vertices(:,i);
                    obj.vertices(:,i) = matrizAzimut*obj.vertices(:,i);

                    % displacement
                    obj.vertices(:,i) = obj.vertices(:,i) + obj.position.';

                end
        obj.vertices = obj.vertices.';
    end
    
    % shows a 3D representation of the obstacle
    function draw(obj)
        for i=1:length(obj)
        hold on
        % transpose for using in patch function
        v = obj(i).vertices;
        % vertex order for connections
        f = [1 2 3 4; 5 6 7 8 ; 1 5 8 4; 2 6 7 3; 1 5 6 2; 4 8 7 3];  
        patch('Faces',f,'Vertices',v,'FaceColor',[0.3 0.3 0.7],'FaceAlpha',0.5);
        end
        hold off
    end
     
end %methods
end %clasedef

