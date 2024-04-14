classdef createModule < handle

%you pass:
% width:    width of the module
% height:   height of the module
% position: 1 vector with components x,y,z
% az:       azimuth angle in deg
% tilt:     tilt angle in deg

properties
    % size
    height
    width
    % vertices
    vertices;
    % angles
    azimuth
    tilt
    %position
    position = [0, 0, 0];
end
  

methods 
    
% constructor
function obj = createModule(module)
    obj.width = module.width;
    obj.height = module.height;
    obj.azimuth = module.azimuth;
    obj.tilt = module.tilt;
    
%     switch obj.tilt
%         case 0
%             obj.tilt = 0.001;
%         case 90
%             obj.tilt = 89.9;
%     end
%     
    %creates each corner of the module
    vertices = zeros(4,3);
    %bottom left corner
    vertices(1,1) = 0;
    %bottom right corner
    vertices(2,2) = obj.width;
    %top right corner
    vertices(3,1) = -obj.height;
    vertices(3,2) = obj.width;
    %top left corner
    vertices(4,1) = -obj.height;
    
    %1- rotate the module with tilt angle (around y axis)
    %2- rotate the module with azimuth angle (around z axis)
    
    vertices = vertices.';
    for i=1:4
    vertices(:,i) = roty(obj.tilt)*vertices(:,i);
    vertices(:,i) = rotz(obj.azimuth)*vertices(:,i);
    end
    vertices = vertices.';
    
    obj.vertices = vertices;
end

% shows 3D representation of the module
function draw(obj)
       v = obj.vertices;
       x = [v(:,1); v(1,1)];
       y = [v(:,2); v(1,2)];
       z = [v(:,3); v(1,3)];

       fill3(x,y,z,'b')
       

%        f = [1 2 3 4];
%        patch('Faces',f,'Vertices',v,'FaceColor','b','FaceAlpha',0.6);
%        axis equal
       v;
end

function move(obj,displacement)
    %moves the module
    obj.position = obj.position+displacement;
    obj.vertices=obj.vertices+displacement;
end

function rotateFromOrigin(obj,angle)
    obj.vertices = obj.vertices*rotz(angle);
end

end
end