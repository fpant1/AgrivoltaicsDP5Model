clc;
clear all; 



test = importGeometry("ellie_cube.STL");

test.rotate(test,10,[0,-5,0],[1,-5,0]);

pdegplot(test, "FaceLabels", "on")





% 
% % Create thermal model
% 
% model = createpde("thermal");
% 
% % Create the geometry and include it in the model
% % figure
% gm = importGeometry(model, "ellie_cube.STL");
% 
% 
% 
% 
% 
% gm.rotate(30,[0,-5,0],[1,-5,0])
% % pdegplot(h,"FaceLabels","on");
% % % hold on
% pdegplot(gm, "FaceLabels","on")
% 
% g=model.Geometry;
% g.Vertices(8,3) = 10;
% [F, V] = g.allDisplayFaces()% V for vertices and F for the facets
% [Ex, Ey, Ez] = g.allDisplayEdges()% Here we have the edges
% pdegplot(g,"FaceLabels","on");
% 


% pdegplot(model);
% 
% % Assign material properties
% thermalProperties(model,"ThermalConductivity",1.5);
% internalHeatSource(model,2e-4)
% 
% % Apply thermal boundaries
% %Bottom layer
% thermalBC(model,"Face",[6],"Temperature",10);
% %Mid Layer
% % thermalBC(thermalmodel,"Face",[2],"Temperature",300);
% %Top Layer
% % thermalBC(thermalmodel,"Face",[2],"Temperature",1);
% 
% % Generate mesh
% mesh = generateMesh(model);
% figure; pdemesh(mesh)
% thermalresults = solve(model)
% 
% pdeplot3D(model,"ColorMapData",thermalresults.Temperature)




% 
% model = createpde;
% importGeometry(model,'Torus.stl');
% mesh = generateMesh(model);
% [p,e,t] = meshToPet(mesh); %convert form of mesh
% figure; pdemesh(mesh);
% figure; pdemesh(p,e,t);

% thermalmodel = createpde('thermal','steadystate');
% importGeometry(thermalmodel,'CombineCubes.stl');
% pdegplot(thermalmodel, 'FaceAlpha',0.5,'FaceLabels','on','CellLabels','on')
% mesh=generateMesh(thermalmodel)
% pdemesh(thermalmodel);
% T1 =90;
% T2 = 20;
% thermalProperties(thermalmodel,'cell',1,'ThermalConductivity',0.15);
% thermalProperties(thermalmodel,'cell',2,'ThermalConductivity',0.3);
% thermalBC(thermalmodel,'Face',5,'Temperature',T1); %Temperature upside
% thermalBC(thermalmodel,'Face',12,'Temperature',T2); %Temperature downside
% results = solve(thermalmodel);
% pdeplot3D(thermalmodel,'ColorMapData',results.Temperature)
% caxis([0 100])