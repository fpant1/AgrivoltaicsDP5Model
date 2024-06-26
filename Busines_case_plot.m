clear;
clc;

% Load the Thursday night sensitivity study results
load('thursdaynight_sensitivity_study_results.mat');
% The table loaded is named "results", renaming for clarity
thursdayResults = results; % Rename the variable for clarity

% Load the Saturday sensitivity study results
load('saturday_sensitivity_study_results.mat');
% The table loaded is also named "results", renaming for clarity
saturdayResults = results; % Rename the variable for clarity

% Vertically concatenate the Saturday results to the Thursday results
combinedResults = [thursdayResults; saturdayResults];

% Multiply the final_site_kWh_output by 7.36613965e-4
combinedResults.final_site_kWh_output = combinedResults.final_site_kWh_output * 7.36613965e-4;

% Store the columns by their heading names from the combined results
total_panels = combinedResults.total_panels;
num_modules = combinedResults.num_modules;
num_rows = combinedResults.num_rows;
num_panels = combinedResults.num_panels;
row_spacing = combinedResults.row_spacing;
gcr = combinedResults.gcr;
installed_capacity_kW = combinedResults.installed_capacity_kW;
final_site_kWh_output = combinedResults.final_site_kWh_output;
final_agricultural_output = combinedResults.final_agricultural_output;

% Proceed with your previous code logic after combining the tables
electrical_limit = 1.5e10 * 7.36613965e-4;



% Calculate the maximum values for agricultural and electrical outputs
max_agricultural_output = max(combinedResults.final_agricultural_output);
max_electrical_output = electrical_limit;

no_structure_yield =22.3;


% Calculate the percentage change for each entry in the agricultural output
combinedResults.percent_change_agricultural_yield = (combinedResults.final_agricultural_output / no_structure_yield);

% Calculate the percentage change for each entry in the electrical output
combinedResults.percent_change_electrical_yield = (combinedResults.final_site_kWh_output / max_electrical_output);

% Sum the two percentage changes to get the land efficiency
combinedResults.land_efficiency = combinedResults.percent_change_agricultural_yield + combinedResults.percent_change_electrical_yield;

% Store the updated results
save('complete_thursdaynight_sensitivity_study_results.mat', 'results');





% Create the data
X = num_panels; % GCR
Y = row_spacing; % tilt
Z1 = final_site_kWh_output; % Yield 1
Z2 = final_agricultural_output; % Yield 2
Z3 = combinedResults.land_efficiency; % Land Efficiency



% Y = Y.*-1;

% Create a grid for the surface plot
xlin = linspace(min(X), max(X), 2000);
ylin = linspace(min(Y), max(Y), 2000);
[Xgrid, Ygrid] = meshgrid(xlin, ylin);

% Interpolate the data for a smoother surface
Z1grid = griddata(X, Y, Z1, Xgrid, Ygrid, 'cubic');
Z2grid = griddata(X, Y, Z2, Xgrid, Ygrid, 'cubic');
Z3grid = griddata(X, Y, Z3, Xgrid, Ygrid, 'cubic');
% Interpolate the percentage change data for agricultural yield over the grid
AgriPercentChangeGrid = griddata(X, Y, combinedResults.percent_change_agricultural_yield, Xgrid, Ygrid, 'cubic');
ElecPercentChangeGrid = griddata(X, Y, combinedResults.percent_change_electrical_yield, Xgrid, Ygrid, 'cubic');

nan_grid = Z1grid > electrical_limit;

Z1grid(nan_grid) = nan;
Z2grid(nan_grid) = nan;
Z3grid(nan_grid) = nan;
AgriPercentChangeGrid(nan_grid) = nan;
ElecPercentChangeGrid(nan_grid) = nan;

contourline70 = zeros(size(AgriPercentChangeGrid,1),size( AgriPercentChangeGrid,2));
elecline70 = zeros(size(ElecPercentChangeGrid,1),size( ElecPercentChangeGrid,2));


% contourline80 = zeros(size(AgriPercentChangeGrid,1),size( AgriPercentChangeGrid,2));
% contourline60 = zeros(size(AgriPercentChangeGrid,1),size( AgriPercentChangeGrid,2));
for i = 1:size(AgriPercentChangeGrid,1)
    for j = 1: size( AgriPercentChangeGrid,2)
        if AgriPercentChangeGrid(i,j) < 0.7001 && AgriPercentChangeGrid(i,j) > 0.6999
            contourline70(i,j) = 1;
        end
        if AgriPercentChangeGrid(i,j) < 0.6001 && AgriPercentChangeGrid(i,j) > 0.5999
            contourline60(i,j) = 1;
        end
        if AgriPercentChangeGrid(i,j) < 0.8001 && AgriPercentChangeGrid(i,j) > 0.7999
            contourline80(i,j) = 1;
        end
        if ElecPercentChangeGrid(i,j) < 0.7003 && ElecPercentChangeGrid(i,j) > 0.6993
            elecline70(i,j) = 1;
        end
%         if ElecPercentChangeGrid(i,j) < 0.6001 && AgriPercentChangeGrid(i,j) > 0.5999
%             contourline60(i,j) = 1;
%         end
%         if ElecPercentChangeGrid(i,j) < 0.8001 && ElecPercentChangeGrid(i,j) > 0.7999
%             contourline80(i,j) = 1;
%         end
    end

end
plotcontourline70 = contourline70.*Z3grid;
% plotcontourline60 = contourline60.*Z3grid;
% plotcontourline80 = contourline80.*Z3grid;

plotcontourline70(plotcontourline70==0)= nan;
% plotcontourline60(plotcontourline60==0)= nan;
% plotcontourline80(plotcontourline80==0)= nan;

plotelecline70 = elecline70.*Z3grid;
plotelecline70(plotelecline70==0) = nan;


% Save the updated combined results to a CSV file
writetable(combinedResults, 'adjust_max_crop_sensitivity_study_results.csv');

% % Create the first surface plot
% figure;
% surf(Xgrid, Ygrid, Z1grid);
% shading interp;
% hold on;
% 
% % Add contour lines to the first surface plot
% contour3(Xgrid, Ygrid, Z1grid, 'k');
% 
% % Customize the colormap
% % colormap parula;
% % Set the font size for axes labels and title
% ax = gca; % Current axes handle
% ax.FontSize = 22; % Set the font size for axis tick labels
% ax.XTick = [3, 4, 5, 6, 7]; % Set specific x-axis ticks
% % Add labels and legend
% xlabel('Number of Panels', 'FontSize', 22);
% ylabel('Row spacing (m)', 'FontSize', 22);
% zlabel('Number of Panels', 'FontSize', 22);
% title('Electrical Yield as a Function of Row Spaciing and Number of Panels', 'FontSize', 22);
% cb = colorbar;
% cb.Label.String = 'Electical Yield (kWh)';
% cb.Label.FontSize = 22; % Set the font size for the color bar label
% cb.FontSize = 22; 
% % legend('Electical Yield (kWh)', 'FontSize', 22);
% view(2);
% set(gcf, 'WindowState', 'maximized');
% 
% % Release the hold on the current figure
% hold off;



% ---------------------------------------------------------

% % Create the second surface plot
% figure;
% surf(Xgrid, Ygrid, Z2grid);
% shading interp;
% hold on;
% 
% % Add contour lines to the second surface plot
% contour3(Xgrid, Ygrid, Z2grid, 'k');
% 
% % % Create a grid with Y = 40 for the specific contour line
% % Ygrid40 = 40 * ones(size(Xgrid));
% % 
% % % Interpolate the data for the specific contour line
% % Z2grid40 = griddata(X, Y, Z2, Xgrid, Ygrid40, 'cubic');
% % 
% % % Add the specific contour line at Y = 40
% % contour3(Xgrid, Ygrid40, Z2grid40, 'r', 'LineWidth', 2);
% 
% % Customize the colormap
% % colormap parula;
% % Set the font size for axes labels and title
% ax = gca; % Current axes handle
% ax.FontSize = 22; % Set the font size for axis tick labels
% ax.XTick = [3, 4, 5, 6, 7]; % Set specific x-axis ticks
% % Add labels and legend
% xlabel('Number of Panels', 'FontSize', 22);
% ylabel('Row spacing (m)', 'FontSize', 22);
% zlabel('Agricultural Yield', 'FontSize', 22);
% title('Agricultural yield as a Function of row spacing and number of panels', 'FontSize', 22);
% cb = colorbar;
% cb.Label.String = 'Agricultural Yield (t/ha)';
% cb.Label.FontSize = 22; % Set the font size for the color bar label
% cb.FontSize = 22; 
% % legend('Agricultural Yield (t/ha)', 'FontSize', 22);
% view(2);
% set(gcf, 'WindowState', 'maximized');
% 
% % Release the hold on the current figure
% hold off;
% Plot for Land Efficiency
figure;
plot3(Xgrid, Ygrid, plotcontourline70, 'ro', 'MarkerSize', 2);
hold on% Red circles for 70% agricultural yield
plot3(Xgrid, Ygrid, plotelecline70, 'bo', 'MarkerSize', 2); % Blue circles for 70% electrical yield

surf(Xgrid, Ygrid, Z3grid); % Plotting the 3D surface
shading interp; % Interpolating shading for surface


% Add contour lines to the land efficiency surface plot


% Plot specific contour lines representing specific yields

contour3(Xgrid, Ygrid, Z3grid, 'k'); % Black contour lines
% Set the font size for axes labels and title
ax = gca; % Current axes handle
ax.FontSize = 26; % Set the font size for axis tick labels
ax.XTick = [3, 4, 5, 6, 7]; % Set specific x-axis ticks

% Add labels and legend
xlabel('Number of Panels', 'FontSize', 32);
ylabel('Row spacing (m)', 'FontSize', 32);
zlabel('Land Efficiency Ratio', 'FontSize', 32);
% title('Land Efficiency as a Function of Number of Panels and Row Spacing', 'FontSize', 26);

% Adding a colorbar for visualising the land efficiency ratio
cb = colorbar;
cb.Label.String = 'Land Efficiency Ratio';
cb.Label.FontSize = 26; % Set the font size for the color bar label
cb.FontSize = 26; % Set the font size for the color bar ticks

% Adding legend to identify plots
legend('70% Agricultural Yield', '70% Electrical Yield', 'FontSize', 32);

view(2); % Change the view to 2D
set(gcf, 'WindowState', 'maximized'); % Maximize the figure window

% Release the hold on the current figure
hold off;


% 
% % Plot for Agricultural Percentage Change
% figure;
% surf(Xgrid, Ygrid, AgriPercentChangeGrid);
% shading interp;  % This makes the surface look smoother
% hold on;
% 
% % Add contour lines to the agricultural percentage change surface plot
% contour3(Xgrid, Ygrid, AgriPercentChangeGrid, 'k');
% % Set the font size for axes labels and title
% ax = gca; % Current axes handle
% ax.FontSize = 22; % Set the font size for axis tick labels
% ax.XTick = [3, 4, 5, 6, 7]; % Set specific x-axis ticks
% % Add labels and legend
% xlabel('Number of Panels', 'FontSize', 22);
% ylabel('Row spacing (m)', 'FontSize', 22);
% zlabel('Agricultural Yield Change (%)', 'FontSize', 22);
% title('Relative Agricultural Yield as a Function Row Spacing and Number of Panels',  'FontSize', 22);
% cb = colorbar;  % Adds a color bar to the right of the plot to indicate value ranges
% % legend('Agricultural Yield Change', 'FontSize', 22);
% % cb = colorbar;  % Adds a color bar to the right of the plot to indicate value ranges
% cb.Label.String = 'Relative Agricultural Yield';
% cb.Label.FontSize = 22; % Set the font size for the color bar label
% cb.FontSize = 22; % Set the font size for the color bar ticks
% view(2);  % Adjust the view to be top-down
% set(gcf, 'WindowState', 'maximized');  % Make the figure fullscreen
% 
% % Release the hold on the current figure
% hold off;
% 
