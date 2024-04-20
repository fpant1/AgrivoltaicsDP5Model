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
electrical_limit = 1.5e10;



% Calculate the maximum values for agricultural and electrical outputs
max_agricultural_output = max(combinedResults.final_agricultural_output);
max_electrical_output = electrical_limit;

% Calculate the percentage change for each entry in the agricultural output
combinedResults.percent_change_agricultural_yield = (combinedResults.final_agricultural_output / max_agricultural_output);

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
xlin = linspace(min(X), max(X), 500);
ylin = linspace(min(Y), max(Y), 500);
[Xgrid, Ygrid] = meshgrid(xlin, ylin);

% Interpolate the data for a smoother surface
Z1grid = griddata(X, Y, Z1, Xgrid, Ygrid, 'cubic');
Z2grid = griddata(X, Y, Z2, Xgrid, Ygrid, 'cubic');
Z3grid = griddata(X, Y, Z3, Xgrid, Ygrid, 'cubic');


nan_grid = Z1grid > electrical_limit;

Z1grid(nan_grid) = nan;
Z2grid(nan_grid) = nan;
Z3grid(nan_grid) = nan;

% Create the first surface plot
figure;
surf(Xgrid, Ygrid, Z1grid);
shading interp;
hold on;

% Add contour lines to the first surface plot
contour3(Xgrid, Ygrid, Z1grid, 'k');

% Customize the colormap
% colormap parula;

% Add labels and legend
xlabel('GCR');
ylabel('Row spacing (m)');
zlabel('Crop yield (t/ha)');
title('Agrivoltaics Yield 1 as a Function of GCR and tilt');
colorbar;
legend('Electical Yield (kWh)');
view(2);
set(gcf, 'WindowState', 'maximized');

% Release the hold on the current figure
hold off;

% Create the second surface plot
figure;
surf(Xgrid, Ygrid, Z2grid);
shading interp;
hold on;

% Add contour lines to the second surface plot
contour3(Xgrid, Ygrid, Z2grid, 'k');

% % Create a grid with Y = 40 for the specific contour line
% Ygrid40 = 40 * ones(size(Xgrid));
% 
% % Interpolate the data for the specific contour line
% Z2grid40 = griddata(X, Y, Z2, Xgrid, Ygrid40, 'cubic');
% 
% % Add the specific contour line at Y = 40
% contour3(Xgrid, Ygrid40, Z2grid40, 'r', 'LineWidth', 2);

% Customize the colormap
% colormap parula;

% Add labels and legend
xlabel('Number of Panels');
ylabel('Row spacing (m)');
zlabel('Agricultural Yield');
title('Agricultural yield as a Function of row spacing and number of panels');
colorbar;
legend('Electrical Yield (kWh/m^2/day)');
view(2);
set(gcf, 'WindowState', 'maximized');

% Release the hold on the current figure
hold off;


% Plot for Land Efficiency
figure;
surf(Xgrid, Ygrid, Z3grid);
shading interp;
hold on;

% Add contour lines to the land efficiency surface plot
contour3(Xgrid, Ygrid, Z3grid, 'k');

% Add labels and legend
xlabel('Number of Panels');
ylabel('Row spacing (m)');
zlabel('Land Efficiency (%)');
title('Land Efficiency as a Function of Number of Panels and Row Spacing');
colorbar;
legend('Land Efficiency');
view(2);
set(gcf, 'WindowState', 'maximized');

% Release the hold on the current figure
hold off;