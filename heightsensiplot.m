% clc;
% clear;
% 
% % Load data
% load('height_validation3.mat');

% Extract necessary variables from results
height = results.Var11;  % Assuming this represents different heights
row_spacings = results.row_spacing;
final_site_kWh_output = results.final_site_kWh_output;
final_agricultural_output = results.final_agricultural_output;

% Constants and calculations
electrical_limit = 1.5e10 * 7.36613965e-4;
max_agricultural_output = max(final_agricultural_output)+0.18;
% max_electrical_output = electrical_limit;
max_electrical_output = max(final_site_kWh_output);

% Calculate percentage changes relative to maximum outputs
percent_change_agricultural_yield = (final_agricultural_output / max_agricultural_output) * 100; % Convert to percentage
percent_change_electrical_yield = (final_site_kWh_output / max_electrical_output) * 100; % Convert to percentage

% Prepare data for plotting
X = row_spacings; % Switched X and Y
Y = height;
Z1 = final_site_kWh_output;  % Final site kWh output
Z2 = percent_change_agricultural_yield;  % Agricultural yield percentage change
Z3 = percent_change_electrical_yield;  % Electrical yield percentage change

% Create grid for surface plots
xlin = linspace(min(X), max(X), 200);
ylin = linspace(min(Y), max(Y), 200);
[Xgrid, Ygrid] = meshgrid(xlin, ylin);

% Interpolate the data
Z1grid = griddata(X, Y, Z1, Xgrid, Ygrid, 'cubic');
Z2grid = griddata(X, Y, Z2, Xgrid, Ygrid, 'cubic');
Z3grid = griddata(X, Y, Z3, Xgrid, Ygrid, 'cubic');

% Plotting the final site kWh output with contours
figure;
surf(Xgrid, Ygrid, Z1grid); % Create the surface plot
shading interp; % Apply interpolated shading for smooth color transitions
hold on;  % Hold the plot to add contour lines
contour3(Xgrid, Ygrid, Z1grid, 20, 'k');  % Add black contour lines for depth
xlabel('Row Spacing (m)', 'FontSize', 40); % Set x-axis label
ylabel('Height', 'FontSize', 40); % Set y-axis label
zlabel('Final Site kWh Output', 'FontSize', 40); % Set z-axis label
% title('Impact of Height and Row Spacing on Site Electrical Output', 'FontSize', 40); % Set title
cb = colorbar; % Create a colorbar
cb.Label.String = 'Electrical Output (kWh)'; % Label the colorbar
cb.Label.FontSize = 40; % Set colorbar label font size
cb.FontSize = 40; % Set colorbar tick label font size
ax = gca; % Get the current axes
ax.FontSize = 40; % Set font size for the x, y, and z tick labels
view(2); % Set the view to 2D for better visibility of contours
hold off; % Release the hold

% Plotting the percentage change in agricultural yield with contours
figure;
surf(Xgrid, Ygrid, Z2grid); % Create the surface plot
shading interp; % Apply interpolated shading
hold on;
contour3(Xgrid, Ygrid, Z2grid, 20, 'k'); % Add contour lines
xlabel('Row Spacing (m)', 'FontSize', 40); % Set x-axis label
ylabel('Height', 'FontSize', 40); % Set y-axis label
zlabel('Percentage Change in Agricultural Yield', 'FontSize', 40); % Set z-axis label
% title('Impact of Height and Row Spacing on Agricultural Yield Percentage Change', 'FontSize', 40); % Set title
cb = colorbar; % Create a colorbar
cb.Label.String = 'Relative Agricultural Yield (%)'; % Label the colorbar
cb.Label.FontSize = 40; % Set colorbar label font size
cb.FontSize = 40; % Set colorbar tick label font size
ax = gca; % Get the current axes
ax.FontSize = 40; % Set font size for the x, y, and z tick labels
view(2); % Set the view to 2D
hold off; % Release the hold


% Save the results table to a file
% save('height_validation3.mat', 'results');
% 
% % Plotting the percentage change in electrical output with contours
% figure;
% surf(Xgrid, Ygrid, Z3grid);
% shading interp;
% hold on;
% contour3(Xgrid, Ygrid, Z3grid, 20, 'k');  % Add black contour lines
% xlabel('Row Spacing (m)', 'FontSize', 40, 'Interpreter', 'latex');
% ylabel('Height', 'FontSize', 40, 'Interpreter', 'latex');
% zlabel('Percentage Change in Electrical Output', 'FontSize', 40, 'Interpreter', 'latex');
% title('Impact of Height and Row Spacing on Electrical Yield Percentage Change', 'FontSize', 40, 'Interpreter', 'latex');
% colorbar;
% view(2);
% hold off;