
% % Data for tomato and wheat growth
% tomato_values = [80.88, 66.3, 47.75, 74.96];
% wheat_values = [8.59, 7.6, 5.1, 6.31];
% categories = {'Reported Data', 'Open Field Modelling', 'Agrivoltaic Modelling', 'Irrigated Agrivoltaic Modelling'};
% 
% % Define the positions for the bars to ensure they are side by side
% tomato_positions = 1:4;
% wheat_positions = tomato_positions + 0.4;  % Reduced offset to make the bars touch
% 
% % Create the bar chart
% figure;
% yyaxis left; % Left Y-axis for tomato growth
% b1 = bar(tomato_positions, tomato_values, 0.4, 'FaceColor', 'b'); % Blue bars for tomato
% 
% % Customize the left Y-axis for tomato
% ylabel('Tomato Growth (t/ha)', 'FontSize', 28); % Y-axis label for tomato growth
% ax = gca;
% ax.YColor = 'b'; % Set the color of the y-axis to blue
% ax.XTick = (tomato_positions + wheat_positions) / 2; % Center x-tick marks between the two sets of bars
% ax.XTickLabel = categories; % Set custom category labels
% ax.XTickLabelRotation = 45; % Rotate labels for better readability
% ax.FontSize = 28; % Set font size for the axis tick labels and axes labels
% ylim([0 100]); % Adjust the y-axis limits based on tomato data
% 
% % Add a second bar graph on the same plot for wheat growth
% yyaxis right; % Right Y-axis for wheat growth
% b2 = bar(wheat_positions, wheat_values, 0.4, 'FaceColor', 'r'); % Red bars for wheat
% 
% % Customize the right Y-axis for wheat
% ylabel('Wheat Growth (t/ha)', 'FontSize', 28); % Y-axis label for wheat growth
% ylim([0 10]); % Adjust the y-axis limits based on wheat data
% ax.YColor = 'r'; % Set the color of the y-axis to red
% ax.FontSize = 28; % Ensure font size for the right axis is consistent
% 
% grid on; % Enable grid for easier reading of values
% 
% % Adjust legend to include both sets of bars
% legend([b1, b2], {'Tomato Growth', 'Wheat Growth'}, 'Location', 'northwest', 'FontSize', 28);
% 
% % Ensure the figure is maximized for better visibility
% set(gcf, 'WindowState', 'maximized');

% Data for tomato and wheat growth
tomato_values = [80.88, 66.3, 47.75, 74.96];
wheat_values = [8.59, 7.6, 5.1, 6.31];
categories = {'Reported Data', 'Open Field Modelling', 'Agrivoltaic Modelling', 'Irrigated Agrivoltaic Modelling'};

% Define the positions for the bars to ensure they are side by side
tomato_positions = 1:4;
wheat_positions = tomato_positions + 0.4;  % Reduced offset to make the bars touch

% Create the bar chart
figure;
yyaxis left; % Left Y-axis for tomato growth
b1 = bar(tomato_positions, tomato_values, 0.4, 'FaceColor', 'b'); % Blue bars for tomato

% Customize the left Y-axis for tomato
ylabel('Tomato Growth (t/ha)', 'FontSize', 28); % Y-axis label for tomato growth
ax = gca;
ax.YColor = 'b'; % Set the color of the y-axis to blue
ax.XTick = []; % Remove all x-tick marks
ax.FontSize = 28; % Set font size for the axis tick labels and axes labels
ylim([0 100]); % Adjust the y-axis limits based on tomato data

% Add a second bar graph on the same plot for wheat growth
yyaxis right; % Right Y-axis for wheat growth
b2 = bar(wheat_positions, wheat_values, 0.4, 'FaceColor', 'r'); % Red bars for wheat

% Customize the right Y-axis for wheat
ylabel('Wheat Growth (t/ha)', 'FontSize', 28); % Y-axis label for wheat growth
ylim([0 10]); % Adjust the y-axis limits based on wheat data
ax.YColor = 'r'; % Set the color of the y-axis to red
ax.FontSize = 28; % Ensure font size for the right axis is consistent

grid on; % Enable grid for easier reading of values

% Adjust legend to include both sets of bars
legend([b1, b2], {'Tomato Growth', 'Wheat Growth'}, 'Location', 'northwest', 'FontSize', 28);

% Ensure the figure is maximized for better visibility
set(gcf, 'WindowState', 'maximized');

