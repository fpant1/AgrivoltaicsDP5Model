% Data from the table for output kWh
categories = {'Power prioritised', 'Agrivoltaic algorithm', 'Crop prioritised'};
kWhValues = [1977, 1917, 1750];
kWhValues= kWhValues*576;

2500*576

% Additional data for percentage agricultural growth
growthPercentages = [70.5, 72, 79];

% Create the bar chart
figure;
yyaxis left; % Left Y-axis for kWh output
b1 = bar(1:length(kWhValues), kWhValues, 0.4, 'FaceColor', 'flat'); % Use a narrower bar width
b1.FaceColor = 'b'; % Set color of kWh output bars

% Customize the left Y-axis
ylabel('Electrical Energy (kWh)', 'FontSize', 28); % Y-axis label for output kWh
ax = gca;
ax.YColor = 'b'; % Set the color of the y-axis to blue to match the kWh bars
ax.XTick = []; % Remove all x-tick marks
ax.FontSize = 28; % Set font size for the axis tick labels and axes labels
ylim([0 1500000]); % Set the y-axis limits for the left axis
% ax.YTick = 0:500:2500; % Set y-ticks in 250 increments

% Customize grid appearance
grid on; % Enable grid for easier reading of values
ax.GridLineStyle = '-'; % Solid grid lines
ax.GridColor = [0.8 0.8 0.8]; % Light grey grid color for better visibility
ax.GridAlpha = 0.5; % Transparency of grid lines

% Add a second bar graph on the same plot
yyaxis right; % Right Y-axis for percentage growth
b2 = bar(1:length(growthPercentages), growthPercentages, 0.4, 'FaceColor', 'flat'); % Use the same bar width
b2.FaceColor = 'r'; % Set color of growth percentage bars
% Offset the X data of the second bar group for alignment
b2.XData = b2.XData + 0.4;

% Customize the right Y-axis
ylabel('Agricultural Yield (%)', 'FontSize', 28); % Y-axis label for percentage growth
ylim([0 100]); % Set the y-axis limits for the right axis
ax = gca;
ax.YColor = 'r'; % Set the color of the y-axis to red to match the growth bars
ax.FontSize = 28; % Ensure font size for the right axis is consistent
ax.YTick = 0:20:100; % Set y-ticks in 20% increments

% Adjust legend and make sure it includes both sets of bars
legend([b1, b2], {'Electrical Energy ', 'Agricultural Yield'}, 'Location', 'northwest', 'FontSize', 28);

% Ensure the figure is maximized for better visibility
set(gcf, 'WindowState', 'maximized');

