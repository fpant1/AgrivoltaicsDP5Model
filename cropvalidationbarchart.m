% Define colors for the bars
colors = [0.8500, 0.3250, 0.0980;  % A shade of red
          0.9290, 0.6940, 0.1250;  % A shade of yellow
          0, 0.4470, 0.7410;       % A shade of blue
          0.4940, 0.1840, 0.5560]; % A shade of purple

% Data for tomato growth in Spain
tomato_values = [80.88, 66.3, 47.75, 74.96];
tomato_categories = {'Reported Data', 'Open Field Modelling', 'Agrivoltaic Modelling', 'Irrigated Agrivoltaic Modelling'};

% Data for wheat growth in the UK
wheat_values = [8.59, 7.6, 5.1, 6.31];
wheat_categories = {'Reported Data', 'Open Field Modelling', 'Agrivoltaic Modelling', 'Irrigated Agrivoltaic Modelling'};

% Create bar chart for tomato growth
figure;
bT = bar(tomato_values, 0.5, 'FaceColor', 'flat'); % Set width to 0.5 for thinner bars
set(gca, 'xticklabel', tomato_categories, 'FontSize', 30);
set(gca, 'ylim', [0 90], 'ytick', 0:10:90); % Set y-axis limits and ticks for tomato
for k = 1:length(tomato_values)
    bT.CData(k,:) = colors(k,:);
end
% title('Tomato Growth in Spain', 'FontSize', 30);
ylabel('Growth (t/ha)', 'FontSize', 30); % Units of tonnes per hectare
grid on;

% Create bar chart for wheat growth
figure;
bW = bar(wheat_values, 0.5, 'FaceColor', 'flat'); % Set width to 0.5 for thinner bars
set(gca, 'xticklabel', wheat_categories, 'FontSize', 30);
set(gca, 'ylim', [0 9], 'ytick', 0:1:9); % Set y-axis limits and ticks for wheat
for k = 1:length(wheat_values)
    bW.CData(k,:) = colors(k,:);
end
% title('Wheat Growth in Britain', 'FontSize', 30);
ylabel('Growth (t/ha)', 'FontSize', 30); % Units of tonnes per hectare
grid on;
