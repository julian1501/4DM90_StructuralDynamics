% exercise 1.2
runs = 5;
nels = [10, 50, 100, 200]; % number of elements to test

evTimes = zeros(runs,size(nels,2));
eigenfrequencies = zeros(6,size(nels,2));

for i = 1:length(nels)
    for j = 1:runs
        % show information about current run
        fprintf("Run %d of %d for %d elements\n",j,runs,nels(i))
        % run the analysis and store the time taken
        [efs, evTimes(j,i)] = discreteBeamAnalysis(nels(i),false,"eig",1); % unitsize = 1 for meters, 1e3 for mm
    end
    eigenfrequencies(:,i) = efs;
end

format short g
% Display individual run times in a table
runTimesTable = array2table(evTimes, ...
    'VariableNames', strcat("Elements ", string(nels)), ...
    'RowNames', strcat("Run ", string(1:runs)));
disp(runTimesTable)

% Calculate average run times
avgTimes = round(mean(evTimes), 2);

% Prepare table data
format short e
T = table(nels(:), avgTimes(:), eigenfrequencies(1,:)', ...
    eigenfrequencies(2,:)', eigenfrequencies(3,:)', eigenfrequencies(4,:)', ...
    eigenfrequencies(5,:)', eigenfrequencies(6,:)', ...
    'VariableNames', {'NumElements', 'AvgRunTime', 'Freq1', 'Freq2', 'Freq3', 'Freq4', 'Freq5', 'Freq6'});

disp(T)