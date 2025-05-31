% This function post process log.overInterDyMFoam files
% input:
%   - directory which contains all log.overInterDyMFoam files (it will combine all data files together if there are multiple log files)
% output:
%   - all the data the log file records

clear; clc;

% filename = 'C:\Users\s2251764\OneDrive - University of Edinburgh\Company Projects\Mocean\OpenFOAM\PostProcess\hr1a\log.overInterDyMFoam.1';

% caseNames = {'hr1a', 'hr1b', 'hr1c', 'hr1d2', 'hr2', 'hr2-3', 'hr3', 'hr3c', 'hr4', 'hr4c', 'hr5c', 'hr5d', 'hr6d'};
% caseNames = {'hr1d_u'};
% caseNames = {'hr1d6-5', 'hr1d6', 'hr1d5-5', 'hr1d5', 'hr1d4-5'};
% caseNames = {'hr1d6-5'};
% caseNames = {'lhr1d', 'lhr1d4-5', 'lhr1d5', 'lhr1d5-5', 'lhr1d6', 'lhr1d6-5'};

caseNames = {'wbst4'};
% caseNames = {'wbst3_nomooring', 'wbst4'};
% caseNames = {'rw9_nomooring'};
% caseNames = {'wbst2_ucap_6dof_yaw'};
% caseNames = {'wbst2_ucap_6dof_yaw'};
% caseNames = {'wbst3', 'wbst3_nomooring'};
% caseNames = {'nnlhr1d4-5', 'nnlhr1d5', 'nnlhr1d5-5', 'nnlhr1d6', 'nnlhr1d6-5', 'nnlhr1d7'};
% caseNames = {'nnhr1d7', 'nnhr1d6-5', 'nnhr1d6', 'nnhr1d5-5', 'nnhr1d5', 'nnhr1d4-5'};
% caseNames = {'nnlhr1d8', 'nnlhr1d10'};
% caseNames = {'hr6b', 'hr8b', 'hr10b'};
% caseNames = {'gapRerun_hr10b', 'gapRerun_hr8b', 'gapRerun_hr6b', 'gapRerun_hr1d', 'gapRerun_hr2', 'gapRerun_hr2-3', 'gapRerun_hr3c', 'gapRerun_hr4c', 'gapRerun_hr5d'};

%% if apply transform

ifTransform = true; %true;
% patternFile = 'log.overInterDyMFoam';
% patternFile = 'log.overWaveDyMFoam';
patternFile = 'log.overWaveDyMFoamUCap';

%% Heave transform

if ifTransform
    heaveOffset = 0.00456;
else
    heaveOffset = 0;
end
% heaveOffset = 0;

%% Reverse ROT MAT 

% this rotation matrix is reverse of Rxyz

% these angles are opposite signs of whats already applied in mesh

%%%%%% fwd


if ifTransform
    rollFwd = 0;
    % pitchFwd = 1.72; % 3.9; 
    % yawFwd = -4.22; % 0;
    
    pitchFwd = -1.72; % 3.9; 
    yawFwd = 4.22; % 4.8; % 4.22; % 0;
else
    rollFwd = 0;
    pitchFwd = 0; % 3.9; 
    yawFwd = 0; % 0;
end

% original orientation matrix
orientMatOrig = [1 0 0;
    0 1 0;
    0 0 1];

% rotation matrices
RxFwd = [1 0 0;
    0 cosd(rollFwd) -sind(rollFwd);
    0 sind(rollFwd) cosd(rollFwd)];

RyFwd = [cosd(pitchFwd) 0 sind(pitchFwd);
    0 1 0;
    -sind(pitchFwd) 0 cosd(pitchFwd)];

RzFwd = [cosd(yawFwd) -sind(yawFwd) 0;
    sind(yawFwd) cosd(yawFwd) 0;
    0 0 1];

% rotMatFwd = RxFwd * RyFwd * RzFwd;
rotMatFwd = inv(RxFwd * RyFwd * RzFwd);

% rotMatFwd = RzFwd * RyFwd * RxFwd;
% rotMatFwd = inv(RzFwd * RyFwd * RxFwd);


%%%%%%%% aft
if ifTransform
    rollAft = 0;
    % pitchAft = -2.18; % 3.9; 
    % yawAft = -4.22; % 0;
    
    pitchAft = 2.18; % 3.9; 
    yawAft = 4.22; % 4.8; % 4.22; % 0;
else
    rollAft = 0;
    pitchAft = 0; % 3.9;
    yawAft = 0; % 0;
end

% original orientation matrix
orientMatOrig = [1 0 0;
    0 1 0;
    0 0 1];

% rotation matrices
RxAft = [1 0 0;
    0 cosd(rollAft) -sind(rollAft);
    0 sind(rollAft) cosd(rollAft)];

RyAft = [cosd(pitchAft) 0 sind(pitchAft);
    0 1 0;
    -sind(pitchAft) 0 cosd(pitchAft)];

RzAft = [cosd(yawAft) -sind(yawAft) 0;
    sind(yawAft) cosd(yawAft) 0;
    0 0 1];

% rotMatAft = RxAft * RyAft * RzAft;
rotMatAft = inv(RxAft * RyAft * RzAft);

% rotMatAft = RzAft * RyAft * RxAft;
% rotMatAft = inv(RzAft * RyAft * RxAft);


%%
for iii = 1:length(caseNames)

    clearvars -except caseNames iii rotMatFwd rotMatAft heaveOffset patternFile;

    caseName = caseNames{iii};
    saveDir = 'C:\Users\s2251764\OneDrive - University of Edinburgh\Company Projects\Mocean\OpenFOAM\GapStudy\PostPropData\';
    
    fileDir = ['C:\Users\s2251764\OneDrive - University of Edinburgh\Company Projects\Mocean\OpenFOAM\GapStudy\PostProcess\', caseName, '\'];
    fileList = dir(fullfile(fileDir, '*.*'));
    allFileNames = {fileList.name};
    
    noLogFiles = 0;
    idxLogFiles = [];
    
    for k = 1 : length(allFileNames)
        % Get this filename.
	    thisFileName = fullfile(fileList(k).folder, allFileNames{k});
	    % See if it contains our required pattern.
        if ~contains(thisFileName, patternFile, 'IgnoreCase', true)
        % Skip this file because the filename does not contain the required pattern.
            continue;
        end
    
        noLogFiles = noLogFiles+1;
        idxLogFiles = [idxLogFiles, k];
    
    end
    
    nn = 1;
    
    for k = 1 : noLogFiles
	    % Get this filename.
	    thisFileName = fullfile(fileList(idxLogFiles(k)).folder, allFileNames{idxLogFiles(k)});
    
        % Read the file content
        fileContent = fileread(thisFileName);
        
        %% Time step
        searchStringTime = [newline, 'Time = '];
        startIndexTime = strfind(fileContent, searchStringTime);
        
        pattern =  '([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)'; 
        
        if ~isempty(startIndexTime) % if the string is found
            % Extract the substring starting from the index of the specified string
            % Find the end of the line after the specified string
            for ii = 1:length(startIndexTime) % record every other, because fwd and aft the same
                endIndex = strfind(fileContent(startIndexTime(ii):startIndexTime(ii)+200), newline); % find the newline in the next 200 (arbirary, enough to cover the next few lines) characters to save time 
                endIndex = endIndex(2) + startIndexTime(ii) - 2;  % Adjust the index to include the end of the line, as Time = starts with new line, so it's the ind 2 new line we looking for
                substring = fileContent(startIndexTime(ii) + length(searchStringTime):endIndex);
            
                % Extract the vector from the substring using regular expression
                matches = regexp(substring, pattern, 'match');
            
                % Convert the matches to a numeric vector
                timeStep(ii, :) = str2double(matches);
            end
        else
            disp('String not found in the file.');
        end
        
    %     timeStep = timeStep(2:end); % the 1st time step recorded haven't got all the values calculated yet
        
        % save in cell array for all log files
        timeStepAllLog{nn} = timeStep;
    
        % clear variable for next iteration
        timeStep = [];

        %% DeltaT
        searchStringDt = [newline, 'deltaT = '];
        startIndexDt = strfind(fileContent, searchStringDt);
        
        pattern =  '([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)'; 
        
        if ~isempty(startIndexDt) % if the string is found
            % Extract the substring starting from the index of the specified string
            % Find the end of the line after the specified string
            for ii = 1:length(startIndexDt) % record every other, because fwd and aft the same
                endIndex = strfind(fileContent(startIndexDt(ii):startIndexDt(ii)+200), newline); % find the newline in the next 200 (arbirary, enough to cover the next few lines) characters to save time 
                endIndex = endIndex(2) + startIndexDt(ii) - 2;  % Adjust the index to include the end of the line, as Time = starts with new line, so it's the ind 2 new line we looking for
                substring = fileContent(startIndexDt(ii) + length(searchStringDt):endIndex);
            
                % Extract the vector from the substring using regular expression
                matches = regexp(substring, pattern, 'match');
            
                % Convert the matches to a numeric vector
                Dt(ii, :) = str2double(matches);
            end
        else
            disp('String not found in the file.');
        end
        
    %     timeStep = timeStep(2:end); % the 1st time step recorded haven't got all the values calculated yet
        
        % save in cell array for all log files
        DtAllLog{nn} = Dt;
        if strcmp(caseNames{1}, 'bst1_2m4_co') % 2nd log file for bst1_2m4_co had fixed time steps of 0.01
            if k == 2
                DtAllLog{nn} = [0.01 .* ones(21, 1); Dt];
            end
        end

    
        % clear variable for next iteration
        Dt = [];
        
        %% CoR (for Heave Surge Sway)
        % Specify the string you want to search for
        searchString = 'Centre of rotation:';  % Replace 'specificString' with the actual string you're looking for
        
        % Find the index of the specified string
        startIndex = strfind(fileContent, searchString);
        
        pattern = '([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)';  % Regular expression for scientific notation
        
        if ~isempty(startIndex) % if the string is found
            % Extract the substring starting from the index of the specified string
            % Find the end of the line after the specified string
            for ii = 1:length(startIndex) 
                endIndex = strfind(fileContent(startIndex(ii):startIndex(ii)+200), newline); % find the newline in the next 200 (arbirary, enough to cover the next few lines) characters to save time 
                endIndex = endIndex(1) + startIndex(ii) - 2;  % Adjust the index to include the end of the line
                substring = fileContent(startIndex(ii) + length(searchString):endIndex);
            
                % Extract the vector from the substring using regular expression
                matches = regexp(substring, pattern, 'match');
            
                % Convert the matches to a numeric vector
                corVector(ii, :) = str2double(matches);
            end
        else
            disp('String not found in the file.');
        end
        
        % record every other, because fwd and aft the same
        corVector = corVector(1:2:end, :);
        
        % match with time step size
        corVector = corVector(1:length(timeStepAllLog{nn}), :);
        
        corVectorAllLog{nn} = corVector;
    
        corVector = [];
        
        %% Orientation
        
        % Specify the string you want to search for
        searchStringOrien = 'Orientation:';  % Replace 'specificString' with the actual string you're looking for
        
        % Find the index of the specified string
        startIndexOrien = strfind(fileContent, searchStringOrien);
        
        pattern = '([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)';  % Regular expression for scientific notation
        
        if ~isempty(startIndexOrien) % if the string is found
            % Extract the substring starting from the index of the specified string
            % Find the end of the line after the specified string
            for ii = 1:length(startIndexOrien) 
                endIndex = strfind(fileContent(startIndexOrien(ii):startIndexOrien(ii)+200), newline); % find the newline in the next 200 (arbirary, enough to cover the next few lines) characters to save time 
                endIndex = endIndex(1) + startIndexOrien(ii) - 2;  % Adjust the index to include the end of the line
                substring = fileContent(startIndexOrien(ii) + length(searchStringOrien):endIndex);
            
                % Extract the vector from the substring using regular expression
                matches = regexp(substring, pattern, 'match');
            
                % Convert the matches to a numeric vector
                orienVector(ii, :) = str2double(matches);
            end
        else
            disp('String not found in the file.');
        end
        
        % extract fwd and aft orientation
        orienVectorFwd = orienVector(1:2:end, :); 

        orienVectorFwd_3x3 = permute(reshape(orienVectorFwd', 3, 3, []), [2, 1, 3]); % Reshape to 3x3 for applying rotation matrix
        orienVectorFwd_R_3x3 = pagemtimes(orienVectorFwd_3x3, rotMatFwd);                 % Multiply each slice by R
        orienVectorFwd = reshape(permute(orienVectorFwd_R_3x3, [2, 1, 3]), 9, [])'; % Back to vector

        orienVectorAft = orienVector(2:2:end, :);

        orienVectorAft_3x3 = permute(reshape(orienVectorAft', 3, 3, []), [2, 1, 3]); % Reshape to 3x3 for applying rotation matrix
        orienVectorAft_R_3x3 = pagemtimes(orienVectorAft_3x3, rotMatAft);                 % Multiply each slice by R
        orienVectorAft = reshape(permute(orienVectorAft_R_3x3, [2, 1, 3]), 9, [])'; % Back to vector
    
        orienVector = [];
        
        % match with time step size
        orienVectorFwd = orienVectorFwd(1:length(timeStepAllLog{nn}), :);
        orienVectorAft = orienVectorAft(1:length(timeStepAllLog{nn}), :);
        
        orienVectorFwdAllLog{nn} = orienVectorFwd;
        orienVectorAftAllLog{nn} = orienVectorAft;
    
        orienVectorFwd = [];
        orienVectorAft = [];
        orienVectorFwd_3x3 = [];
        orienVectorAft_3x3 = [];
        orienVectorFwd_R_3x3 = [];
        orienVectorAft_R_3x3 = [];
        
        %% Linear Velocity
        
        % Specify the string you want to search for
        searchStringLinVel = 'Linear velocity:';  % Replace 'specificString' with the actual string you're looking for
        
        % Find the index of the specified string
        startIndexLinVel = strfind(fileContent, searchStringLinVel);
        
        pattern = '([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)';  % Regular expression for scientific notation
        
        if ~isempty(startIndexLinVel) % if the string is found
            % Extract the substring starting from the index of the specified string
            % Find the end of the line after the specified string
            for ii = 1:length(startIndexLinVel) 
                endIndex = strfind(fileContent(startIndexLinVel(ii):startIndexLinVel(ii)+80), newline); % find the newline in the next 200 (arbirary, enough to cover the next few lines) characters to save time 
                endIndex = endIndex(1) + startIndexLinVel(ii) - 2;  % Adjust the index to include the end of the line
                substring = fileContent(startIndexLinVel(ii) + length(searchStringLinVel):endIndex);
            
                % Extract the vector from the substring using regular expression
                matches = regexp(substring, pattern, 'match');
            
                % Convert the matches to a numeric vector
                linVelVector(ii, :) = str2double(matches);
            end
        else
            disp('String not found in the file.');
        end
        
        % extract fwd and aft orientation
        linVelVectorFwd = linVelVector(1:2:end, :); 
        linVelVectorAft = linVelVector(2:2:end, :);
        
        linVelVector = [];
    
        % match with time step size
        linVelVectorFwd = linVelVectorFwd(1:length(timeStepAllLog{nn}), :);
        linVelVectorAft = linVelVectorAft(1:length(timeStepAllLog{nn}), :);
        
        linVelVectorFwdAllLog{nn} = linVelVectorFwd;
        linVelVectorAftAllLog{nn} = linVelVectorAft;
    
        linVelVectorFwd = [];
        linVelVectorAft = [];
        
        %% Angular Velocity
        
        % Specify the string you want to search for
        searchStringAngVel = 'Angular velocity:';  % Replace 'specificString' with the actual string you're looking for
        
        % Find the index of the specified string
        startIndexAngVel = strfind(fileContent, searchStringAngVel);
        
        pattern = '([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)';  % Regular expression for scientific notation
        
        if ~isempty(startIndexAngVel) % if the string is found
            % Extract the substring starting from the index of the specified string
            % Find the end of the line after the specified string
            for ii = 1:length(startIndexAngVel)
                if ii ~= length(startIndexAngVel)
                    endIndex = strfind(fileContent(startIndexAngVel(ii):startIndexAngVel(ii)+80), newline); % find the newline in the next 200 (arbirary, enough to cover the next few lines) characters to save time 
                else
                    endIndex = strfind(fileContent(startIndexAngVel(ii):end), newline); % not enough char to set arbitrary number to find in next 80
                end
        
                endIndex = endIndex(1) + startIndexAngVel(ii) - 2;  % Adjust the index to include the end of the line
                substring = fileContent(startIndexAngVel(ii) + length(searchStringAngVel):endIndex);
            
                % Extract the vector from the substring using regular expression
                matches = regexp(substring, pattern, 'match');
            
                % Convert the matches to a numeric vector
                angVelVector(ii, :) = str2double(matches);
            end
        else
            disp('String not found in the file.');
        end
        
        % extract fwd and aft orientation
        angVelVectorFwd = angVelVector(1:2:end, :); 
        angVelVectorAft = angVelVector(2:2:end, :);
    
        angVelVector = [];
        
        % match with time step size
        angVelVectorFwd = angVelVectorFwd(1:length(timeStepAllLog{nn}), :);
        angVelVectorAft = angVelVectorAft(1:length(timeStepAllLog{nn}), :);
        
        angVelVectorFwdAllLog{nn} = angVelVectorFwd;
        angVelVectorAftAllLog{nn} = angVelVectorAft;
    
        angVelVectorFwd = [];
        angVelVectorAft = [];
        
        %% Pressure Forces and Moments
        
        % Specify the string you want to search for
        searchStringPressure = 'Pressure : ';  % pressure force
        
        % Find the index of the specified string
        startIndexPressure = strfind(fileContent, searchStringPressure);
        
        pattern = '([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)';  % Regular expression for scientific notation
        
        if ~isempty(startIndexPressure) % if the string is found
            % Extract the substring starting from the index of the specified string
            % Find the end of the line after the specified string
            for ii = 1:length(startIndexPressure)
                if ii ~= length(startIndexPressure)
                    endIndex = strfind(fileContent(startIndexPressure(ii):startIndexPressure(ii)+80), newline); % find the newline in the next 200 (arbirary, enough to cover the next few lines) characters to save time 
                else
                    endIndex = strfind(fileContent(startIndexPressure(ii):end), newline); % not enough char to set arbitrary number to find in next 80
                end
        
                endIndex = endIndex(1) + startIndexPressure(ii) - 2;  % Adjust the index to include the end of the line
                substring = fileContent(startIndexPressure(ii) + length(searchStringPressure):endIndex);
            
                % Extract the vector from the substring using regular expression
                matches = regexp(substring, pattern, 'match');
            
                % Convert the matches to a numeric vector
                pressureVector(ii, :) = str2double(matches);
            end
        else
            disp('String not found in the file.');
        end
        
        % extract fwd and aft orientation
        pressureForceVectorFwd = pressureVector(1:4:end, :); 
        pressureMomentVectorFwd = pressureVector(2:4:end, :);
        pressureForceVectorAft = pressureVector(3:4:end, :);
        pressureMomentVectorAft = pressureVector(4:4:end, :);
    
        pressureVector = [];
        
        pressureForceVectorFwdAllLog{nn} = pressureForceVectorFwd;
        pressureMomentVectorFwdAllLog{nn} = pressureMomentVectorFwd;
        pressureForceVectorAftAllLog{nn} = pressureForceVectorAft;
        pressureMomentVectorAftAllLog{nn} = pressureMomentVectorAft;
    
        pressureForceVectorFwd = [];
        pressureMomentVectorFwd = [];
        pressureForceVectorAft = [];
        pressureMomentVectorAft = [];
        
        %% Viscous forces and moments
        
        % Specify the string you want to search for
        searchStringViscous = 'Viscous  : ';  % viscous forces and moments
        
        % Find the index of the specified string
        startIndexViscous = strfind(fileContent, searchStringViscous);
        
        pattern = '([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)';  % Regular expression for scientific notation
        
        if ~isempty(startIndexViscous) % if the string is found
            % Extract the substring starting from the index of the specified string
            % Find the end of the line after the specified string
            for ii = 1:length(startIndexViscous)
                if ii ~= length(startIndexViscous)
                    endIndex = strfind(fileContent(startIndexViscous(ii):startIndexViscous(ii)+80), newline); % find the newline in the next 200 (arbirary, enough to cover the next few lines) characters to save time 
                else
                    endIndex = strfind(fileContent(startIndexViscous(ii):end), newline); % not enough char to set arbitrary number to find in next 80
                end
        
                endIndex = endIndex(1) + startIndexViscous(ii) - 2;  % Adjust the index to include the end of the line
                substring = fileContent(startIndexViscous(ii) + length(searchStringViscous):endIndex);
            
                % Extract the vector from the substring using regular expression
                matches = regexp(substring, pattern, 'match');
            
                % Convert the matches to a numeric vector
                viscousVector(ii, :) = str2double(matches);
            end
        else
            disp('String not found in the file.');
        end
        
        % extract fwd and aft orientation
        viscousForceVectorFwd = viscousVector(1:4:end, :); 
        viscousMomentVectorFwd = viscousVector(2:4:end, :);
        viscousForceVectorAft = viscousVector(3:4:end, :);
        viscousMomentVectorAft = viscousVector(4:4:end, :);
    
        viscousVector = [];
        
        viscousForceVectorFwdAllLog{nn} = viscousForceVectorFwd;
        viscousMomentVectorFwdAllLog{nn} = viscousMomentVectorFwd;
        viscousForceVectorAftAllLog{nn} = viscousForceVectorAft;
        viscousMomentVectorAftAllLog{nn} = viscousMomentVectorAft;
    
        viscousForceVectorFwd = [];
        viscousMomentVectorFwd = [];
        viscousForceVectorAft = [];
        viscousMomentVectorAft = [];
    
        nn = nn+1;
    end
    
    %% Concat all data
    cleanedTimeStep = [];
    cleanedDt = [];
    cleanedCorVector = [];
    cleanedOrienVectorFwd = [];
    cleanedOrienVectorAft = [];
    cleanedLinVelVectorFwd = [];
    cleanedLinVelVectorAft = [];
    cleanedAngVelVectorFwd = [];
    cleanedAngVelVectorAft = [];
    
    cleanedPressureForceVectorFwd = [];
    cleanedPressureMomentVectorFwd = [];
    cleanedPressureForceVectorAft = [];
    cleanedPressureMomentVectorAft = [];
    
    cleanedViscousForceVectorFwd = [];
    cleanedViscousMomentVectorFwd = [];
    cleanedViscousForceVectorAft = [];
    cleanedViscousMomentVectorAft = [];
    
    for k = 1:noLogFiles
    
        cleanedTimeStep = [cleanedTimeStep; timeStepAllLog{k}];
        cleanedDt = [cleanedDt; DtAllLog{k}];
        cleanedCorVector = [cleanedCorVector; corVectorAllLog{k}];
        cleanedOrienVectorFwd = [cleanedOrienVectorFwd; orienVectorFwdAllLog{k}];
        cleanedOrienVectorAft = [cleanedOrienVectorAft; orienVectorAftAllLog{k}];
        cleanedLinVelVectorFwd = [cleanedLinVelVectorFwd; linVelVectorFwdAllLog{k}];
        cleanedLinVelVectorAft = [cleanedLinVelVectorAft; linVelVectorAftAllLog{k}];
        cleanedAngVelVectorFwd = [cleanedAngVelVectorFwd; angVelVectorFwdAllLog{k}];
        cleanedAngVelVectorAft = [cleanedAngVelVectorAft; angVelVectorAftAllLog{k}];
        
        cleanedPressureForceVectorFwd = [cleanedPressureForceVectorFwd; pressureForceVectorFwdAllLog{k}];
        cleanedPressureMomentVectorFwd = [cleanedPressureMomentVectorFwd; pressureMomentVectorFwdAllLog{k}];
        cleanedPressureForceVectorAft = [cleanedPressureForceVectorAft; pressureForceVectorAftAllLog{k}];
        cleanedPressureMomentVectorAft = [cleanedPressureMomentVectorAft; pressureMomentVectorAftAllLog{k}];
        
        cleanedViscousForceVectorFwd = [cleanedViscousForceVectorFwd; viscousForceVectorFwdAllLog{k}];
        cleanedViscousMomentVectorFwd = [cleanedViscousMomentVectorFwd; viscousMomentVectorFwdAllLog{k}];
        cleanedViscousForceVectorAft = [cleanedViscousForceVectorAft; viscousForceVectorAftAllLog{k}];
        cleanedViscousMomentVectorAft = [cleanedViscousMomentVectorAft; viscousMomentVectorAftAllLog{k}];
    
        if k ~= noLogFiles % if not the last log file, check for next log start time
            nextStartTime = timeStepAllLog{k+1}(1); 
            [closestTimeStep, replace_id] = min(abs(cleanedTimeStep - nextStartTime));
    
            cleanedTimeStep(replace_id:end, :) = []; % get rid of repeated rows before adding new data
            cleanedDt(replace_id:end, :) = [];
            cleanedCorVector(replace_id:end, :) = [];
            cleanedOrienVectorFwd(replace_id:end, :) = [];
            cleanedOrienVectorAft(replace_id:end, :) = [];
            cleanedLinVelVectorFwd(replace_id:end, :) = [];
            cleanedLinVelVectorAft(replace_id:end, :) = [];
            cleanedAngVelVectorFwd(replace_id:end, :) = [];
            cleanedAngVelVectorAft(replace_id:end, :) = [];
    
            cleanedPressureForceVectorFwd(replace_id:end, :) = [];
            cleanedPressureMomentVectorFwd(replace_id:end, :) = [];
            cleanedPressureForceVectorAft(replace_id:end, :) = [];
            cleanedPressureMomentVectorAft(replace_id:end, :) = [];
    
            cleanedViscousForceVectorFwd(replace_id:end, :) = [];
            cleanedViscousMomentVectorFwd(replace_id:end, :) = [];
            cleanedViscousForceVectorAft(replace_id:end, :) = [];
            cleanedViscousMomentVectorAft(replace_id:end, :) = [];
    
        end
    
    end
    
    %% Save data as struct files for easier access/load later on
    
    % save data according to the time steps of pressure calculation (if last time step not completed, will be missing pressure data)
    idx2Save = size(cleanedPressureForceVectorFwd, 1);
    
    %WEC OpenFoam
    WECOF.TimeStep = cleanedTimeStep(1:idx2Save, :);
    WECOF.Dt = cleanedDt(1:idx2Save, :);
    WECOF.CoR = cleanedCorVector(1:idx2Save, :);
    
    % adjust for heave offset
    WECOF.CoR(:, 3) = WECOF.CoR(:, 3) + heaveOffset;
    
    % fwd
    WECOF.Fwd.Orientation = cleanedOrienVectorFwd(1:idx2Save, :);
    % WECOF.Fwd.RollAng = atand(WECOF.Fwd.Orientation(:, 8)./WECOF.Fwd.Orientation(:, 9));
    % WECOF.Fwd.PitchAng = -atand(WECOF.Fwd.Orientation(:, 7)./(sqrt(WECOF.Fwd.Orientation(:, 1).^2+WECOF.Fwd.Orientation(:, 4).^2))); % negative because y is going into page
    % WECOF.Fwd.YawAng = atand(WECOF.Fwd.Orientation(:, 4)./WECOF.Fwd.Orientation(:, 1));

    % roll pitch yaw RzRyRx order, putting negative signs in front of all 
    WECOF.Fwd.PitchAng = asind(WECOF.Fwd.Orientation(:, 7)); % negative because y is going into page
    WECOF.Fwd.RollAng = -asind(WECOF.Fwd.Orientation(:, 8)./cosd(WECOF.Fwd.PitchAng));
    WECOF.Fwd.YawAng = -asind(WECOF.Fwd.Orientation(:, 4)./cosd(WECOF.Fwd.PitchAng));

    % yaw pitch roll RxRyRz order, putting negative signs in front of all
    % WECOF.Fwd.PitchAng = -asind(WECOF.Fwd.Orientation(:, 3)); % negative because y is going into page
    % WECOF.Fwd.RollAng = -asind(WECOF.Fwd.Orientation(:, 6)./-cosd(WECOF.Fwd.PitchAng));
    % WECOF.Fwd.YawAng = -asind(WECOF.Fwd.Orientation(:, 2)./-cosd(WECOF.Fwd.PitchAng));

    WECOF.Fwd.LinVel = cleanedLinVelVectorFwd(1:idx2Save, :);
    WECOF.Fwd.AngVel = cleanedAngVelVectorFwd(1:idx2Save, :);
    WECOF.Fwd.Force.Pressure = cleanedPressureForceVectorFwd(1:idx2Save, :);
    WECOF.Fwd.Force.Viscous = cleanedViscousForceVectorFwd(1:idx2Save, :);
    WECOF.Fwd.Force.Total = WECOF.Fwd.Force.Pressure + WECOF.Fwd.Force.Viscous;
    WECOF.Fwd.Moment.Pressure = cleanedPressureMomentVectorFwd(1:idx2Save, :);
    WECOF.Fwd.Moment.Viscous = cleanedViscousMomentVectorFwd(1:idx2Save, :);
    WECOF.Fwd.Moment.Total = WECOF.Fwd.Moment.Pressure + WECOF.Fwd.Moment.Viscous;
    
    WECOF.Fwd.AngAccel = (WECOF.Fwd.AngVel(:, 2) - [0; WECOF.Fwd.AngVel(1:end-1, 2)])./WECOF.Dt;

    % aft
    WECOF.Aft.Orientation = cleanedOrienVectorAft(1:idx2Save, :);
    % WECOF.Aft.RollAng = atand(WECOF.Aft.Orientation(:, 8)./WECOF.Aft.Orientation(:, 9));
    % WECOF.Aft.PitchAng = -atand(WECOF.Aft.Orientation(:, 7)./(sqrt(WECOF.Aft.Orientation(:, 1).^2+WECOF.Aft.Orientation(:, 4).^2))); % negative because y is going into page
    % WECOF.Aft.YawAng = atand(WECOF.Aft.Orientation(:, 4)./WECOF.Aft.Orientation(:, 1));

    % roll pitch yaw RzRyRx order, putting negative signs in front of all
    WECOF.Aft.PitchAng = asind(WECOF.Aft.Orientation(:, 7)); % negative because y is going into page
    WECOF.Aft.RollAng = -asind(WECOF.Aft.Orientation(:, 8)./cosd(WECOF.Aft.PitchAng));
    WECOF.Aft.YawAng = -asind(WECOF.Aft.Orientation(:, 4)./cosd(WECOF.Aft.PitchAng));

    % yaw pitch roll RxRyRz order, putting negative signs in front of all
    % WECOF.Aft.PitchAng = -asind(WECOF.Aft.Orientation(:, 3)); % negative because y is going into page
    % WECOF.Aft.RollAng = -asind(WECOF.Aft.Orientation(:, 6)./-cosd(WECOF.Aft.PitchAng));
    % WECOF.Aft.YawAng = -asind(WECOF.Aft.Orientation(:, 2)./-cosd(WECOF.Aft.PitchAng));
    
    WECOF.Aft.LinVel = cleanedLinVelVectorAft(1:idx2Save, :);
    WECOF.Aft.AngVel = cleanedAngVelVectorAft(1:idx2Save, :);
    WECOF.Aft.Force.Pressure = cleanedPressureForceVectorAft(1:idx2Save, :);
    WECOF.Aft.Force.Viscous = cleanedViscousForceVectorAft(1:idx2Save, :);
    WECOF.Aft.Force.Total = WECOF.Aft.Force.Pressure + WECOF.Aft.Force.Viscous;
    WECOF.Aft.Moment.Pressure = cleanedPressureMomentVectorAft(1:idx2Save, :);
    WECOF.Aft.Moment.Viscous = cleanedViscousMomentVectorAft(1:idx2Save, :);
    WECOF.Aft.Moment.Total = WECOF.Aft.Moment.Pressure + WECOF.Aft.Moment.Viscous;

    WECOF.Aft.AngAccel = (WECOF.Aft.AngVel(:, 2) - [0; WECOF.Aft.AngVel(1:end-1, 2)])./WECOF.Dt;
    
%     WECOF.Hinge.Moment = WECOF.Fwd.Force.Total .*  - WECOF.Aft.Moment.Total;
    WECOF.Hinge.AngVelMinus = WECOF.Fwd.AngVel - WECOF.Aft.AngVel;
    WECOF.Hinge.AngVelPlus = WECOF.Fwd.AngVel + WECOF.Aft.AngVel;
    WECOF.Hinge.FlexAng = WECOF.Aft.PitchAng - WECOF.Fwd.PitchAng; % deg
    
    save([saveDir caseName], '-struct', "WECOF")

    disp(['Case ' caseName ' saving complete.'])
end