%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Readme %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This file is used to explain the purpose of all functions in the Folder %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There six functions to generate every array structure % 
% Three for Billboard & Three for T-shaped %
% Two more functions to generate the underlaying 1D super nested array %
% The function URA_segment, supplemented with small_ura, is used to find the uDOF % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% When running the code "Derive_Lc_Lu_vs_N" % 
% The user could select the intended 2D array %
% Select the number of elements using either:
    % A predefined range of elements, 
    % User defined range, 
    % Random numbers 

% The software will simulate the selected array with the selected number of
% elements, simulate the uDOF, calculate the uDOF based on Table II % 
% Plot the simulated and the analytical expressions for uDOF %
% Plot the FODC 
% Draw a rectangle on the uDOF of the last N 

% The code also calculate the
    % Number of unique lags in the FODC
    % Aperture size as per Table III

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For 2D-DOA Estimation, you need just to
    % Extract the uDOF
    % Extract the corresponding measurements
    % Do spatial smoothing
    % Perform the estimation

