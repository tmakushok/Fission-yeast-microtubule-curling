%% The task of the program is to interact with the user: ask for the time
%% point when catartrophy happens and for the exact position of the bundle
%% tip at that moment. 
%% Then the function tracking the position of MT tip in 3D is called, results are saved 
clear;
close all;
%-------------------------------------------------------------------------- 
%!!!--!!! Basis for the name of Z-projection image
ProjImPathBase = '_InputImages/ZProj/MAX_';
%!!!--!!! Path base for outputs
Path_MTimage = '_Output/MT_';
Path_UserTimeAndPos = '_Output/User_CatastrTimeAndPos.mat';
Path_MTsEnds = '_Output/MTsEndPoints.mat';
%--------------------------------------------------------------------------
User_TimeAndPos = [];
MTsEnds = [];
Nb = 1;
% !!! User determines the time when MT catastrophy takes place 
% (last moment before the catastrophy) before executing this program
% UserInput = f_ShowMovie();
UserInput = input('Time point of MT catastrophy, please ("0" for exit):\n');

while UserInput ~= 0
    % Opening projection image corresponding to the time point specified by
    % the user
    ProjImage = load([ProjImPathBase int2str(UserInput) '.mat']); 
    ProjImage = ProjImage.MaxProj;      % Already of type 'double'   
    imshow(ProjImage, []); 
%% Asking the user to press on the end-point of the MT undergoing cat-phy    
    % !!! Manual zoom can be done here 
    % using dragging_an_area_with_mouse zooming tool        
    zoom    
    % !!! To finish the pause, press any key on the keyboard
    pause
    % 2 points can be entered: one specifying the end of MT, 
    % the other- direction where MT goes
    [X_Input, Y_Input] = ginput(2);      
    % Adding 0s if only one point was entered by the user
    if length(X_Input) == 1
        X_Input = [X_Input; 0];
        Y_Input = [Y_Input; 0];
    end  
    % Keeping the values the user gave as input 
    % (to be able then to go back from results to images)
    TotalInput = [UserInput, X_Input(1), Y_Input(1), X_Input(2), Y_Input(2)];
    User_TimeAndPos = [User_TimeAndPos; TotalInput];    
    % Executing the function that finds 3D coordinates of points at the end
    % part of the MT. 
    [MTpts, h_Res] = f_MTpart3D(TotalInput, ProjImage);     % Output is MT points and final result image
    % Store the results for first MT tip in MTsEnds(1, 1), in three columns: [x1, y1, z1; x2, y2, z2;...]. 
    % Next MT tip is in MTsEnds(2, 1)
    % The first point corresponds to the tip of MT, as determined by the user.          
    MTsEnds = [MTsEnds; MTpts];
    % Save figure where the position of detected MT points is shown
    saveas(h_Res, [Path_MTimage int2str(Nb) '.fig']);
    Nb = Nb + 1;        
    % !!! User determines the next time when MT catastrophy takes place 
    % (last moment before the catastrophy)
%     UserInput = f_ShowMovie();
    UserInput = input('Time point of catastrophy, please ("0" for exit):\n');
end
% Saving time points and locations of catastrophies analysed by the user with this program
save(Path_UserTimeAndPos, 'User_TimeAndPos');
% Saving final result: MT points
save(Path_MTsEnds, 'MTsEnds');
