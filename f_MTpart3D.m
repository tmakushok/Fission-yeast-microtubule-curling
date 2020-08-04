%% The task of the function is to follow points in 3D corresponding to 
%% MT end (the tip of the MT is parameter of the function)
function [MTs, h_Res] = f_MTpart3D(TotalInput, AverImage)    
%--------------------------------------------------------------------------
%!!!--!!! The kernel size for median filtering of the initial stack
KernSize = 5;
%!!!--!!! Number of points detected on MT from its tip inwards (not counting the first one)
% Now: the number of points is defined by the distance between 2 points defined by the user 
% NbMT_Pts = 14;
%!!!--!!! Diameter of the circle that is used to find the next point along
%- MT bundle (diameter a bit bigger than MT bundle width)
CircleDiam = 5;     % 11 for RS microscope
%!!!--!!! Coefficient showing how much will be taken off from both sides of
% intensity profiles, to avoid coming back of the traced curve along the
% same MT
ProfCutCoeff = 1/4;         % 1/4;     % 1/6;
%!!!--!!! To find optimal Z plane (the one with biggest intensity at that pixel)
% only 'Z_RangePart' adjacent planes are considered, not to get to
% another bundle at completely different z-level in the cell
Z_RangePart = 3;        % 2 was before
%--------------------------------------------------------------------------
CellNb = 0;
AllMTs = [];
MTs = {};

Stack = load(['_InputImages/STACK_' int2str(TotalInput(1)) '.mat']);
Stack = Stack.CurrentStack;
NbSlices = length(Stack(1,1,:));
%% Visualise all slices of the stack and projection image
% figure(); 
% for i = 1:NbSlices   
%     figure();
%     imshow(Stack(:,:,i), []);
%     pause(0.1);
% end
% imshow(AverImage, []);
%% Filtering of all slices of the stack      
for i_Slice = 1:NbSlices
%     figure, imshow(Stack(:, :, i_Slice), []);
    Stack(:, :, i_Slice) = medfilt2(Stack(:, :, i_Slice), [KernSize KernSize]);
% 	  figure, imshow(Stack(:, :, i_Slice), []);  
end                
%% Projection of nearby z slices (3 slices for every z position)
LocalProj{1} = (Stack(:,:,1) + Stack(:,:,2)) / 2;
for i = 2:NbSlices - 1
    LocalProj{i} = (Stack(:,:,i - 1) + Stack(:,:,i) + Stack(:,:,i + 1)) / 3;
end
LocalProj{NbSlices} = (Stack(:,:,NbSlices - 1) + Stack(:,:,NbSlices)) / 2; 
% for i = 1:NbSlices   
%     figure();
%     imshow(LocalProj{i}, []);
%     pause(0.1);
% end
%% Getting the MT tip and the second point that defines the direction of the MT
MTTipsI = uint16(TotalInput(3));    % Corresponds to Y
MTTipsJ = uint16(TotalInput(2));    % Corresponds to X
MTSecondPt = uint16(TotalInput(4:5));
%% Filling 'MTs' array with MT start points
% Structure : MTs{ImageNumber / CellNumber, MT_Number}(i_InsideMT, :) = [i_Track, j_Track, Value];
MTs{1}(1, 1) = MTTipsI;
MTs{1}(1, 2) = MTTipsJ;
%% Going into 3D to find according Z coordinate of the tip
TipInZ = zeros(1, NbSlices);
for i_Slice = 1:NbSlices
    TipInZ(i_Slice) = Stack(MTTipsI, MTTipsJ, i_Slice);
end
% figure, grid on;
% line(1:NbSlices, TipInZ, 'Color', [.8 0 0], 'Marker', 'o'); 

[a, MTs{1, 1}(1, 3)] = max(TipInZ');           
%% Finding continuous MT points in 3D, from starting point defined by the user                
% NbPts = uint16(sqrt((TotalInput(2) - TotalInput(4)) .^ 2 + (TotalInput(3) - TotalInput(5)) .^ 2) / 2);    % Take every second pixel, approximately
LineX = [TotalInput(2:2:4)];
LineY = [TotalInput(3:2:5)];
% Positions on the line where circle centers will be, random image is chosen there
[LineX, LineY, a] = improfile(LocalProj{1}, LineX, LineY);   
% [LineX, LineY, a] = improfile(LocalProj{1}, LineX, LineY, NbPts);   
for i_Pt = 1:length(LineX)   % Loop on MT bundle points (number is defined)
%     close all;     % ifnot, images accumulate for each point
%% Define a circle of diameter a bit bigger than MT bundle width   
    I_LastPoint = length(MTs{1}(:, 1));     % The number of the last row in MTs{1, i_MT}
    % !!! Circle center is going to be taken at the line defined by the
    % user, not at the last point found.     
    x0 = LineX(i_Pt);    % Circle center coordinates
    y0 = LineY(i_Pt);
%     x0 = double(MTs{1}(I_LastPoint, 2));    % Circle center coordinates
%     y0 = double(MTs{1}(I_LastPoint, 1));
    CircleAngles = 0:5:360;         % Angles for points on the circle
    X_Prof = (x0 + CircleDiam * cosd(CircleAngles) / 2)';
    Y_Prof = (y0 + CircleDiam * sind(CircleAngles) / 2)'; 
%% Intens. profiles will be taken from an image that is projection of 3 nearby z slices
    ProjImage = LocalProj{MTs{1}(I_LastPoint, 3)};
    
%     figure; imshow(ProjImage, []);
%     hold on; plot(X_Prof, Y_Prof, '-bo', 'MarkerSize', 3);

    % 'bicubic' makes the profile more smooth
    Profile = improfile(ProjImage, X_Prof, Y_Prof, 'bicubic');
%                 figure; plot(1:length(Profile), Profile, '-o', 'MarkerSize', 3); grid on; 

    if I_LastPoint ~= 1         % Take off one third of the circle only for 2nd and on points
%% Take off from intens. profile part of the circle
        % Find the point on current circle closest to the previous circle
        % center. This point is the main direction to follow next. 
        % A sector of the cercle all around this point is conserved for the
        % analysis.
        Dist = sqrt(double((X_Prof - LineX(i_Pt - 1)).^2 + (Y_Prof - LineY(i_Pt - 1)).^2));
        [a, I_CirclePoint] = min(Dist);
%         [a, I_CirclePoint] = min(abs(Dist - CircleDiam / 2));
%         % 2) Find point on the circle correspoding to it
%         Dist = sqrt((X_Prof - double(MTs{1}(I_PrevPoint, 2))).^2 + (Y_Prof - double(MTs{1}(I_PrevPoint, 1))).^2);
%         [a, I_CirclePoint] = min(Dist);     

%                     figure; imshow(ProjImage, []);
%                     hold on; plot(X_Prof, Y_Prof, '-bo', 'MarkerSize', 3);    
%                     hold on; plot(X_Prof(I_CirclePoint), Y_Prof(I_CirclePoint), '-ro', 'MarkerSize', 5);

        % Circle profile is 'opened' the point closest to the old circle center
        Shift = - I_CirclePoint;         % + 1;
    %     Dist3 = circshift(Dist2, Shift);        
        X_Prof = circshift(X_Prof, Shift);
        Y_Prof = circshift(Y_Prof, Shift);
        Profile = circshift(Profile, Shift);
%         figure; plot(1:length(Profile), Profile, '-o', 'MarkerSize', 3);
%         grid on; 
        % Cut one 'ProfCutCoeff' part of intens. profile at each profile end
        % It will provide directionality to MT bundles detection
        L = length(Profile);
        Cut = uint16(L * ProfCutCoeff);
        Profile = Profile(Cut:L - Cut);
        X_Prof = X_Prof(Cut:L - Cut);
        Y_Prof = Y_Prof(Cut:L - Cut);
%                     figure; plot(1:length(Profile), Profile, '-o', 'MarkerSize', 3); grid on;                 
%% Position of next MT segment end : positions of central pics maxima
%         DiffProfile = [];
%         for i_DP = 1:length(Profile) - 1
%             DiffProfile(i_DP) = Profile(i_DP + 1) - Profile(i_DP); 
%         end
% %                 figure; plot(1:length(DiffProfile), DiffProfile, '-o',
%                 'MarkerSize', 3); grid on; 
        % Linear approx. 
        Pts = 4;
        Ps = [];                
        x = (1:Pts)';
        y = Profile(1:Pts);        
        p = polyfit(x, y, 1);        
        Ps(1) = p(1);    
%                     figure; plot(1:length(Profile), Profile, '-o', 'MarkerSize', 3); grid on; 
        for i = 2:length(Profile) - Pts                        
            x = (i:i + Pts - 1)';       
            y = Profile(x);        
            p = polyfit(x, y, 1);
            Ps(i) = p(1);  
%             hold on;
%             plot(x, p(1) * x + p(2), '-ro', 'MarkerSize', 5);
        end    
%                     figure; plot(3:length(Ps) + 2, Ps, '-o', 'MarkerSize', 3); grid on;                 
        % find roots; see if they are max
        Ind = [];
        for i_Ps = 2:length(Ps)
            if (Ps(i_Ps - 1) * Ps(i_Ps) < 0) && (Ps(i_Ps - 1) > Ps(i_Ps))                        
                Ind = [Ind, i_Ps];
            end
        end                
        % take first max from center
        [a, ind] = min(abs(Ind - length(Ps) / 2));
        MaxPos = Ind(ind) + 1;      % ??? not always ???
        % There should be a maximum in the chosen area of the profile
        % ifnot - just take the maximum of the profile, even if it is not the pic
        if isempty(MaxPos)
            [a I_MaxProf] = max(Profile);
            continue    
        end
        x_MP = max(1, MaxPos - 3):MaxPos + 3;
        [MaxProf, ind] = max(Profile(x_MP));
        I_MaxProf = x_MP(ind);  
    end     % From here on things are done for all points, including the first one
    % For first point position of max of the profile is position of next MT segment end
    if i_Pt == 1
        [a I_MaxProf] = max(Profile);
    end

    MTs{1}(I_LastPoint + 1, 1) = uint16(Y_Prof(I_MaxProf));     %  + CenterCorr   
    MTs{1}(I_LastPoint + 1, 2) = uint16(X_Prof(I_MaxProf));     %  + CenterCorr
%% Find optimal Z plane (the one with biggest intensity at that pixel)
% Only some adjacent planes are considered, not to get to
% another bundle at completely different z-level in the cell
    TipInZ = zeros(NbSlices, 1);
    OldZ = MTs{1}(I_LastPoint, 3);
    for i_Slice = max(1, OldZ - Z_RangePart):min(NbSlices, OldZ + Z_RangePart)
        TipInZ(i_Slice) = Stack(MTs{1}(I_LastPoint + 1, 1), MTs{1}(I_LastPoint + 1, 2), i_Slice);
    end
%     figure, grid on;
%     line(1:NbSlices, TipInZ, 'Color', [.8 0 0], 'Marker', 'o');   
    
    [a, Ind] = max(TipInZ);    
    MTs{1}(I_LastPoint + 1, 3) = Ind;           %Ind + OldZ - Z_RangePart - 1;
end
%% Visualise one MT end region found
h_Res = figure(); imshow(AverImage, []);
hold on; plot(MTs{1}(:, 2), MTs{1}(:, 1), '-ro', 'MarkerSize', 5);       
%% Reorganize 'MTs' for output so that first column is X, then Y, then Z
MTs{1} = [MTs{1}(:,2), MTs{1}];
MTs{1}(:,3) = [];

 