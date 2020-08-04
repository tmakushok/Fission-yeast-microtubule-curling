%% The task of the program is 1) to smooth MTs final bits detected in previous
%% step; 2) find tangent vector at the MT end
clear;
close all;
%% Parameters
%!!!--!!! Paths
Path_MTsEnds = '_Output/MTsEndPoints.mat';
%!!!--!!! Length of 2D vector that helps visualise final vectors on proj. images
VectLen2D = 30;
%!!!--!!! Parameter for smoothing spline to smooth MT bundles
p_SmSpline = 0.01;
p_SmSplineXY = 0.01;
%!!!--!!! Number of points to do linear approximation in Oxy
Pts_LinApp = 5;
%!!!--!!! The length of the stretch of z-smoothing to use for finding
%-- z-direction
% StretchForZ = 0.5; % in x-y pixel units
%!!!--!!! Ratio between Z-slice-distance (here in microns) and XY-resolution (here in microns)
UmPerPx = 0.08;
ResolRatio = 0.2 / UmPerPx;  
ScaleBarLen = 5;    % Scale bar on control images will be this number of microns long
%!!!--!!! Paths
UserInputFile = '_Output/User_CatastrTimeAndPos.mat';
% ProjImPathBase = '_InputImages/ZProj/AVG_';   % Image files for finding cell parameters
BrightFieldImPath = '_InputImages/Tea1del36C#13_BrtField0002.tif';
PathBase2DVisual = '_Output/MT_';
%Result output
ResFile_PtCoord = '_Output/_ResultPtCoord.mat';
ResFile_Vector = '_Output/_ResultTVector.mat';
ResFile_CellParams = '_Output/CellParams.mat';
ResFile_ControlIm = '_Output/ControlImagesList.mat';
%% Finding cell parameters: cell ends(1,2,3,4 (x1,y1,x2,y2)), cell axis angle in degrees(5), 
%% cell width (better: cell width close to cell end) (6), cell length(7)
% In the beginning, not to do it for every point
BrFieldImBad = imread(BrightFieldImPath);    
% Inverting the image: it's bottom becomes it's top
s_BrF = size(BrFieldImBad);
BrFieldIm = zeros(s_BrF);
BrFieldIm(1:s_BrF(1), :) = BrFieldImBad(s_BrF(1):-1:1, :);
[AllCellsParams, CellsPixels, CellsOutlines] = f_CellParams_BrightField(BrFieldIm);  % for bright field images  
%% Initialisations
load(UserInputFile);    % 'User_TimeAndPos'
MTs = load(Path_MTsEnds);
MTs = MTs.MTsEnds;
% Initialisations
MTsSmoothed = MTs;     
MTendCoordCell = zeros(length(MTs), 6);
TanVect_CoordCell = zeros(length(MTs), 6);
AllControlImages = [];
for i_MTs = 1:length(MTs)
    close all;
    if isempty(MTs{i_MTs, 1})      % Badly detected MTs were put to []
        continue
    end
    TangentVector = zeros(1, 3);
    % Changing to 'double' data type
    MTs{i_MTs,1} = double(MTs{i_MTs,1});
    MTtheEnd = [MTs{i_MTs,1}(1, 1), MTs{i_MTs,1}(1, 2), MTs{i_MTs,1}(1, 3)];      % (x, y, z) of the very end of the MT
%% Changing Z coordinates from slice number units to 
%% same distance units as for XY (XY pixel size).    
    MTs{i_MTs,1}(:, 3) = (MTs{i_MTs,1}(:, 3) - 1) * ResolRatio;     % '- 1' is to have Z starting from 0 and not from 'ResolRatio'    
    MTsSmoothed{i_MTs, 1} = MTs{i_MTs,1};
%% Visualisation of one MT_end in 3 D
%     figure; plot3(MTs{i_MTs,1}(:, 2), MTs{i_MTs,1}(:, 1), MTs{i_MTs,1}(:, 3), '-o', 'MarkerSize', 3);
%     grid on; xlabel('X'); ylabel('Y'); zlabel('Z');  
%% Ordering MT_end points so that X or Y coordinate always grows,
%% not to have zig-zag lines when smoothing
    len = length(MTsSmoothed{i_MTs, 1}(:, 1));
    if range(MTsSmoothed{i_MTs, 1}(:, 1)) > range(MTsSmoothed{i_MTs, 1}(:, 2))
        if MTsSmoothed{i_MTs, 1}(1, 1) < MTsSmoothed{i_MTs, 1}(len, 1) 
            % If X grows, sort in ascending order
            MTsSmoothed{i_MTs, 1} = sortrows(MTsSmoothed{i_MTs, 1}, 1);
        else
            % Ifnot, sort in descending order 
            MTsSmoothed{i_MTs, 1} = sortrows(MTsSmoothed{i_MTs, 1}, -1);
        end
    else
        if MTsSmoothed{i_MTs, 1}(1, 2) < MTsSmoothed{i_MTs, 1}(len, 2) 
            % If Y grows, sort in ascending order
            MTsSmoothed{i_MTs, 1} = sortrows(MTsSmoothed{i_MTs, 1}, 2);
        else
            % Ifnot, sort in descending order 
            MTsSmoothed{i_MTs, 1} = sortrows(MTsSmoothed{i_MTs, 1}, -2);
        end                
    end
%% Smoothing (only z- coordinate) with smoothing spline           
    ToSmooth = double(MTsSmoothed{i_MTs, 1}(:, 3));                            
    x = 1:length(ToSmooth);  
    FitType = fittype('smoothingspline');
    FitOptions = fitoptions('Method', 'SmoothingSpline', 'SmoothingParam', p_SmSpline);
    cfunZ = fit(x', ToSmooth, FitType, FitOptions);
    
    figure();
    plot(x, ToSmooth, '-'); 
    hold on; grid on;
    x_MorePts = 1:0.1:length(ToSmooth);  
    plot(x_MorePts, cfunZ(x_MorePts), '-r');  
        
    MTsSmoothed{i_MTs, 1}(:, 3) = cfunZ(x);
%% Smoothing (in x-y plane) with smoothing spline      
    X_ToSmooth = double(MTsSmoothed{i_MTs, 1}(:, 1));                        
    Y_ToSmooth = double(MTsSmoothed{i_MTs, 1}(:, 2));
    % Choosing which of the coordinates has bigger range of values
    % The chosen one is going to be the basis for smoothing
    if range(X_ToSmooth) > range(Y_ToSmooth)
        FitType = fittype('smoothingspline');
        FitOptions = fitoptions('Method', 'SmoothingSpline', 'SmoothingParam', p_SmSplineXY);
        cfunXY = fit(X_ToSmooth, Y_ToSmooth, FitType, FitOptions);
        MTsSmoothed{i_MTs, 1}(:, 2) = cfunXY(MTsSmoothed{i_MTs, 1}(:, 1));   
    else 
        a = X_ToSmooth;
        X_ToSmooth = Y_ToSmooth;
        Y_ToSmooth = a;   
        
        FitType = fittype('smoothingspline');
        FitOptions = fitoptions('Method', 'SmoothingSpline', 'SmoothingParam', p_SmSplineXY);
        cfunXY = fit(X_ToSmooth, Y_ToSmooth, FitType, FitOptions);
        MTsSmoothed{i_MTs, 1}(:, 1) = cfunXY(MTsSmoothed{i_MTs, 1}(:, 2));           
    end
    
    figure();
    plot(X_ToSmooth, Y_ToSmooth, '-o'); 
    hold on;
    plot(X_ToSmooth, cfunXY(X_ToSmooth), '-r'); grid on  
%% To find the direction of the tangent vector in Oxy
%% we take the direction of a number of last points (linear approx.) 
%% taken from smoothed MT positions on Oxy
    Ind = Pts_LinApp:-1:1;
    x_LinApp = MTsSmoothed{i_MTs, 1}(Ind, 1);
    y_LinApp = MTsSmoothed{i_MTs, 1}(Ind, 2);
    p_XY = polyfit(x_LinApp, y_LinApp, 1);  
    LinApp = p_XY(1) * x_LinApp + p_XY(2);
    % Visualisation
    figure; plot(MTs{i_MTs, 1}(:, 1), MTs{i_MTs, 1}(:, 2), '-o', 'MarkerSize', 3);
    hold on; grid on;
    plot(x_LinApp, LinApp, '-ro', 'MarkerSize', 5);
%% To find the direction of the tangent vector in Z-direction
%% we take the direction of the last stretch of the Z-smoothing
    StretchForZ = 1;
    while true
       if (MTsSmoothed{i_MTs, 1}(1, 1) == MTsSmoothed{i_MTs, 1}(1 + StretchForZ, 1)) && (MTsSmoothed{i_MTs, 1}(1, 2) == MTsSmoothed{i_MTs, 1}(1 + StretchForZ, 2)) 
           StretchForZ = StretchForZ + 1;
       else
           break
       end
    end            
    DistXY = sqrt((MTsSmoothed{i_MTs, 1}(1, 1) - MTsSmoothed{i_MTs, 1}(1 + StretchForZ, 1)) .^ 2 + (MTsSmoothed{i_MTs, 1}(1, 2) - MTsSmoothed{i_MTs, 1}(1 + StretchForZ, 2)) .^ 2);
    p_Z = (cfunZ(1) - cfunZ(1 + StretchForZ)) / DistXY;       
%% Visualisation in 3D          
    % And the final vector is NOT ... (1, p_XY(1), p_Z)    
    if MTsSmoothed{i_MTs, 1}(1, 1) < MTsSmoothed{i_MTs, 1}(len, 1) 
        Incr = -6;
    else
        Incr = 6;
    end                  
    x1 = double(MTsSmoothed{i_MTs,1}(1, 1));
    x2 = double(x1 + Incr);
    
    y1 = p_XY(1) * (x1) + p_XY(2);
    y2 = p_XY(1) * (x2) + p_XY(2);
    
    Incr_Z = sqrt((x1 - x2) .^ 2 + (y1 - y2) .^ 2); % Distance along the MT, in XY- distance units
    z1 = MTsSmoothed{i_MTs,1}(1, 3);
    z2 = MTsSmoothed{i_MTs,1}(1, 3) + p_Z * Incr_Z;          % / ResolRatio;  
    % Visualisation
    figure; plot3(MTs{i_MTs,1}(:, 1), MTs{i_MTs,1}(:, 2), MTs{i_MTs,1}(:, 3), '-o', 'MarkerSize', 3);
    grid on; xlabel('X'); ylabel('Y'); zlabel('Z'); 
    hold on;
    plot3([x1 x2], [y1 y2], [z1 z2], '-ro', 'MarkerSize', 5);       
    rotate3d on
    % A breakpoint can be set here to be able to interact with the plot
    rotate3d off
%% Visualisation on 2D image
    % Creation of 2D vector of length 1
    x_TV = x2 - x1;
    y_TV = y2 - y1;    

    R_TV = sqrt(x_TV .^ 2 + y_TV .^ 2);
    % Normalizing vector length to VectLen2D
    % 'x_TV' is independent variable
    a1_TV = y_TV / x_TV;    % (y = a1 * x)    
    % We have to be careful with the sign of x_TV:
    if x_TV > 0
        x_TV = sqrt(1 / (1 + a1_TV .^ 2)); 
    else
        x_TV = - sqrt(1 / (1 + a1_TV .^ 2)); 
    end
    y_TV = a1_TV * x_TV;    
    % Making vector longer
    x_TV = x_TV * VectLen2D; 
    y_TV = y_TV * VectLen2D;
    
    h = figure(); 
    ContrFigNb = h + 1;     % '+1' as 'open' opens in nest image after 'figure'
    open([PathBase2DVisual int2str(i_MTs) '.fig']);
    hold on;
    plot([x1 - x_TV, x1], [y1 - y_TV, y1], '-go', 'MarkerSize', 3); 
    plot([x1, x1 + x_TV], [y1, y1 + y_TV], '-go', 'MarkerSize', 3); 
    plot([x1 + x_TV], [y1 + y_TV], 'go', 'MarkerSize', 6);     
    
    CurrDir = pwd;
    k = strfind(CurrDir, '\');
    if i_MTs < 10
        Nb = ['0' int2str(i_MTs)];
    else
        Nb = int2str(i_MTs);
    end
    Path2DControlImage = ['_Output\' CurrDir(k(length(k)) + 1:length(CurrDir)) '_MT' Nb '.fig'];    
    AllControlImages = [AllControlImages; Path2DControlImage];    
%% Creation of the tangent vector of length 1
    x_TV = x2 - x1;
    y_TV = y2 - y1;
    z_TV = z2 - z1;

    R_TV = sqrt(x_TV .^ 2 + y_TV .^ 2 + z_TV .^ 2);
    % Normalizing vector length to 1: "x_TV .^ 2 + y_TV .^ 2 + z_TV .^ 2 = 1"
    % 'x_TV' is independent variable
    a1_TV = y_TV / x_TV;    % (y = a1 * x)
    a2_TV = z_TV / x_TV;
    % We have to be careful with the sign of x_TV:
    if x_TV > 0
        x_TV = sqrt(1 / (1 + a1_TV .^ 2 + a2_TV .^ 2)); 
    else
        x_TV = - sqrt(1 / (1 + a1_TV .^ 2 + a2_TV .^ 2)); 
    end
    y_TV = a1_TV * x_TV;
    z_TV = a2_TV * x_TV;
    TangentVector = [x_TV, y_TV, z_TV];
%% Finding which cell is analysed right now
    PosInCell = User_TimeAndPos(i_MTs, :);
    TheCellNb = 0;
    for i_Pix = 1:length(CellsPixels)
        LinCl = find(CellsPixels{i_Pix}(:, 1) == int16(PosInCell(4)) & CellsPixels{i_Pix}(:, 2) == int16(PosInCell(5)));
        if ~isempty(LinCl)
            TheCellNb = i_Pix;
            break
        end
    end    
    if TheCellNb == 0
        for i_Pix = 1:length(CellsPixels)
            LinCl = find(CellsPixels{i_Pix}(:, 1) == int16(PosInCell(2)) & CellsPixels{i_Pix}(:, 2) == int16(PosInCell(3)));
            if ~isempty(LinCl)
                TheCellNb = i_Pix;
                break
            end
        end    
    end
    if TheCellNb == 0
        MTendCoordCell(i_MTs, :) = [0 0 0 0 0 0]
        continue
    end
    CellParams = AllCellsParams(TheCellNb, :);
%% Transforming the MT end coordinates and tangent vector into one cell's coordinates:
%% Find closest to the MT cell end (which will be origin of coordinates)
    D1 = sqrt((MTtheEnd(1) - CellParams(1)) .^ 2 + (MTtheEnd(2) - CellParams(2)) .^ 2);
    D2 = sqrt((MTtheEnd(1) - CellParams(3)) .^ 2 + (MTtheEnd(2) - CellParams(4)) .^ 2);
    IsFirstCellEnd = 0;
    if D1 < D2
        Origin = [CellParams(1), CellParams(2)];        % (x, y) of the new origin of coordinates
        if Origin(1) < CellParams(3)
            IsFirstCellEnd = 1;     % First cell end (the one from which Matlab takes angles) is the one with smaller x
        end
    else
        Origin = [CellParams(3), CellParams(4)];
        if Origin(1) < CellParams(1)
            IsFirstCellEnd = 1;     % First cell end (the one from which Matlab takes angles) is the one with smaller x
        end        
    end    
%% One_cell's coordinates for MT end
    % Putting the origin (0,0) into the cell end 
    % and inverting direction of Oy axis
    MTtheEnd_C = MTtheEnd;
    % !!!!! Axis Y is inverted in comparison with original Matlab's Y !!!!!
    MTtheEnd_C(1) = MTtheEnd(1) - Origin(1);   % For X 
    MTtheEnd_C(2) = -MTtheEnd(2) + Origin(2);   % For Y   
    % Rotating the coordinate system: going into polar coordinates
    R = sqrt(MTtheEnd_C(1) .^ 2 + MTtheEnd_C(2) .^ 2);      % Distance from (0,0) to the point
    % Find initial polar angle of the point
    % inverse tangent of dy / dx; result in degrees
    TempAngle = atand(MTtheEnd_C(2) / MTtheEnd_C(1));   % 'Temporary' angle, suits only for 1st and 4th quadrants
    if MTtheEnd_C(1) > 0        % Direction towards the point is in first or fourth quadrants
        PtAngle = TempAngle;       
    else                        % Second or third quadrant           
        PtAngle = TempAngle + 180;
    end
    % Angle in one_cell's axes depends on the angle of the cell and on
    % which cell end it is
    CellAngle = CellParams(5);
    
    if IsFirstCellEnd   % If this cell end is concidered by Matlab as being the first one
        AngleDiff = - CellAngle;     % CellAngle can be positive or negative 
    else                % If this cell end is concidered by Matlab as being the second one
        if CellAngle > 0        % Cell 'leans' towards right            
            AngleDiff = + (180 - CellAngle);    
        else                    % Cell 'leans' towards left                     
            AngleDiff = - (180 + CellAngle);     
        end            
    end    
    NewAngle = PtAngle + AngleDiff;
    % Back from polar to cartesian coordinates
    MTtheEnd_C(1) = R * cosd(NewAngle);
    MTtheEnd_C(2) = R * sind(NewAngle);
    
    % Save new coordinates
    MTendCoordCell(i_MTs, 1:3) = MTtheEnd_C; 
    % Save cell width
    MTendCoordCell(i_MTs, 4) = CellParams(6);
    % Save cell length
    MTendCoordCell(i_MTs, 5) = CellParams(7);
    % Save cell number (line number in 'AllCellsParams')
    MTendCoordCell(i_MTs, 6) = TheCellNb;
%% One_cell's coordinates for the tangent vector   
    TangentVector_C = TangentVector;        
    % !!!!! Axis Y is inverted in comparison with original Matlab's Y !!!!!
    TangentVector_C(2) = -TangentVector_C(2);   % For Y   
    % Rotating the coordinate system: going into polar coordinates
    R = sqrt(TangentVector_C(1) .^ 2 + TangentVector_C(2) .^ 2);      % Distance from (0,0) to the point
    % Find initial polar angle of the point
    % inverse tangent of dy / dx; result in degrees        
    TempAngle = atand(TangentVector_C(2) / TangentVector_C(1));     % 'Temporary' angle, suits only for 1st and 4th quadrants
    if TangentVector_C(1) > 0        % If vector is in first or fourth quadrants
        PtAngle = TempAngle;       
    else                        % Second or third quadrant           
        PtAngle = TempAngle + 180;
    end
    % Angle in one_cell's axes depends on the angle of the cell and on
    % which cell end it is
    CellAngle = CellParams(5);
    
    if IsFirstCellEnd   % If this cell end is concidered by Matlab as being the first one
        AngleDiff = - CellAngle;     % CellAngle can be positive or negative 
    else                % If this cell end is concidered by Matlab as being the second one
        if CellAngle > 0        % Cell 'leans' towards right            
            AngleDiff = + (180 - CellAngle);    
        else                    % Cell 'leans' towards left                     
            AngleDiff = - (180 + CellAngle);     
        end            
    end    
    NewAngle = PtAngle + AngleDiff;
    % Back from polar to cartesian coordinates
    TangentVector_C(1) = R * cosd(NewAngle);
    TangentVector_C(2) = R * sind(NewAngle);
    
    % Save new coordinates
    TanVect_CoordCell(i_MTs, 1:3) = TangentVector_C;   
    % Save cell width
    TanVect_CoordCell(i_MTs, 4) = CellParams(6); 
    % Save cell length
    TanVect_CoordCell(i_MTs, 5) = CellParams(7);
    % Save date of movie making, then '0', then number of the movie, then '0',
    % then number of MT
    k1 = strfind(CurrDir, '_Nb');
    k2 = strfind(CurrDir, '01_');     
    TanVect_CoordCell(i_MTs, 6) = str2num(CurrDir(k2 + 3:k2 + 4)) * 1000000 + str2num(CurrDir(k1 + 3:length(CurrDir))) * 1000 + i_MTs;    
    
    figure(ContrFigNb);  
    % Writing results on top of the control image
    ToWrite1 = ['Pos: ' num2str(MTtheEnd(1), 2) ' ' num2str(MTtheEnd(2), 2) ' ' num2str(MTtheEnd(3), 2) ...
        ' Vect: ' num2str(TangentVector(1), 2) ' ' num2str(TangentVector(2), 2) ' ' num2str(TangentVector(3), 2) ';'];
    ToWrite2 = ['Pos conv.: ' num2str(MTtheEnd_C(1), 2) ' ' num2str(MTtheEnd_C(2), 2) ' ' num2str(MTtheEnd_C(3), 2) ...
        ' Vect conv.: ' num2str(TangentVector_C(1), 2) ' ' num2str(TangentVector_C(2), 2) ' ' num2str(TangentVector_C(3), 2)];
    title([ToWrite1 ToWrite2]);
    text(20, 20, ToWrite1, 'Color', 'w'); 
    text(20, 50, ToWrite2, 'Color', 'w'); 
    % Putting scale bar on top of the control image 
    hold on    
    line([20, 20 + ScaleBarLen / UmPerPx], [s_BrF(1) - 20, s_BrF(1) - 20], 'LineWidth', 3, 'Color', 'w');
    % Overlay cells outline on top of the control image
    [Outlines_Y Outlines_X] = find(CellsOutlines);
    line(Outlines_X, Outlines_Y, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 2, 'Color', 'w');
    saveas(ContrFigNb, Path2DControlImage);
end
save(ResFile_PtCoord, 'MTendCoordCell');
save(ResFile_Vector, 'TanVect_CoordCell');
save(ResFile_ControlIm, 'AllControlImages');
save(ResFile_CellParams, 'AllCellsParams', 'CellsPixels');

















