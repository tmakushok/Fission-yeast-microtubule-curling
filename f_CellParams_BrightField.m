%% The task of the program is to find cell cortex border and determine cell
%% parameters (cell ends, cell length, angle, profile of cell width)
% Result array contains: cell ends(1,2,3,4 (x1,y1,x2,y2)), 
% cell axis angle in degrees(5), cell width close to cell end(6), cell length(7)
function [CellParams, CellsPixels, BWoutline] = f_CellParams_BrightField(InitImage, PosInCell)
%--------------------------------------------------------------------------
% %!!!--!!! Parameter for the artificial taking off of light halo around cells
% AntiHalo = 5;
%!!!--!!! Parameter used for filling in the space between the edges of the
%-- cells determined automaticly (in pixels)
MaxEdgeDistance = 100;
%!!!--!!! Maximal distance between centers of two detected areas when they
%-- are still considered as the same cell detected twice
DistDeDoubling = 45;
%!!!--!!! Max of badly segmented cells on an image (for 'ginput' function)
MaxBadCells = 30;
%!!!--!!! Minimal area of a cell, in pixels
MinArea = 1000;       
MinAreaCont = 100;       % Min area of a contour to be kept
%!!!--!!! Defining the size limits (in pixels) for objects to be recognized as being
%-- S.pombe cells
MinCellWidth = 30;  %15; for bin 2
MaxCellWidth = 100; %50;% for bin 2
MinCellLength = 90; %45;% for bin 2
MaxCellLength = 320;    %140;% for bin 2
%!!!--!!! Bigger the coefficient, brighter cells outline
Coeff_OutlineGreyLevel = 2;
%!!!--!!! Bigger the coefficient, less bright the cells 
Coeff_Imshow = 1;
%!!!--!!! Figure number
AllBinCellsFigNb = 100;
%!!!--!!! To add this number of points on top of the image to make
%-- correspondence with fluo image
AddTop = 3;
AddLeft = 5;
%!!!--!!! Take off this number of pixels from detected cell width
CellWidthCorr = 6;
%--------------------------------------------------------------------------
% Initialisations
% For checking detected cell border against fluo image
load('_InputImages/ZProj/MAX_1.mat');
Nb_Cell = 1;    
CellEnds = [];
CellWidthEnds = [];
CellsPixels = cell(0,0);
GoodCells = []; 
FigNb = 1;
[m, n] = size(InitImage);        
figure, imshow(InitImage, []);
FullImWithOneCell = zeros(m, n);        % To have an image of the full size, but with only one filled cell on it	
%% Add some lines to the image to make better correspondence to fluo image
InitImage = [InitImage(1:AddTop, :); InitImage(1:(m - AddTop), :)];
InitImage = [InitImage(:, 1:AddLeft), InitImage(:, 1:(n - AddLeft))];
%% Filtering of the initial image to smooth the background random changes in intensity    
AverFiltered = medfilt2(InitImage, [7 7]);

h = fspecial('average', 5);        %7);
AverFiltered = imfilter(AverFiltered, h, 'replicate');      
figure; imshow(AverFiltered, []);
% -- For Gain 100 on ERS:
% AverFiltered = medfilt2(InitImage, [3 3]);
% h = fspecial('average', 5);        %7);
% AverFiltered = imfilter(AverFiltered, h, 'replicate');      
% figure; imshow(AverFiltered, []);
%% ----- Determining the edges of the cells using 'edge'       
[CannyResult, ThresCanny] = edge(AverFiltered, 'canny');
figure;  imshow(CannyResult);         
%% Taking off very small labelled regions
Labels = bwlabel(CannyResult);
Stats = regionprops(Labels, 'Area'); 

[Condition] = find([Stats.Area] > MinAreaCont);
BigContours = ismember(Labels, Condition);
figure; imshow(BigContours, []); 

Labels = bwlabel(BigContours);
StatsCont = regionprops(Labels, 'Area', 'FilledImage', 'BoundingBox', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');

figure(AllBinCellsFigNb);
GoodContoursImage = zeros(m, n);
for i = 1:length(StatsCont)            
    % loop: looking at all edges longer than a min and not close to image border              
    if (StatsCont(i).Area > MinAreaCont) & (StatsCont(i).BoundingBox(1) > 2) & (StatsCont(i).BoundingBox(2) > 2)
        % (StatsCont(i).BoundingBox(1) ~= ???) & (StatsCont(i).BoundingBox(1) ~= ???)                        
        %close all;            
        [m_OneCont, n_OneCont] = size(StatsCont(i).FilledImage);
        BinaryImage = StatsCont(i).FilledImage;            
%             figure(), imshow(BinaryImage);            
        % Putting some 0s all around the edge (ifnot, edge analysed is too close to the border of the image)
        AddCanvas = 2.5*MaxEdgeDistance;
        BinaryImage = [zeros(AddCanvas, n_OneCont + 2*AddCanvas); 
            zeros(m_OneCont, AddCanvas), BinaryImage, zeros(m_OneCont, AddCanvas); 
            zeros(AddCanvas, n_OneCont + 2*AddCanvas)];            
%             figure(), imshow(BinaryImage);             
%% Extracting coordinates of the upper left corner and width of the
% bounding box that contained the image with one single edge
        lc_x = StatsCont(i).BoundingBox(1) - 0.5; 
        lc_y = StatsCont(i).BoundingBox(2) - 0.5; 
        wid_x = StatsCont(i).BoundingBox(3); 
        wid_y = StatsCont(i).BoundingBox(4); 
        % Extracting the filled cell from BinaryImage inside the box 
        % corresponding to the initial BoundingBox
        CroppedIm = BinaryImage((AddCanvas + 1):(AddCanvas + wid_y), (AddCanvas + 1):(AddCanvas + wid_x));
%% Adding next 'layer' to the segmented image: next cell on
% the black BkGd added to the ones accumulated previously            
        OneCellInNature = [zeros(lc_y - 1, n); 
            zeros(wid_y, lc_x - 1), CroppedIm, zeros(wid_y, n - (lc_x - 1) - wid_x); 
            zeros(m - (lc_y - 1) - wid_y, n)];            
        FullImWithOneCell = FullImWithOneCell + OneCellInNature;            
%% Dilate and smooth the image to get to the true cell border (as seen on fluorescence images)
%             se =  strel('disk', 2, 0);
%             OneCellInNature = imdilate(OneCellInNature, se);    % [se0 se90]);            
%             h = fspecial('average', 5);
%             CroppedIm = imfilter(CroppedIm, h, 'replicate');                                                  
        figure(AllBinCellsFigNb); 
        imshow(FullImWithOneCell, []); 
%% Measuring parameters of the filled cells (before it was done for the cell contour only)   
        Labels = bwlabel(OneCellInNature);
        Stats = regionprops(Labels, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'PixelList');           
%% Discarding cells that have an area that is too small
        if Stats.Area < MinArea
            continue
        end    
%% Finding the four cell tips                        
        [OneCell_CellEnds, OneCell_CellWidthEnds, OneCell_CellsPixels, OneCell_GoodCells] = f_Cell4TipsDetect(OneCellInNature, AllBinCellsFigNb, MinCellWidth, MaxCellWidth, MinCellLength, MaxCellLength);             
        if isempty(OneCell_GoodCells)
            continue
        end
% Accumulation of data for all the cells in the image
        CellEnds = [CellEnds; OneCell_CellEnds];
        CellWidthEnds = [CellWidthEnds; OneCell_CellWidthEnds];
        CellsPixels = [CellsPixels; OneCell_CellsPixels];
        GoodCells = [GoodCells; Nb_Cell, OneCell_GoodCells(2:6)];

        CellCenter = OneCell_GoodCells(5:6);
        Nb_Cell = Nb_Cell + 1;  
%% Visualisation of all contours found so far on the fluorescence image            
%             if isempty(CellsPixels)
%                 continue
%             end                        
%             for i_OnePix = 1:length(CellsPixels{1})     % Filling with 1s the places where 'final good' cells are detected   
%                 GoodContoursImage(CellsPixels{1}(i_OnePix, 2), CellsPixels{1}(i_OnePix, 1)) = 1;
%             end        
%             BWoutline = bwperim(GoodContoursImage);   % To find perimeter pixels in binary image            
%             Segout = FluoImage; 
%             [OutlineGreyLevel, a, a] = simple_max2D(FluoImage);
%             Segout(BWoutline) = OutlineGreyLevel;     
%             figure, imshow(Segout, []);   

%             AllCellEnds{i_Frame, Nb_Cell} = CellEnds;              
%             AllWidthEnds{i_Frame, Nb_Cell} = CellWidthEnds;                                                 
%% Solving the problem of having two detected areas overlapping 
%% (steming out from inner and outer edges of a cell) with slightly different widths 
%             AddedFlag = 0;            
        [LenWE, a] = size(CellWidthEnds);  
        if LenWE < 2
            continue
        end
%             switch LenWE    
%                 case 0
%                     continue            
%                 case 1                                 
%                     Result = [Result; GoodCells(2), GoodCells(3), GoodCells(5), GoodCells(6)];
% %                     AddedToRes = AddedToRes + 1;
%                     %CellsPixels = {CellsPixels; cell(1, 1)};
%                     %i_CP = length(CellsPixels);
%                     CellsPixels = {Stats.PixelList};                     
%                     continue
%             end
        for i_UnDoubl = (LenWE - 1):-1:1  % back count as more chances to have an overlap with one of the 'closely previous' areas             
            % Calculation of the distance between current cell and the
            % cell number i_DeDoubling
            CentersDist = sqrt((CellWidthEnds(i_UnDoubl, 5) - CellCenter(1))^2 + (CellWidthEnds(i_UnDoubl, 6) - CellCenter(2))^2);                
            if CentersDist <= double(DistDeDoubling)    % if centers of the two cells are very close
                CellWidthOld = GoodCells(i_UnDoubl, 3);     %sqrt((CellWidthEnds(i_UnDoubl, 3) - CellWidthEnds(i_UnDoubl, 1))^2 + (CellWidthEnds(i_UnDoubl, 4) - CellWidthEnds(i_UnDoubl, 2))^2);
                CellWidth = OneCell_GoodCells(3);
                if CellWidth < CellWidthOld                        
%                         i_ToRemove = find(Result((OldLenRes + 1):ResLen, 3) == CellWidthEnds(i_UnDoubl, 5) & Result((OldLenRes + 1):ResLen, 4) == CellWidthEnds(i_UnDoubl, 6));       
%                         i_ToRemove = i_ToRemove_CP + OldLenRes;
%                         Result(i_ToRemove, :) = [];   

                    % Remove those lines also in geometry matrixes                        
                    CellWidthEnds(i_UnDoubl, :) = [];
                    CellEnds(i_UnDoubl, :) = [];
                    GoodCells(i_UnDoubl, :) = [];
                    % Taking away from the pixels set
                    CellsPixels(i_UnDoubl) = [];
                    % Adding of the good information to 'Result'
%                         Result = [Result; CellLength, CellWidth, CellCenter];                         
                    % Adding to the pixels set
%                         CellsPixels = [CellsPixels; {Stats.PixelList}];                                                
%                         AddedFlag = 1;
                    break
                else                        
%                         AddedFlag = 1;   
                    % Taking away the last cell's information from
                    % geometry arrays
                    [LinG, a] = size(CellEnds);
                    CellWidthEnds(LinG, :) = [];
                    CellEnds(LinG, :) = []; 
                    GoodCells(LinG, :) = []; 
                    CellsPixels(LinG) = [];
                    break
                end
            end
        end
%             if AddedFlag == 0   % There was no cells with cell centers so close that it could have been an overlap
% %                 Result = [Result; CellLength, CellWidth, CellCenter];                                
% %                 CellsPixels = [CellsPixels; {Stats.PixelList}];		
% %                 AddedToRes = AddedToRes + 1;
%             end             
    end             
end     
%% Visualisation of cell objects after the un-overlapping procedure
FinalImage = zeros(m, n);
% Filling with 1s the places where 'final good' cells are detected
for i_vis = 1:length(CellsPixels)
    for i_OnePix = 1:length(CellsPixels{i_vis})     % Better to use 'size', but here it doesn't matter
%             i_OnePix = 1:length(CellsPixels{i_vis});
        FinalImage(CellsPixels{i_vis}(i_OnePix, 2), CellsPixels{i_vis}(i_OnePix, 1)) = 1;
    end
end
%     figure, imshow(FinalImage, []);
% Representing borders of the segmented image overlaid with initial image 
BWoutline = bwperim(FinalImage);   % To find perimeter pixels in binary image
Segout = MaxProj;     %InitImage;     
Segout(BWoutline) = max(max(Segout));     
figure, imshow(Segout, [430 550]);    
% Putting geometrical parameters on this overlaid image
Len_Res = length(CellsPixels);     
for i_vis = 1:Len_Res    
    % Visualisation of the cell length 
    line([CellEnds(i_vis, 1), CellEnds(i_vis, 3)], [CellEnds(i_vis, 2), CellEnds(i_vis, 4)]);            
    % Visualisation of the cell width            
    line([CellWidthEnds(i_vis, 1), CellWidthEnds(i_vis, 3)], [CellWidthEnds(i_vis, 2), CellWidthEnds(i_vis, 4)]);                    
    % Visualisation of the cell center        
%         line(GoodCells(i_vis, 5), GoodCells(i_vis, 6), 'Color', [.8 0 0], 'Marker', 'o');  
end   
%% Taking off the cells that are not well segmented (with mouse clicks)
[CellOff_x, CellOff_y] = ginput(MaxBadCells);  % User mouse-clicks-in coordinates, 'Enter' to continue    
% Find actually existing cells
%     GoodCellsLabels = bwlabel(FinalImage);         
for i_Off = 1:length(CellOff_x)     % Loop on the cells clicked        
    i_Pix = 1;
    while i_Pix < (length(CellsPixels) + 1)            
        % Finding lines in CellsPixels{i_Pix} in which we have x =
        % clicked_x and y = clicked_y              
        LinCl = find(int16(CellsPixels{i_Pix}(:, 1)) == int16(CellOff_x(i_Off)) & int16(CellsPixels{i_Pix}(:, 2)) == int16(CellOff_y(i_Off)));                                    
        if ~isempty(LinCl) % if it is the cell number i_Pix that contain clicked coordinates   
%                 GoodCellsLabels(find(GoodCellsLabels(:,:) == i_Pix)) = 0;               
            CellEnds(i_Pix, :) = [];  
            CellWidthEnds(i_Pix, :) = [];
            GoodCells(i_Pix, :) = [];
            CellsPixels(i_Pix) = [];                
            break                          
        end
        i_Pix = i_Pix + 1;
    end
end  
%     figure, imshow(BigRegions, []);    
%     GoodCellsLabels(find(GoodCellsLabels)) = 1;     
%     figure, imshow(GoodCellsLabels, []);       
%% Visualisation of the 'good cells' centers     
%figure(ResFigNb);
[Lin_CPix, Col_CPix] = size(CellEnds);
for i_GC = 1:Lin_CPix     % Loop on the good cells        
%         Fig = figure(); 
%         imshow(InitImage, []);
    line(CellEnds(i_GC, 5), CellEnds(i_GC, 6), 'Color', [.8 0 0], 'Marker', 'o');  
%         saveas(Fig, ['Output/CellNumber' num2str(i_GC)]);
end  
%% Final output
CellParams = [CellEnds(:, 1:4), GoodCells(:, 4), GoodCells(:, 3) - CellWidthCorr, GoodCells(:, 2)];



