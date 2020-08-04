%% The task of the program is to find cell cortex border and determine cell
%% parameters (cell ends, cell length, angle, profile of cell width)
% Result array contains: cell ends(1,2,3,4 (x1,y1,x2,y2)), 
% cell axis angle in degrees(5), cell width close to cell end(6), cell length(7)
function [CellParams] = f_CellParams(InitImage, PosInCell)
%--------------------------------------------------------------------------
%!!!--!!! Parameter for the artificial taking off of light halo around cells
AntiHalo = 5;
%!!!--!!! Minimal area of a cell, in pixels
MinArea = 1000;       
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
%--------------------------------------------------------------------------
FigNb = 1;
           
[m, n] = size(InitImage);        
figure, imshow(InitImage, []);
%% Filtering of the initial image to smooth the background random changes in intensity    
AverFiltered = medfilt2(InitImage, [3 3]);
% h = fspecial('average', 3);        %7);
% AverFiltered = imfilter(InitImage, h, 'replicate');      
% figure; imshow(AverFiltered, []); 
%% Determining the edges of the cells using 'edge'       
[CannyResult, ThresCanny] = edge(AverFiltered, 'canny');
% figure;  imshow(CannyResult);         
%% Dilate the image using structuring elements    
se = strel('disk', 3); 
BWsdil = imdilate(CannyResult, se);
%     se90 = strel('line', 7, 90);    % vertical structuring element 
%     se0 =  strel('line', 7, 0);      % horizontal structuring element                        
%     BWsdil = imdilate(CannyResult, [se90 se0]); 		%[se90 se0]); - for both vert. and horis. 
% figure, imshow(BWsdil);    
%% Fill Interior Gaps
BWdfill = imfill(BWsdil, 'holes');    
% figure, imshow(BWdfill), title('');                    
%% Filtering the Object (to make the segmented object look natural) 
h = fspecial('average', 5);
BWdfill = imfilter(BWdfill, h, 'replicate');   
%     BWdfill = medfilt2(BWdfill, [13 13]);    
% figure, imshow(BWdfill);
%% Erode
se = strel('disk', 3); 
BWdfill = imerode(BWdfill, se);
% figure, imshow(BWdfill); 
%% Filtering the Object (to make the segmented object look natural) 
% h = fspecial('average', 5);
% BWdfill = imfilter(BWdfill, h, 'replicate');   
BWdfill = medfilt2(BWdfill, [7 7]);    
% figure, imshow(BWdfill);
%% Getting some region statistics from images  
%-- 2) Converting of the binary image to a label matrix   
Labels = bwlabel(BWdfill);
Stats = regionprops(Labels, 'Area'); 
%% Taking off very small labelled regions
[Condition] = find([Stats.Area] > MinArea);
BigRegions = ismember(Labels, Condition);
Labels = bwlabel(BigRegions);
%imshow(BigRegions);                                            
% %% Erode the Object (to make the segmented object look natural)
%         seD = strel('diamond',1);
%         BWfinal = imerode(BigRegions,seD);
%         BWfinal = imerode(BWfinal,seD);
%         figure, imshow(BWfinal), title('After erosion');
%% Place an outline around segmented cells (for visual control over the segmentation process)
BWoutline = bwperim(BigRegions);   % To find perimeter pixels in binary image
Segout = InitImage; 
[OutlineGreyLevel, a] = max(max(InitImage));    
Segout(BWoutline) = OutlineGreyLevel * Coeff_OutlineGreyLevel;    
figure(FigNb), ResFigNb = FigNb; FigNb = FigNb + 1;
imshow(Segout, []);
% imagesc(Segout, [1, OutlineGreyLevel * Coeff_Imshow]); %title('outlined original image');
%% Preparing an image with only the concidered cell
Stats = regionprops(Labels, 'PixelList', 'BoundingBox');     % 'FilledImage', 
% Finding our cell amongst all cells (using 'PosInCell')
for i_Pix = 1:length(Stats)
    LinCl = find(Stats(i_Pix).PixelList(:, 1) == int16(PosInCell(3)) & Stats(i_Pix).PixelList(:, 2) == int16(PosInCell(4)));
    if ~isempty(LinCl)
        break
    end
end
% Extracting coordinates of the upper left corner and width of the
% bounding box that contained the image 
lc_x = Stats(i_Pix).BoundingBox(1) - 0.5; 
lc_y = Stats(i_Pix).BoundingBox(2) - 0.5; 
wid_x = Stats(i_Pix).BoundingBox(3); 
wid_y = Stats(i_Pix).BoundingBox(4); 
% Extracting the filled cell from BinaryImage inside the box 
% corresponding to the initial BoundingBox
AddCanvas = 0;
CroppedIm = BigRegions((lc_y + 1):(lc_y + wid_y), (lc_x + 1):(lc_x + wid_x));
figure, imshow(BigRegions, []);
figure, imshow(CroppedIm, []);
% Adding next 'layer' to the segmented image: next cell on
% the black BkGd added to the ones accumulated previously     
OneCellInNature = [zeros(lc_y, n); 
    zeros(wid_y, lc_x), CroppedIm, zeros(wid_y, n - lc_x - wid_x); 
    zeros(m - lc_y - wid_y, n)];  
% OneCellInNature = [zeros(lc_y - 1, n); 
%     zeros(wid_y, lc_x - 1), CroppedIm, zeros(wid_y, n - (lc_x - 1) - wid_x); 
%     zeros(m - (lc_y - 1) - wid_y, n)];      
figure, imshow(OneCellInNature, []);
%% Calculating coordinates of the tips and of the 'width tips' for the cell
% 'GoodCells' contains [Nb in 'Stats', cell length, cell width, cell
% angle, cell center (x,y)]
[CellEnds, CellWidthEnds, CellsPixels, GoodCells] = f_Cell4Tips_OneCurvedCell(OneCellInNature, ResFigNb, MinCellWidth, MaxCellWidth, MinCellLength, MaxCellLength);         
%% Visualisation of the cell center
%-- (for cells that were not taken off as being mitotic or PAA)
%TotalIntens = [];       % Accumulation of values for all cells on the image
figure(ResFigNb);            
line(GoodCells(1, 5), GoodCells(1, 6), 'Color', [.8 0 0], 'Marker', 'o'); 
%% Final output
CellParams = [CellEnds(1:4), GoodCells(4), GoodCells(3), GoodCells(2)];
