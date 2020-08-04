%% The task of the program is to show the movie sequence so that the user
%% can choose the time point when catastrophy takes place
function [CellParams] = f_CellParams(InitImage, PosInCell)
clear;     
close all;
%--------------------------------------------------------------------------
ImageFolder = '_InputImages\ZProj\MAX_';
%--------------------------------------------------------------------------
figure;
ImFiles = dir([ImageFolder, '*.mat']);   % Obtaining the list of files 
for i_Frame = length(ImFiles):-1:1     % Loop on the image files to analyse
    FilePath = strcat(ImageFolder, int2str(i_Frame), '.mat');  
    InitImage = load(FilePath);     
    InitImage = InitImage.MaxProj;           
    figure(i_Frame), 
    imshow(InitImage, []); 
end