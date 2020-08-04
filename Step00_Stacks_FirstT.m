%% The task of the program is to create Z-stacks for each time point
%% for the case where files are organised such that first come 
%% all time-points for one z-positions for one 
clear;
close all;
%--------------------------------------------------------------------------
%!!!--!!! Number of Z-slices
SlicesTotal = 22;
Directory = '_InputImages/';
%--------------------------------------------------------------------------
Files = dir([Directory '*.tif']);     % Takes only images !!!
NbTimePoints = length(Files)/SlicesTotal;

% Reading the first slice to know the size of the images
FileName = [Directory Files(1).name];     
InitImage = double(imread(FileName)); 
s = size(InitImage);
CurrentStack = zeros(s(1), s(2), SlicesTotal);    % Preallocation

TimePoint = 1;
for i_Time = 1:NbTimePoints        % Loop on time points of the movie  
    for i_Zslice = 1:SlicesTotal
        FileNb = (i_Zslice - 1) * NbTimePoints + i_Time;
        FileName = [Directory Files(FileNb).name];     
        InitImage = double(imread(FileName)); 
        CurrentStack(:, :, i_Zslice) = InitImage;                      
    end            
%% Visualise all slices of the stack
%     figure();
%     for i = 1:SlicesTotal          
%         imshow(CurrentStack(:,:,i), []);
%     end    
    save([Directory 'STACK_' int2str(TimePoint) '.mat'], 'CurrentStack'); 
    
    CurrentStack = zeros(s(1), s(2), SlicesTotal);
    TimePoint = TimePoint + 1;
end



