%% The task of the program is to show the movie sequence so that the user
%% can choose the time point when catastrophy takes place
function [i_Frame] = f_ShowMovie()
clear;     
close all;
%--------------------------------------------------------------------------
ImageFolder = '_InputImages\ZProj\MAX_';
%--------------------------------------------------------------------------
i_Frame = 1;
while true     % Loop until the user chooses the plane of the movie
    FilePath = strcat(ImageFolder, int2str(i_Frame), '.mat');  
%     InitImage = double(imread(FilePath)); 
    InitImage = load(FilePath);     
    InitImage = InitImage.MaxProj;           
    figure, 
    imshow(InitImage, []); 
    
    UserEntry = input('"+" to go forward; "-" to go backwards; \nanything for return, "0" for global exit\n', 's');
    switch UserEntry
        case '+'
            i_Frame = i_Frame + 1;
        case '-'
            i_Frame = i_Frame - 1;
            if i_Frame < 1
                i_Frame = 1;
            end
        otherwise
            break
    end
end