%% The task of the program is to read Z-stacks for each time point
%% and produce their maximum (and potentially average) projections
clear;
close all;
%--------------------------------------------------------------------------
InputDir = '_InputImages/';
OutputDir = '_InputImages/ZProj/';
%--------------------------------------------------------------------------
Files = dir([InputDir 'STACK_*.mat']);     % Takes only stacks
NbTimePoints = length(Files);

for i_Time = 1:NbTimePoints        % Loop on all stacks 
    % Read the stack
    FileName = [InputDir 'STACK_' int2str(i_Time) '.mat'];     
    Stack = load(FileName); 
    Stack = Stack.CurrentStack;
    NbSlices = length(Stack(1,1,:));    
    % Performing projections of the stack    
    MaxProj = Stack(:,:,1);      % Initialisation with the first z-slice
%     SumProj = Stack(:,:,1); 
    for i = 2:NbSlices 
        MaxProj = max(MaxProj, Stack(:,:,i));
%         SumProj = SumProj + Stack(:,:,i);
%         imshow(MaxProj, []);
    end  
%     AverProj = SumProj / NbSlices;
%     imshow(MaxProj, []);
%     figure, imshow(AverProj, []);    
    % Saving projection images    
    save([OutputDir 'MAX_' int2str(i_Time) '.mat'], 'MaxProj'); 
%     save([OutputDir 'AVG_' int2str(i_Time) '.mat'], 'AverProj'); 
end



