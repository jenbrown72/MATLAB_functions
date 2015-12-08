function [DATA, loadedFiles] = JB_loadTxtFile(MouseID)

startDir = 'C:\Users\adesniklab\Documents\BehaviorRawData\MouseIDs\';

if nargin==1; % if a filename was inputted into function
    
    currfile = MouseID;
    
    tempdataDir = char(strcat(startDir, currfile, '\txtFiles'));
    cd(tempdataDir);
    
    display(' ')
    display(' ')
    disp(['________________Analysing file: ' currfile ' __________________________________'])
    display(' ')
    display(' ')
    
end

%Check to see if files have already been loaded
DATA.loadedFiles = [];

%Initialise variables
% number of angles used in behavior
angleNumber = {'Auto','S1','S2','S6','S12'};
stimAngles = [225, 241, 254, 263, 266, 268, 270, 272, 274, 277, 284, 299, 315];

for i = 1:length(angleNumber)
    DATA.(char((angleNumber{i}))) = [];
end


currDir = cd;
dataDir = strcat(currDir,'\txtFiles\data');
cd(dataDir);
olderDir = dir('DATA.mat');

if ~isempty(olderDir)
    load('DATA.mat')
    loaded = 1;
    
else
    loaded = 0;
end

cd(currDir);
D = dir('*.txt');

for i=1:size((D),1);
    
    filed = 0;
    %    for i=1;
    
    if (loaded==1)
        
        if(any(strcmp(D(i).name, DATA.loadedFiles))) % if already loaded this data, move on
            disp(['DATA file: ' D(i).name ' Already Loaded'])
            filed = 1;
            
            %   continue
            
        end
        
    end
    
    if(filed==0)
        
        tempDATA = csvread(D(i).name,1,0); % Dont read first line - sometimes contains a string 'A'
        tempDATA(length(tempDATA),:) = []; %Delete the last line incase not fully read
        
        disp(' ')
        disp(['Loading file: ' D(i).name])
        
        C = unique(tempDATA(:,11));
        stimNumber = size(find(ismember(stimAngles, C)),2);
        
        if tempDATA(1,8) ==1; %% auto trial
            
            Key = 'Auto';
            Index = strfind(angleNumber,Key);
            
            for f=1:length(Index);
                
                if Index{f}==1
                    currFileNo = length(DATA.Auto)+1;
                    
                    DATA.(angleNumber{f}){currFileNo} = D(i);
                    DATA.(angleNumber{f}){currFileNo}.rawData = tempDATA;
                    
                    DATA.loadedFiles{i} = D(i).name;
                    
                    disp([(angleNumber{f}) ': ' D(i).name ])
                    disp(' ')
                    filed = 1;
                    
                else
                    continue
                end
                
            end
            
        else
            
            for k=1:length(angleNumber)
                
                Key = 'S';
                Index = strfind(angleNumber{k},Key);
                
                if ~isempty(Index)
                    
                    Value = sscanf(angleNumber{k}(Index(1) + length(Key):end), '%g', 1);
                    
                    if (Value==stimNumber)
                        
                        disp([angleNumber{k} ': ' D(i).name ])
                        disp(' ')
                        
                        currFileNo = length(DATA.(angleNumber{k}))+1;
                        
                        DATA.(angleNumber{k}){currFileNo} = D(i);
                        DATA.(angleNumber{k}){currFileNo}.rawData = tempDATA;
                        
                        DATA.loadedFiles{i} = D(i).name;
                        filed = 1;
                        
                    end
                    
                else
                    continue
                end
                
            end
            
        end
        
    end
    
    if (filed==0)
        
        countNumber = max(tempDATA(:,6));
        
        if (countNumber<12)
            
                msgbox(strcat(D(i).name,' was not uploaded - count < 12'));
        disp(' ')
        disp(' ')
        disp(['WARNING..... NO FOLDER FOR: ' D(i).name ])
        disp(' ')
        disp(' ') 
        
        else
        msgbox(strcat(D(i).name,' was not uploaded - stimNumber = ' , num2str(stimNumber)));
        disp(' ')
        disp(' ')
        disp(['WARNING..... NO FOLDER FOR: ' D(i).name ])
        disp(' ')
        disp(' ')
        
        end
    end
end

disp(['Upload complete'])
curr = cd;
cd(strcat(curr, '\txtFiles\data'))
save('DATA.mat', 'DATA', '-v7.3')

end