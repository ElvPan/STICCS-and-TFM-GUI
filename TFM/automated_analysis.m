%%%%%%%%%%% Matlab script written in Grenoble and changed by Lina in 2022

%Before starting:
%change number of cells (nb_cells) in automated analysis.m
%change parameter in line 75 in track_film_iteratif_2cells.m

%% initialize stock of data paths
clear all;
close all;
nb_cells = 1;
beads=cell(2,nb_cells);
initial=cell(2,nb_cells);
brightfield=cell(2,nb_cells);
mask=cell(2,nb_cells);
addpath(genpath('iteratif_2cells_film'));
addpath(genpath('visualisation'));

%% stock data paths
% variables for the type of cells

filespec=["*CROP*bk beads*.tif","*CROP*ak beads*.tif","*CROP*bk bf*.tif","*CROP*mask*.tif"];
  
path=uigetdir('C:', 'Select folder containing data to be treated');
    % create list with subfolders only
    d = dir(path);
    dfolders = d([d(:).isdir]);                                     % remove all files (isdir property is 0)
    dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));  % remove '.' and '..' 
    selpath = {dfolders(:).name};                                   % generate list with cell folders
    
        
for i=1:nb_cells

    % different paths to the respective images of one cell
    d1=dir(fullfile(path,selpath(i),filespec(1)));
    d2=dir(fullfile(path,selpath(i),filespec(2)));
    d3=dir(fullfile(path,selpath(i),filespec(3)));
    d4=dir(fullfile(path,selpath(i),filespec(4)));

    beads{1,i}=d1.name;
    beads_file=beads{1,i};
    beads{2,i}=fullfile(d1.folder,'/');
    beads_path=beads{2,i};

    initial{1,i}=d2.name;
    initial{2,i}=fullfile(d2.folder,'/');
   
    brightfield{1,i}=d3.name;
    brightfield{2,i}=fullfile(d3.folder,'/'); 
    
    mask{1,i}=d4.name;
    mask{2,i}=fullfile(d4.folder,'/'); 

end

clear selpath;
clear filespec;
clear dfolders;
clear d;
clear d1;
clear d2;
clear d3;
clear d4;

%% run analyses
for i=1:nb_cells
    close all;
    track_film_iteratif_2cells(beads{1,i},beads{2,i},initial{1,i},initial{2,i}, brightfield{1,i}, brightfield{2,i}, mask{2,i}, mask{1,i}, mask{2,i}, mask{1,i})
    make_movies(beads{2,i});
end

%%
close all 
clear all