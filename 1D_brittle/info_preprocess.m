%===Path for main.m=================================
cpath=[pwd,'/'];
% Name of the folder for mesh
namepre='I2c';%input('Give name of the folder for preprocessing :');
addpath([cpath,'Post_processing'])
addpath([cpath,'solver'])
addpath([cpath,'Pre_processing/',namepre,'/mesh']) 
addpath([cpath,'Pre_processing/',namepre,'/fea']) 

