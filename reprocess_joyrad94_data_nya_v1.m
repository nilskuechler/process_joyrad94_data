clearvars;
clear mex;
close all; 
clc;

% create netcdf4 files from lv0 and/or lv1 files for joyrad94
% Author: Nils Küchler
% created: 8 February 2017
% modified: 21 May 2017, Nils Küchler 

% add functions in subfolders
addpath(genpath('/home/hatpro/scripts/joyrad94/data_processing/current_version/'))

% set processing variables
%   compact-flag: If compact_flag = 0: only general file is created
%                    compact_flag = 1: only compact file is generated
%                    compact_flag = 2: both files are gerenrated
compact_flag = 2;


% highest moment to be calculated
% ='off' moments taken from files if available, ='dealias' spectra are dealiased first, 
% then moments are calculated, ='spec' moments calculated from spectra
moments_cal = 'dealias';


for ii = 2015:2017
    for iii = 1:12
        for iv = 1:31
            
            if ii == 2017 && iii > 1 && iv > 6
                continue
            end
            
            if ii == 2016 && iii == 6 && iv < 12
                continue
            end
            
           
            path_lv0 = ['/data/obs/site/nya/joyrad94/l0/' num2str(ii)...
                '/' num2str(iii,'%02d') '/' num2str(iv,'%02d') '/'];
            
            path_lv1 = ['/data/obs/site/nya/joyrad94/l1/' num2str(ii)...
                '/' num2str(iii,'%02d') '/' num2str(iv,'%02d') '/'];
            
            files_lv0 = dir([path_lv1 '*nya_2*nc']);
            if isempty(files_lv0)
                continue
            end
            
            % create new subfolder if does not yet exist
            if ~exist(path_lv0,'dir')
                mkdir(path_lv0)
            end
                
            for v = 1:numel(files_lv0)                
               
                % backup file from
                infile = [path_lv0 files_lv0(v).name];
                % processed output file
                outfile = [path_lv1 files_lv0(v).name];
                
                % first copy infile to get backup
                if ~exist(infile,'file') % backup does not yet exist
                    copyfile(outfile, infile);    
                end
                                
                disp(infile);
                
                fid = fopen(infile, 'r', 'l');
                if fid == -1
                    disp(['could not open ' infile])
                    continue
                end
                filecode = int32(fread(fid,1,'int'));
                fclose(fid);
                
                process_joyrad94_data(infile,outfile,filecode,compact_flag,'moments',moments_cal); 
                     
            end % v = 1:numel(files)

                        
            clearvars -except ii iii iv  moments_cal compact_flag
        end % iv
    end % iii
end % ii 

% call plot routine for status data
% try
%     status_plot_nya_v2();
% catch
% end

% close matlab
% quit
