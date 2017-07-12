clearvars;
% clear mex;
close all; 
clc;

% create netcdf4 files from lv0 and/or lv1 files for joyrad94
% Author: Nils Küchler
% created: 8 February 2017
% modified: 21 May 2017, Nils Küchler 

% add functions in subfolders
addpath(genpath('/home/hatpro/scripts/joyrad94/data_processing/current_version/'))

% check last N days
N = 3;

% set processing variables
%   compact-flag: If compact_flag = 0: only general file is created
%                    compact_flag = 1: only compact file is generated
%                    compact_flag = 2: both files are gerenrated
compact_flag = 2;


% highest moment to be calculated
% ='off' moments taken from files if available, ='dealias' spectra are dealiased first, 
% then moments are calculated, ='spec' moments calculated from spectra
moments_cal = 'dealias';


for g = 1:N
    
    time_today = datevec( today-(N-g) );
    path_lv0 = ['/data/obs/site/nya/joyrad94/l1/' num2str(time_today(1))...
        '/' num2str(time_today(2),'%02d') '/' num2str(time_today(3),'%02d') '/*lv0'];
    
    files_lv0 = dir(path_lv0);
    if isempty(files_lv0)
        continue
    end
    

    for h = 1:numel(files_lv0)
 
        
        % ######################## start with level 0 (lv0) files
        infile = [path_lv0(1:end-4) files_lv0(h).name];
        disp(infile);
        
        outfile = [infile(1:end-3) 'nc'];
        % outfile = strrep(outfile,'joyce','nya');
        outfile = strrep(outfile,'joyrad94_','joyrad94_nya_');
        
        if exist(outfile,'file') == 2
            disp([outfile 'nc already exists']); continue
        end

        
        fid = fopen(infile, 'r', 'l');
        if fid == -1
            disp(['could not open ' infile])
            continue
        end
        filecode = int32(fread(fid,1,'int'));
        fclose(fid);       

        process_joyrad94_data(infile,outfile,filecode,compact_flag,'moments',moments_cal);
        
        if exist(outfile,'file') == 2
            delete(infile)
        end

         outfile = strrep(outfile,'_20','_compact_20');
         infile = strrep(infile,'lv0','lv1');
        if exist(outfile,'file') == 2
            delete(infile)
        end     
        

    
    end % h = 1:numel(files)

    status_plot_nya_v2(time_today);
    
clearvars -except g N moments_cal compact_flag
end % g = 1:2

% call plot routine for status data
% try
%     status_plot_nya_v2();
% catch
% end

% close matlab
quit
