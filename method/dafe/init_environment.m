%% set machine and directories
OS = computer;
switch OS
    case 'GLNXA64' % 64-bit linux system
        [status,cmdout] = system('hostname');
        if status==0 && strcmp(cmdout(1:4), 'x270')
            % cliang's laptop
            machine_type = 'x270';
            curr_working_dir = '/home/cliang/work/code/test/reid/DAFE-git/';
            slash ='/';
            if ~exist('run_mode', 'var'), run_mode = 'serial'; end
            if ~exist('debug_flag', 'var'), debug_flag = true; end
        elseif status==0 && strcmp(cmdout(1:5), 'n0255')
            % whu-hpc
            machine_type = 'whuhpc';
            setenv('TZ','Asia/Shanghai');
            curr_working_dir = '/data/liangchao/DAFE/';
            slash ='/';
            if ~exist('run_mode', 'var'), run_mode = 'parallel'; end
            if ~exist('debug_flag', 'var'), debug_flag = false;  end
        else
            error('found unregistered machine');
        end
        
    case 'PCWIN64' % 64-bit windows system
        [status,cmdout] = dos('ECHO %COMPUTERNAME%');
        if status==0 && strcmp(cmdout(1:15), 'PC-20170623DBSD');
            % cliang's mmap-pc
            machine_type = 'mmap-pc';            
            curr_working_dir = 'D:\work\code\test\DAFE\';
            slash ='\';
            if ~exist('run_mode', 'var'), run_mode = 'parallel'; end
            if ~exist('debug_flag', 'var'), debug_flag = false;  end
        elseif status==0 && strcmp(cmdout(1:15), 'LAPTOP-V7EOT95U');
            % zye's x1
            machine_type = 'x1';            
            curr_working_dir = 'D:\work\code\test\DAFE\';
            slash ='\';
            if ~exist('run_mode', 'var'), run_mode = 'serial'; end
            if ~exist('debug_flag', 'var'), debug_flag = true;  end
        else
            error('found unregistered machine');
        end
        
    case 'MACI64'
        error('found unregistered machine');
end
cd(curr_working_dir);
addpath(genpath(fullfile(curr_working_dir)));



%% open parallel/serial
if strcmp(run_mode, 'parallel')==1 
    if ~isempty(gcp('nocreate'))>0
        delete(gcp('nocreate'))
    end

    poolobj = parpool(feature('NumCores')); %parpool;
    batch_size = poolobj.NumWorkers;   
else
    batch_size = 1;
end