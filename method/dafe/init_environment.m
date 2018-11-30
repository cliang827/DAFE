%% set machine and directories
OS = computer;
switch OS
    case 'GLNXA64' % 64-bit linux system
        [status,cmdout] = system('hostname');
        if status==0 && strfind(cmdout, 'x270')
            % cliang's laptop
            machine_type = 'x270';
            cd('/home/cliang/work/code/test/reid/DAFE-git');
            slash ='/';
            if ~exist('run_mode', 'var'), run_mode = 'serial'; end
            if ~exist('debug_flag', 'var'), debug_flag = true; end
        elseif status==0 && strfind(cmdout, 'n0255')
            % whu-hpc
            machine_type = 'whuhpc';
            setenv('TZ','Asia/Shanghai');
            cd('/data/liangchao/DAFE');
            slash ='/';
            if ~exist('run_mode', 'var'), run_mode = 'parallel'; end
            if ~exist('debug_flag', 'var'), debug_flag = false;  end
        else
            error('found unregistered machine');
        end
        
    case 'PCWIN64' % 64-bit windows system
        [status,cmdout] = dos('ECHO %COMPUTERNAME%');
        if s==0 && strfind(cmdout, 'PC-20170623DBSD');
            % cliang's mmap-pc
            machine_type = 'mmap-pc';            
            cd('D:\work\code\test\DAFE');
            slash ='\';
            if ~exist('run_mode', 'var'), run_mode = 'parallel'; end
            if ~exist('debug_flag', 'var'), debug_flag = false;  end
        else
            error('found unregistered machine');
        end
        
    case 'MACI64'
        error('found unregistered machine');
end

%% open parallel/serial
if strcmp(run_mode, 'parallelization')==1 
    if ~isempty(gcp('nocreate'))>0
        delete(gcp('nocreate'))
    end

    poolobj = parpool(feature('NumCores')); %parpool;
    batch_size = poolobj.NumWorkers;   
else
    batch_size = 1;
end