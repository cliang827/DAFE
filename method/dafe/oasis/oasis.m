function model = oasis(data, class_labels, parms)
% model = oasis(data, class_labels, parms)
%
% Code version 1.3 May 2011 Fixed random seed setting
% Code version 1.2 May 2011 added call to oasis_m.m
% Code version 1.1 May 2011 handle gaps in class_labels
% 
%  Input:
%   -- data         - Nxd sparse matrix (each instance being a ROW)
%   -- class_labels - label of each data point  (Nx1 integer vector)
%   -- parms (do sym, do_psd, aggress etc.)
% 
%  Output: 
%   -- model.W - dxd matrix
%   -- model.loss_steps - a binary vector: was there an update at
%         each iterations
%   -- modeo.parms, the actual parameters used in the run (inc. defaults)
% 
%  Parameters:
%   -- aggress: The cutoff point on the size of the correction
%         (default 0.1) 
%   -- rseed: The random seed for data point selection 
%         (default 1)
%   -- do_sym: Whether to symmetrize the matrix every k steps
%         (default 0)
%   -- do_psd: Whether to PSD the matrix every k steps, including
%         symmetrizing them (defalut 0)
%   -- do_save: Whether to save the intermediate matrices. Note that
%         saving is before symmetrizing and/or PSD in case they exist
%         (default 0)
%   -- save_path: In case do_save==1 a filename is needed, the
%         format is save_path/part_k.mat
%   -- num_steps - Number of total steps the algorithm will
%         run (default 1M steps)
%   -- save_every: Number of steps between each save point
%         (default num_steps/10)
%   -- sym_every: An integer multiple of "save_every",
%         indicates the frequency of symmetrizing in case do_sym=1. The
%         end step will also be symmetrized. (default 1)
%   -- psd_every: An integer multiple of "save_every",
%         indicates the frequency of projecting on PSD cone in case
%         do_psd=1. The end step will also be PSD. (default 1)
%   -- use_matlab: Use oasis_m.m instead of oasis_c.c
%      This is provided in the case of compilation problems.
% 
%
% See Chechik et al. Large scale online learning of image
% similarity through ranking, J. Machine learning Research 2010.
% 
% (C) 2008-2010 Gal Chechik, Uri Shalit. 

  [N,dim] = size(data);
  W = eye(dim);
  [unused_values, inds] = sort(class_labels); %#ok
  class_labels = class_labels(inds);
  data = data(inds,:);
  classes = sort(unique(class_labels));
  num_classes = length(classes);
  
  % Translate class labels to serial numbers 1,2,...
  new_class_labels = zeros(size(class_labels));
  for i=1:length(classes)
    new_class_labels(class_labels==classes(i)) = i;
  end
  class_labels = new_class_labels;
  class_sizes = zeros(num_classes,1);
  class_start = zeros(num_classes,1);
  for k=1:num_classes
    class_sizes(k) = sum(class_labels==k);
    class_start(k) = find(class_labels==k, 1, 'first');
  end    
  
  % Initialize
  aggress = take_from_struct(parms, 'aggress', 0.1);
  rseed = take_from_struct(parms, 'rseed', 1);
  num_steps = take_from_struct(parms, 'num_steps', 10^6);
  do_sym = take_from_struct(parms, 'do_sym', 0);
  sym_every = take_from_struct(parms, 'sym_every', 1);
  do_psd = take_from_struct(parms, 'do_psd', 0);
  psd_every = take_from_struct(parms, 'psd_every', 1);  
  do_save = take_from_struct(parms, 'do_save', 1);
  save_every = take_from_struct(parms, 'save_every', ceil(num_steps/10));
  use_matlab = take_from_struct(parms, 'use_matlab', 0);
%   save_path = take_from_struct(parms, 'save_path', fullfile(pwd,'./temp/oasis_saves'));
%   mkdir(save_path);
  rand('seed', rseed); %#ok

  % Optimize 
  loss_steps = zeros(1, num_steps);
  num_batches = ceil(num_steps / save_every);
  steps_vec = ones(1,num_batches-1)*save_every;
  steps_vec(end+1) = num_steps - (num_batches-1)*save_every;
  
  fprintf('num_batches = %d, total num_steps = %d\n', num_batches, num_steps);

  data = full(data);
  for i_batch = 1:num_batches
    if use_matlab
      fprintf('run batch %d/%d, for %d steps ("." = 100 steps)\n', ...
              i_batch, num_batches, save_every);
      [W,loss_steps_batch] = oasis_m(W, data, class_labels, class_start, ...
                                     class_sizes, steps_vec(i_batch), ...
                                     aggress);
    else
      fprintf('run batch %d/%d, for %d steps ("." = 1K steps)\n', ...
              i_batch, num_batches, save_every);
      [W,loss_steps_batch] = oasis_c(W, data', class_labels, class_start, class_sizes, ...
                                     steps_vec(i_batch), aggress, rseed, i_batch);
    end
    fprintf('\n');
    loss_steps_batch(end:save_every) = 0;
    loss_steps((i_batch-1)*save_every+1:i_batch* save_every) = loss_steps_batch;
    if do_save
      temp_filename = sprintf('part_%03.0f.mat', i_batch);
      temp_fullname = fullfile(save_path, temp_filename);
      save(temp_fullname, 'W');
    end
    if do_sym
      if (mod(i_batch,sym_every) == 0) || i_batch == num_batches
        W = 0.5*(W+W');
      end
    end
    if do_psd
      if (mod(i_batch,psd_every) == 0) || i_batch == num_batches
        [V,D] = eig(0.5*(W+W'));
        W = V*max(D,0)*V';
        clear V D;
      end
    end
  end

  model.W = W; 
  model.loss_steps = loss_steps;
  model.num_loss_steps = num_batches*save_every;
  model.parms = parms;
  
end

% ====================================================================
function [val, out_options]=take_from_struct(options,fieldname,default)
%
% [VAL, OUT_OPTIONS] = TAKE_FROM_STRUCT(OPTIONS, FIELDNAME, DEFAULT)  
% 
% Take values from the options structure, use default if fiekd does
% not exist. Provide meaningful error messages. The function
% also updates the structure when the default is used.
% If default is not given, program aborts if field does not exist.
%
% Examples: 
% 
%  1. get the values of n_restarts, use default of 20 if field isnt set
%  [n_restarts, options] = take_from_struct(options,'n_restarts',20);
% 
%  2. get the values of n_restarts, abort if doesnt exist
%  [n_restarts, options] = take_from_struct(options,'n_restarts');
%
% (C) GAl Chechik, 2004 Software available for academic use. Other
% uses require explicit permission.
% 
  
  out_options = options;    

  try
    val = out_options.(fieldname);
  catch %#ok
    if(exist('default','var')==0)
      fprintf('\n\nError:\n');
      fprintf('    Field "%s" does not exist in structure\n\n', ...
	      fieldname);
      error('Trying to read from a field that does not exist');
    end    
    val=default;
    out_options.(fieldname) = val;
  end
  
  return
end