function  [W,loss_steps_batch] = oasis_m(W, data, class_labels, ...
                                         class_start, class_sizes, ...
                                         num_steps, aggress)
% Code version 1.3 May 2011: removed reset random seed  
% Code version 1.2 May 2011
%
% function  [W,loss_steps_batch] = oasis_m(W, data, class_labels, ...
%                                          class_start, class_sizes, ...
%                                          num_steps, aggress, rseed)
% 
% See Chechik et al. Large scale online learning of image
% similarity through ranking, J. Machine learning Research 2010. 
% (C) 2008-2010 Gal Chechik, Uri Shalit. 
%

  n_iters = num_steps;
  loss_steps_batch = zeros(num_steps,1);
  num_objects = length(class_labels);
  if size(data,1) ~= num_objects
    error('dimension mismatch');
  end
  if size(W,1) ~= size(data,2)
    error('dimension mismatch');
  end
  
  for i_iter = 1: n_iters
    if mod(i_iter, 100) == 0
      fprintf('.');
    end
    if mod(i_iter, 1000) == 0
      fprintf('%d', i_iter);
    end    
    if mod(i_iter, 10000) == 0
      fprintf('[%s]', datestr(now));
    end
    
    % Sample a query image
    p_ind   = ceil(rand(1)*num_objects);
    class = class_labels(p_ind);
    pos_ind = class_start(class) - 1 + ceil(rand(1)* class_sizes(class));
    
    neg_ind = ceil(rand(1)*num_objects);
    while  class_labels(neg_ind) == class
      neg_ind = ceil(rand(1)*num_objects);
    end
    p   = data(p_ind,:);
    samples_delta = [+1, -1] * data([pos_ind, neg_ind],:);
    
    loss = 1 - p*W*samples_delta';
    
    if loss > 0 % Update W

        grad_W = p'*samples_delta;
        loss_steps_batch(i_iter) = 1;

        norm_grad_W = sum(p .* p) * sum(samples_delta .*samples_delta);
        
        
        % constraint on the maximal update step size
        tau_val = loss/norm_grad_W; %loss / (V*V');
        tau = min(aggress, tau_val); 
        
        W = W + tau*grad_W;
        
     end
  end

  loss_steps_batch = sparse(loss_steps_batch);
