function [ err ] = nntest( hyperParas, modelparas, test_x, test_y )
%NNTEST 
%   

[netState, ~] = nnfp(hyperParas, modelparas, test_x, test_y);

probs = netState.activity{hyperParas.numLayer};

% [~, predict] = max(probs);
% [~, target] = max(test_y);

% results_test = abs(probs-test_y);
% disp(probs);
size_i = size(test_y, 1);
size_j = size(test_y, 2);
cnt = 0;
for i=1:1:size_i
   for j=1:1:size_j
      
      if round(probs(i,j)) ~= test_y(i,j)
         cnt = cnt + 1; 
      end
   end
end

% bad = find(predict ~= target);
% err = numel(bad) / size(test_x, 2);
err = cnt / (size_i*size_j);

end