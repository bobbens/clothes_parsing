% run test on the full dataset using the bigram model
function [cmx, score] = run_test( ...
   train, ...
   test, ...
   ftype, ...
   num_nn, ...
   item_limit)

% compute total number of classes
n_test    = length(test.filenames);
n_classes = n_test/test.n_test_per_class;

% initialize results
cmx       = zeros([n_classes n_classes]);
n_correct = 0;
n_total   = 0;
score     = 0;

% classify test images
for n = 1:n_test
   % load exemplar
   [features, x_pos, y_pos] = load_exemplar(test.filenames{n}, ftype);
   % classify exemplar 
   id_vec = query_db(features, x_pos, y_pos, num_nn, item_limit);
   id = id_vec(1);
   % update confusion matrix
   cmx(test.class_id(n), id) = cmx(test.class_id(n), id) + 1;
   % update score
   if (test.class_id(n) == id)
      n_correct = n_correct + 1;
      msg = '[correct]';
   else
      msg = ['[incorrect - ' num2str(id) ']'];
   end
   n_total = n_total + 1;
   score = n_correct/n_total;
   % display result
   disp(['test image ' num2str(n_total) ' of ' num2str(n_test) ': ' test.filenames{n} ' ' msg ' [score = ' num2str(score) ']']);
end
