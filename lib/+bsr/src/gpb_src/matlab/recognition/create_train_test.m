% split the dataset into train and test sets
function [train, test] = create_train_test( ...
   dirname, ...
   dirname_img, ...
   n_train_per_class, ...
   n_test_per_class)

classdirs = dir(dirname);
classdirs = classdirs(3:end);
classdirs = {classdirs.name};
n_classes = length(classdirs);

train_filenames = cell([1 n_train_per_class*n_classes]);
train_filenames_img = cell([1 n_train_per_class*n_classes]);
train_class_id  = zeros([1 n_train_per_class*n_classes]);
train_model_id  = zeros([1 n_train_per_class*n_classes]);

test_filenames  = cell([1 n_test_per_class*n_classes]);
test_filenames_img = cell([1 n_test_per_class*n_classes]);
test_class_id   = zeros([1 n_test_per_class*n_classes]);
test_model_id   = zeros([1 n_test_per_class*n_classes]);

train_num = 1;
test_num = 1;

for n = 1:n_classes
    files = dir([dirname '/' classdirs{n}]);
    files = files(3:end);
    files = {files.name};
    inds = randperm(length(files));
    train_files = files(inds(1:n_train_per_class));
    test_files  = files(inds((n_train_per_class+1):(n_train_per_class+n_test_per_class)));
    for f = 1:length(train_files)
        train_filenames{train_num} = [dirname '/' classdirs{n} '/' train_files{f}];
        img_prefix = train_files{f};
        img_prefix = img_prefix(1:10);
        train_filenames_img{train_num} = [dirname_img '/' classdirs{n} '/' img_prefix '_r.jpg'];
        train_class_id(train_num) = n;
        train_model_id(train_num) = f;
        train_num = train_num + 1;
    end
    for f = 1:length(test_files)
        test_filenames{test_num} = [dirname '/' classdirs{n} '/' test_files{f}];
        img_prefix = test_files{f};
        img_prefix = img_prefix(1:10);
        test_filenames_img{test_num} = [dirname_img '/' classdirs{n} '/' img_prefix '_r.jpg'];
        test_class_id(test_num) = n;
        test_model_id(test_num) = f;
        test_num = test_num + 1;
    end
end

train.filenames = train_filenames;
train.filenames_img = train_filenames_img;
train.class_id = train_class_id;
train.model_id = train_model_id;
train.n_train_per_class = n_train_per_class;

test.filenames = test_filenames;
test.filenames_img = test_filenames_img;
test.class_id = test_class_id;
test.model_id = test_model_id;
test.n_test_per_class = n_test_per_class;
