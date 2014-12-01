num_nn = 1; item_limit = 1000;
[train, test] = create_train_test('/home/mmaire_data/caltech101_gb_pb_pca','/home/mmaire_data/caltech101_resized',15,15);
build_db(train,'pca');
[cmx, score] = run_test(train,test,'pca',num_nn,item_limit);
