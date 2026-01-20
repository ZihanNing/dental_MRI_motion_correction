%% Example script for clustering function
dim = 3; % dimensions (coordinates)
v = rand(500,dim);
idx = 0;
n_pc = 34; 
eq_n = 1;
split_clust = 1.3;
init_iter = 20;
realign = 1; 
drawfig = 1;
dot_size = [10,40];
[C,v_c,n_c] = centroid_fct(v,idx,n_pc,eq_n,split_clust,init_iter,realign,drawfig,dot_size);

[stateSample, XClustered, dClusters, dSamples] = kmeans(X, KTemp,'Options', options,'Distance', distMetric,'MaxIter',numIt,'Replicates',nReplicas);
