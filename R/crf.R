list.of.packages <- c("foreach", "doParallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
library(foreach)
library(doParallel)

rf.getNodes = function(rf, data, y_name, c_name) {
  pred = predict(rf, data[ ,!(names(data) %in% c(y_name, c_name)), drop=F], nodes = T)
  nodes = attributes(pred)$nodes

  return(nodes)
}

rf.getWeights = function(rf, traindata, testdata, y_name, c_name) {
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  cat("number of cores: ", cores[1]-1)
  registerDoParallel(cl)

  # retrieve training nodes
  nodes = rf.getNodes(rf, traindata, y_name, c_name) # [train.samples, trees]

  # retrieve nodes for test data
  test.nodes = rf.getNodes(rf, testdata, y_name, c_name) # [test.samples, trees]
  nodesize = test.nodes # [test.samples, trees]

  # loop over trees
  print("loop begins...\n")
  ntrees = rf$ntree
  weights <- foreach(k=1:ntrees, .combine='+') %dopar% {
      in.leaf = 1*outer(test.nodes[,k],nodes[,k],'==') # [test.samples, train.samples] indicates whether a test and a train data are in the same leaf node
      mapping=plyr::count(nodes[,k])
      nodesize[,k] = plyr::mapvalues(test.nodes[,k], from=mapping$x, to=mapping$freq, warn_missing=F)
      weight = in.leaf / nodesize[,k]
      weight
  }
  stopCluster(cl)

  return(weights/ntrees)
}
