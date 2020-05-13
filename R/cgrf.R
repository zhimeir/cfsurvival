list.of.packages <- c("grf")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
library(grf)

grf.getWeights = function(grf, testdata, y_name, c_name, x_name=NULL) {
  
  if (is.null(x_name)){
    weights <- get_sample_weights(grf, newdata=testdata[ ,!(names(testdata) %in% c(y_name, c_name)), drop=F]) 
  } else {
    weights <- get_sample_weights(grf, newdata=testdata[ ,x_name])
  }
  
  return (weights)
}

# example
#weights <- grf.getWeights(grf_qf, data_test, 'y', 'status')
#dim(weights)
#class(weights)
