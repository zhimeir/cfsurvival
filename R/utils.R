####################
### plot the data ##
####################
#' plot_raw
#'
#' An internal function to plot the raw data
#'
#' @param data The data frame to be plotted. Consists (X,R,censored_T).
#' @param dir outpu directory.
#' @param plot a logical value indicating if plotting the figure in teh interactive session.
#' @param hist a logical value indicating if plotting the histogram.
#' @export

plot_raw <- function(data,dir,plot =FALSE,hist=FALSE){

# plot the uncensored T
if(hist){
phist <- ggplot(data,aes(x=T))+
  geom_histogram()+
  theme_bw()+
  ggtitle("Histogram of uncensored survival time")+
  theme(
    axis.text = element_text(size=15),
    title = element_text(size=15,face="bold")
  )
if(plot){print(phist)}
filename <- sprintf("%s/histT.pdf",dir)
ggsave(filename,phist,width=8,height=3)

# plot censored T
pchist <- ggplot(data,aes(x=censored_T))+
  geom_histogram()+
  theme_bw()+
  ggtitle("Histogram of (randomly) censored survival time")
  theme(
    axis.text = element_text(size=15),
    title = element_text(size=15,face="bold")
  )
if(plot){print(pchist)}
filename <- sprintf("%s/histT_censored.pdf",dir)
ggsave(filename,pchist,width=8,height=3)
}
## plot the true data and the censored data
data <- data%>%
  rename(type = event)%>%
  mutate(type=ifelse(type==TRUE,
                     "censored",
                     "not censored"))
data$title <- "Censored T"
pp <- ggplot(data)+geom_point(aes(x=X,y=censored_T,col = type),size=.5,alpha=.3)+
  theme_bw()+
  xlab("X")+
  ylab("Censored_T")+
  ylim(c(0,max(data$censored_T)*2))+
  scale_color_manual(values = c("gray","red"))+
  theme(
    legend.position = c(0.9,0.8)
  )+
facet_grid(.~title)
if(plot){print(pp)}
filename <- sprintf("%s/censoredT.pdf",dir)
ggsave(filename,plot=pp,width=7,height=3)

data$title <- "Uncensored T"
pp <- ggplot(data)+geom_point(aes(x=X,y=T),size=.3)+
  theme_bw()+
  xlab("X")+
  ylab("T")+
  facet_grid(.~title)
if(plot){print(pp)}
filename <- sprintf("%s/T.pdf",dir)
ggsave(filename,plot=pp,width=7,height=3)
}


##############################
### plot result v.s. truth ###
##############################
plot_res_truth_weibull <- function(x_base,r,alpha,a,b,beta,plot=FALSE,version = "new"){
get_true_quant <- function(x,alpha,a,b,beta){return(qweibull(alpha,a,b*exp(-x*beta/a)))}

true_lo <- sapply(x_base,get_true_quant,a=a,b=b,alpha=alpha,beta = beta)
## true_up <- sapply(x_base,get_true_quant,a=a,b=b,alpha=(1+alpha)/2,beta = beta)

if(version == "new"){
  cf <- cfsurv(x_base,r=r,data,alpha)
  cf_lo <- unlist(cf[1,])
}else{
  cf <- cfsurv_quantreg(x_base,r=r,data,1-alpha)
  cf_lo <- unlist(cf[1,])
}
  ##   cf_up <- unlist(cf[2,])
  df <- data.frame(x=rep(x_base,2),y=c(true_lo,cf_lo),Curve = rep(c("True lower bnd","cf lower bnd"),each=length(x_base)))
pp <- ggplot(df)+geom_line(aes(x=x,y=y,col=Curve))+
  theme_bw()+
  ylab("Confidence bound")+
  geom_hline(yintercept = r, linetype="dashed",col="black")+
  ggtitle(sprintf("alpha=%.2f",alpha))
if(plot){print(pp)}
filename <- sprintf("%s/true_ci_comparison_alpha_%.2f_%s.pdf",dir,alpha,version)
ggsave(filename,pp,width=7,height=3)

}


#########################
### plot the result ####
#########################
plot_res <- function(x_base,r,alpha_list,
                            upp_ci=NULL,
                            low_ci,
                            method_name,
                            lower_only = FALSE,
                            includeR=NULL,
                            dir,plot){
if(lower_only){
  ##   y_max <- includeR*(r+5)+(1-includeR)*r
  df <- data.frame(x = rep(x_base,times=length(alpha_list)),
                   y_max=r,
                   y_min = low_ci,
                   CI_level = as.factor(rep(alpha_list,each=length(x_base))),
                   title = sprintf("%s:r=%d",method_name,r))

  pp <- ggplot(data = df,aes(x=x,ymin=y_min,ymax=y_max,group = CI_level,col = CI_level))+
  #geom_point(size = 1,shape=17)+
  #geom_line()+
  geom_ribbon(alpha=0.1)+
  theme_bw()+
  ylab("Survival time confidence interval")+
  theme(
    axis.text = element_text(size=15),
    strip.text = element_text(size=15),
    title = element_text(size=15)
  )+
  facet_grid(.~title)+
  ylim(c(0,r))
}else{
df <- data.frame(x = rep(x_base,times=length(alpha_list)),y_max=upp_ci,y_min = low_ci,CI_level = as.factor(rep(alpha_list,each=length(x_base))))
pp <- ggplot(data = df,aes(x=x,ymin=y_min,ymax=y_max,group = CI_level,col = CI_level))+
  #geom_point(size = 1,shape=17)+
  #geom_line()+
  geom_ribbon(alpha=0.1)+
  theme_bw()+
  ggtitle(paste0(method_name,":r=",as.character(r)))+
  ylab("Survival time confidence interval")+
  theme(
    axis.text = element_text(size=15),
    title = element_text(size=15,face="bold")
  )
}
if(plot){print(pp)}
filename <- sprintf("%s/%s_%d.pdf",dir,method_name,r)
ggsave(filename,pp,width=10,height=5)
}

##############################
### generate weibull T #######
##############################
#' generate_weibull
#'
#' A generic function to generate survival data generated from a cox model
#'
#' @export
generate_weibull <- function(n,a,b,beta,Rmax,seed,info=FALSE){
  set.seed(seed)
  X      <- runif(n,0,0.5)       # covariate
  Xbeta  <- beta*X 
  U <- runif(n, 0, 1)       # uniformly distributed RV
  T <- ceiling((-log(U)*b^a*exp(-Xbeta))^(1/a))   # simulated event time
  maxT <- max(T)
  R  <- sample(1:Rmax,n,replace = TRUE)# observing time

  # generate a data frame
  data <- data.frame(T=T,R=R,X=X,event = (T<R),censored_T = pmin(T,R))
  mortality <- sum(data$event)/n
  if(info) cat(sprintf("The mortality rate of this data set is: %.2f.\n",mortality))
  return(data)
}


##############################
### Check coverage     #######
##############################
check_coverage <-  function(res,T,R,limit=TRUE){
  cover <- 0
  n <- length(R)
  for(i in 1:n){
    if(limit){ 
      statement  <-  (T[i]<=res[[i]]$ci_upp)&(T[i]>=res[[i]]$ci_low)
    }else{
      statement <- (T[i]%in%res[[i]]$CI)
    }

    if(statement){
        cover = cover+1
    }else{
      if(res[[i]]$includeR == 1 & (T[i]>=R[i])){
        cover = cover+1
      }
    }
  }

  cover <- cover/n
  return(cover)

}

#' check_lower_coverage
#'
#' @param 
check_coverage_lower <-  function(res,T,R){
  cover <- 0
  n <- length(R)
  for(i in 1:n){
    statement  <-  (T[i]>=res[,i]$ci_low)*(T[i]<=R[i])
    if(statement){
        cover = cover+1
    }else{
      if(res[,i]$includeR == 1 & (T[i]>R[i])){
        cover = cover+1
      }
    }
  }

  cover <- cover/n
  return(cover)

}


plot_coverage <- function(cover,
                          alpha_list,
                          dir,
                          plotname,
                          plot = FALSE){
df <- data.frame(x = alpha_list,
                 y = cover)

pp <- ggplot(df,aes(x=x,y=y))+
  geom_point(shape=17,alpha=0.8,col="red")+
  geom_line(alpha=0.5,col = "red")+
  geom_abline(slope = 1,intercept=0,linetype="dashed",col="gray",alpha=.4)+
  theme_bw()+
  ##   scale_color_manual(values = c("blue","red"))+
  xlim(c(min(alpha_list),max(alpha_list)))+
  ylim(c(min(min(cover),min(alpha_list)),max(max(cover),max(alpha_list))))+
  ylab("Covering prob.")+
  xlab("Target level")
if(plot){print(pp)}
ggsave(paste0(dir,"/coverrate_",plotname,".pdf"),pp,width=5,height=3)
}


