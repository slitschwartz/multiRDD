
#Generate  Simulation Data Set

#Functions for Generating the Sample Data

#' Generate clusters
#'
#' @param pars are the sample mean and sd's for the underlying regression
#'   parameters. pars cotains means, sds, and names
#' @param K is the total clusters in the sample
#' @param ICC is the ICC of the running variable
#' @param njbar is the mean total observations in each cluster
#'
#' @export
gen.clusters<-function (pars,K,ICC,njbar){
    #draw parameters for each cluster
    clusters<-with( pars, data.frame( mvrnorm(n=K,mu=means,Sigma=diag(sds^2))) )
    #names comes from pars
    names(clusters)<-pars$names

    #draw the running variable cluster level change in the mean
    clusters$rj<-rnorm(K,0,sqrt(ICC))
    #puts the ICC into the data set
    clusters$ICC<-ICC
    #draw the total observations in each cluster, +1 so each cluster has at least one student
    clusters$Nj<-(rpois(K,njbar))+1
    #assign a cluster id for each school
    clusters$cl.id<-seq(1:K)
    #put total schools into the data set
    clusters$TC<-K
    #put school mean of total students into the data set
    clusters$nmean<-njbar
    return(clusters)
}

#' generate a sample of students using the school parameters
#'
#' assume a cut point of 0, treatment if below cut
#' @param pars are the parameters for each school from gen.schools
#' @param se is the student level standard error
#' @param rmean is the grand mean of the running variable
#'
#' @export
gen.obs<-function (a0, b0, a1, b1, rj, ICC, Nj, cl.id, TC, nmean, se, rmean){
    #print(pars)
    #get distribution of running variable in the school
    running<-rmean+rj+rnorm(Nj,0,sqrt((1-ICC)))
    #student id
    obs.id<-seq(1:Nj)
    #observation level error
    epsilon<-rnorm(Nj,0,se)
    #treatment status, assume cut point of 0, treated if less than cut
    treat<-1*(running<(0))
    #outcome value for control
    Y0<-a0+b0*running+epsilon
    #treament value
    Y1<-Y0+a1+b1*running
    #assign actual Y based on treatment status
    Y<-ifelse(treat==1,Y1,Y0)
    #combine observation id, running var, student error, treatment status,
    #potential and actual outcomes, and regression parameters into an output data set
    observations<-data.frame(obs.id,running,epsilon,treat,Y0,Y1,Y)
    return(observations)
    #print(pars)
}


#' Generate set of clusters with RDD in each clustrer.
#'
#' @param TC  total clusters (TC) to 150
#' @param nmean  average observations per cluster (nmean) to 130
#'
#' @export
gen.data<-function( TC = 150, nmean = 130,
                    parameters.reg,
                    obs.se, Rmean ) {

    #generate the sample as a data frame
    clusters<-gen.clusters(parameters.reg,
                           TC,
                           parameters.gen$ICCR,
                           nmean)

    #Generate the full data sample from the cluster samples
    #gen.students is the function that generates the students
    #stu.se is the student level standard error
    #Rmean is the grand mean of the running variable

    sim.sample = as_tibble( clusters )
    sim.sample$dat <-pmap(clusters, gen.obs,
                          se=obs.se, rmean=Rmean)
    sim.sample = unnest( sim.sample, cols="dat" )
    sim.sample

}
