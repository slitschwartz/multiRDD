
# Multisite RDD function library
#
# (C) 2021 Sophie Litschwartz, Luke Miratrix
#




#' LLR
#'
#' Before running this method, your data needs to be subsetted to just in
#' bandwidth obs and small clusters dropped before.
#'
#' @param df The dataset to analyze.
#'
#' @export
reg.LLR<-function(df, clusterid, cutpoint=0, runvar, outcomevar, less.than=FALSE) {

    #fix data frame to be the right format if using map_df
    if(typeof(df)=="list"){df<-df[]}

    df<- df %>% dplyr::mutate(obsid=row_number())
    #turn the data to a pdata frame to cluster the errors
    df.pl<-pdata.frame(df,index = c(clusterid,"obsid"))
    #drop all the variables we don't need
    df.pl <- df.pl %>% dplyr::select(all_of(clusterid),all_of(runvar),all_of(outcomevar),obsid)
    df.pl<- df.pl %>% dplyr::rename(cl.id=all_of(clusterid),running=all_of(runvar),Y=all_of(outcomevar))

    #create treatment variable
    if(less.than==T){df.pl$treat=as.numeric(df.pl$running<cutpoint)}else{
        df.pl$treat=as.numeric(df.pl$running>cutpoint)}

    reg<-plm(Y ~0+running+treat+I(treat*running)+factor(cl.id),data=df.pl,model="pooling")

    #this gets the clustered errors
    er.cl<-coeftest(reg,vcov=vcovHC(reg,type="sss",cluster="group"))

    #output data frame
    data.frame(Tr.Effect=as.numeric(coef(reg)["treat"]),
               Tr.Effect.se=er.cl["treat","Std. Error"],
               Tr.Effect.ll=confint(er.cl)["treat",1],
               Tr.Effect.ul=confint(er.cl)["treat",2],
               sample.size=nrow(df.pl),
               total.clusters=length(unique(df.pl$cl.id))
    )

}

#' The Meta approach
#'
#' Function to estimate homoskedastic standard errors for the individual site level effects
#'
#' @export
se.meta<-function(X,sig.bar){
    #print(X)
    X<-as.matrix(X)
    #print(X)
    XX<-t(X) %*% X
    #use mean squared error to get school level se's on the coefs
    er.var<-solve(XX)*(sig.bar^2)
    as.data.frame(t(sqrt(diag(er.var))))
}

#' Q-Stat Inverstion Function
#'
#' @param bj List of estimated impacts
#' @param vj List of associated standard errors.
#'
#' @export
Q.Stat.Inv<-function(bj,vj){

    stopifnot( length(bj) == length(vj) )

    #test tau values from 0 to 1
    #test tau is sd and not the variance
    tau_test <- seq(0,1,.00001)

    #degrees freedom
    deg.fr=length(bj)-1

    ## test two-sided version
    lowbound <- qchisq(0.025,deg.fr)
    low80 <- qchisq(0.10,df=deg.fr)

    highbound <- qchisq(0.975,df=deg.fr)
    high80 <- qchisq(0.90,df=deg.fr)

    #get a precision weight average of site level treatment effects

    #inverse precision weights
    wj <- 1/vj
    #precision weighted average of site level effects
    bbar <- sum(wj*bj)/(sum(wj))

    if( lowbound < highbound ){
        ## intialize values
        q_invert <- c()
        #print(vj)
        for (i in 1:length(tau_test)){
            #calculate new denominator that takes into account tau^2
            denom <- vj + tau_test[i]^2
            #new q-stat
            q_invert[i] <- sum((bj - bbar)^2/denom)
        }

        #compare q.stat to lower and upperbounds (Weiss et al, JREE p. 55)
        inCI_95 <- (q_invert>=lowbound & q_invert<=highbound)
        inCI_80 <- (q_invert>=low80 & q_invert<=high80)

        #get p-values to find variance with highest p-value
        pvals<-.5-abs(pchisq(q_invert,deg.fr)-.5)

        #get 95 upper limit and lower limit and
        #80 upper limit and lower limit and
        #create flags for if no variation or too much variation

        if ( max( q_invert ) < lowbound ) {
            # no evidence of any variation
            tau.est.ll<- 0
            tau.est.ul <- 0
            tau.est.ll80<- 0
            tau.est.ul80<- 0
            flag_novar <- 1
            flag_toovar <- 0
        } else if ( min( q_invert ) > highbound ) {
            # evidence of extreme variation, tau grid not large enough
            tau.est.ll<- 1
            tau.est.ul <- 1
            tau.est.ll80<- 1
            tau.est.ul80<- 1
            flag_novar <- 0
            flag_toovar <- 1
        } else {
            tau.est.ll<-min(tau_test[inCI_95])
            tau.est.ul<-max(tau_test[inCI_95])
            tau.est.ll80<-min(tau_test[inCI_80])
            tau.est.ul80<-max(tau_test[inCI_80])
            flag_novar <- 0
            flag_toovar <- 0
        }
    }else{ tau.est.ll<-NA
    tau.est.ul<-NA
    tau.est.ll80<- NA
    tau.est.ul80<- NA
    flag_novar <- NA
    flag_toovar <- NA }

    output<-c(tau.est.ll,tau.est.ul,tau.est.ll80,tau.est.ul80,flag_novar,flag_toovar)
    names(output)<-c("tau.est.ll","tau.est.ul","tau.est.ll80","tau.est.ul80","flag_novar","flag_toovar")
    return(output)
}


#' Conduct meta-analysis of multi RDD
#'
#' drop small clusters and subset to only data in bw before running data
#'
#' @inheritParams reg.LLR
#'
#' @export
reg.Meta<-function(df,cutpoint=0,less.than=FALSE,clusterid,runvar,outcomevar){

    #fix data frame to be the right format if using map_df
    if(typeof(df)=="list"){df<-df[]}

    df.meta.1<- df %>% dplyr::rename(cl.id=all_of(clusterid),running=all_of(runvar),
                                     Y=all_of(outcomevar))
    #create treatment variable
    if(less.than==T){df.meta.1$treat=as.numeric(df.meta.1$running<cutpoint)}else{
        df.meta.1$treat=as.numeric(df.meta.1$running>cutpoint)}

    #drop any phantom clusters
    df.meta.1$cl.id<-droplevels(df.meta.1$cl.id)

    #turn data into a list of schools and then run a linear regression in each school
    df.out.1<-df.meta.1 %>%
        #turn data set into a list of schools
        split(.$cl.id) %>%
        #run the regression by school, treat, running, and treat by running interaction in the regression
        map(~lm(Y~I(treat*running)+treat+running,data=.x)) %>%
        #pull out into a data frame for each school  the all
        #coefficients, standard errors, residuals, and number of obs
        map_df(~ data.frame(t(as.matrix(coef(.))),se=t(as.matrix(se.coef(.))),
                            sighat=sigma.hat(.),nobs=nobs(.)),.id="regnum")

    #recalculate the se's assuming homoskedacity to improve precision of the precisions
    #get weighted mean residual standard error
    mrse<-sqrt(sum(with(df.out.1,(sighat)^2*nobs))/sum(df.out.1$nobs))

    #this recalculates the se's
    df.out.2<-df.meta.1 %>%
        dplyr::select(treat,running,cl.id) %>%
        #this creates an X matrix
        mutate(Intercept=1,treat.running=I(treat*running)) %>%
        #this splits by school
        split(.$cl.id) %>%
        #this calculates the erros with the X matrices,
        #but uses mrse to get the se's for each school
        map_df(~se.meta(subset(.x,select=-cl.id),mrse),.id="regnum")

    #this subsets to just grab the homoskedastic se for the treatment effect
    df.out.3<-df.out.2 %>%
        dplyr::rename("treat.se.ho"=treat) %>%
        dplyr::select(treat.se.ho,regnum)

    #Calculate the treatment effect variance using the DL estimator

    #calculate Q-statistic
    #site level treatment effects
    bj <- df.out.1$treat
    #site level se
    vj <- (df.out.3$treat.se.ho)^2

    #inverse precision weights
    wj.1 <- 1/vj
    #precision weighted average of site level effects
    bbar.1 <- sum(wj.1*bj)/(sum(wj.1))

    #Q stat
    Q.1 <- sum((bj - bbar.1)^2/vj)


    #degrees freedom
    deg.fr=length(bj)-1
    #quartic term in the DL ormula
    w4j.1<-1/(vj^2)

    #DL Numerator
    numerator.1<-Q.1-deg.fr
    #DL Demoninator
    denominator.1<-sum(wj.1)-(sum(w4j.1)/sum(wj.1))

    #The variance estimate from DL formula, if less than 0 than estimate is 0
    tau.est.Q=sqrt(max(0,numerator.1/denominator.1))

    #now do test inversion to get the confident interval

    ## get confidence interval
    tau.est.CI<-Q.Stat.Inv(bj,vj)

    #use tau estimate to get LATE estimate
    #get new weights that incorporate the variance
    wj.tre=1/(vj+tau.est.Q^2)
    #treatment effect
    tr.effect=sum(bj*wj.tre)/sum(wj.tre)

    #get average treatment effect se
    tr.effect.se=sqrt(1/sum(wj.tre))
    #output data frame
    data.frame(Tr.Effect=tr.effect, #this is the LATE
               Tr.Effect.se=tr.effect.se, #this is the standard error on the LATE
               Tr.Effect.ll=tr.effect - 1.96*tr.effect.se, #lower limit of 95% CI on LATE
               Tr.Effect.ul=tr.effect + 1.96*tr.effect.se, #upper limit of 95% CI on LATE
               Sigma.Est=tau.est.Q, #sd of LATE, DL estimator
               Sigma.LL.Q=tau.est.CI["tau.est.ll"], #lower limit of 95% CI on sd of LATE
               Sigma.UL.Q= tau.est.CI["tau.est.ul"], #upper limit of 95% CI on sd of LATE
               Sigma.LL.Q80= tau.est.CI["tau.est.ll80"], #lower limit of 80% CI on sd of LATE,
               Sigma.UL.Q80= tau.est.CI["tau.est.ul80"], #upper limit of 80% CI on sd of LATE
               flag_novar.Q=tau.est.CI["flag_novar"], #no variance flag
               flag_toovar.Q=tau.est.CI["flag_toovar"], #too much variance flag
               sample.size=nrow(df),
               total.clusters=length(unique(df.meta.1$cl.id)))
}

#' Conduct FIRC analysis of multi RDD
#'
#' @inheritParams reg.LLR
#' @export
reg.FIRC<-function(df,cutpoint=0,less.than=FALSE,clusterid,runvar,outcomevar,model){

    #inititialize no Q-stat variable
    noQ=F

    #save the data in a data set with correct variables and names
    df.reg<- df %>% dplyr::rename(cl.id=all_of(clusterid),running=all_of(runvar),
                                  Y=all_of(outcomevar)) %>%
        dplyr::select(cl.id,running,Y)

    #create treatment variable
    if(less.than==T){df.reg$treat=as.numeric(df.reg$running<cutpoint)}else{
        df.reg$treat=as.numeric(df.reg$running>cutpoint)}

    #break up the parts of the regression into saved objects,
    #so that the function works for different variable combinations in the model

    #pooled regression fixed effects
    regvars<-paste(0,"running","treat","I(treat*running)",sep="+")

    #unpooled indicators
    if(model["intercept"]!="up"){ups<-paste(1,sep="")}
    if(model["intercept"]=="up"){ups<-paste("factor(cl.id)",sep="")}
    if(model["running"]=="up"){ups<-paste(ups,"factor(cl.id):running",sep="+")}
    if(model["treat"]=="up"){ups<-paste(ups,"factor(cl.id):treat",sep="+")}
    if(model["tr.rn"]=="up"){ups<-paste(ups,"factor(cl.id):I(treat*running)",sep="+")}

    #random effects
    if(model["intercept"]!="pp"){res<-paste(0,sep="")}
    if(model["intercept"]=="pp"){res<-paste(1,sep="")}
    if(model["treat"]=="pp"){res<-paste(res,"treat",sep="+")}
    if(model["running"]=="pp"){res<-paste(res,"running",sep="+")}
    if(model["tr.rn"]=="pp"){res<-paste(res,"I(treat*running)",sep="+")}

    #print(regvars)
    #print(ups)
    #print(res)

    #run regression
    reg<-tryCatch(lmer(paste("Y~",regvars,"+",ups,"+(",res,"|cl.id)",sep=""),data = df.reg,REML=T),error=function(e) NA)
    if(is.na(reg)==F){
        #only get random effect info if there is a random effect on treatment
        if(model["treat"]=="pp"){
            sigma.tr<-data.frame(sigma.hat(reg)$sigma$cl.id["treat"])
            names(sigma.tr)<-c("sigma.est")
        }else{sigma.tr<-NULL}

        #get Q-Stat inversion CI's, this requires it's own model
        #becuase the site level effect estimates from the FIRC are over shrunk
        #no Q variable is true if there is no random effect on treatment
        if(model["treat"]=="pp"){
            #interacted variables (unpooled and pp)
            if(!(unlist(model["intercept"]) %in% c("up","pp"))){inter<-paste(1,sep="")}
            if(unlist(model["intercept"]) %in% c("up","pp")){inter<-paste("factor(cl.id)",sep="")}
            if(unlist(model["running"]) %in% c("up","pp")){inter<-paste(inter,"factor(cl.id):running",sep="+")}
            if(unlist(model["treat"]) %in% c("up","pp")){inter<-paste(inter,"factor(cl.id):treat",sep="+")}
            if(unlist(model["tr.rn"]) %in% c("up","pp")){inter<-paste(inter,"factor(cl.id):I(treat*running)",sep="+")}

            #uninteracted variables (pooled)
            #random effects
            if(model["intercept"]!="pool"){pools<-paste(0,sep="")}
            if(model["intercept"]=="pool"){pools<-paste(1,sep="")}
            if(model["treat"]=="pool"){pools<-paste(pools,"treat",sep="+")}
            if(model["running"]=="pool"){pools<-paste(pools,"running",sep="+")}
            if(model["tr.rn"]=="pool"){pools<-paste(pools,"I(treat*running)",sep="+")}
            if(model["tr.rn"]=="quad"){pools<-paste(pools,"I(running*running)",sep="+")}

            ols<-lm(paste("Y~",pools,"+",inter,sep=""),data = df.reg)

            #get treatment effects
            ols.coef<-as.data.frame(coef(summary(ols)))
            ols.coef$varname<-rownames(ols.coef)

            #get site level treatment effects and standard errors to run
            #Q-statistic inversion CIs
            bj<-as.numeric(unlist(ols.coef %>% filter(str_detect(varname, ':treat'))
                                  %>% dplyr::select("Estimate")))
            vj<-as.numeric(unlist((ols.coef %>% filter(str_detect(varname, ':treat')) %>% dplyr::select("Std. Error"))^2))

            #now run test inversion to get CI
            if(length(bj)>0 | length(vj)>0){
                sigma.tr.Q<-Q.Stat.Inv(bj,vj)}else{noQ=T}}else{noQ=T}
        if(noQ==T){
            sigma.tr.Q<-rep(NA,6)
            names(sigma.tr.Q)<-c("tau.est.ll","tau.est.ul","tau.est.ll80","tau.est.ul80",
                                 "flag_novar","flag_toovar" )
        }

        #get the CI for the treatment effect
        confint.fe<-tryCatch(confint<-confint(reg,c("treat"),method="Wald"),error=function(e) NA)

        #this onl creates CI variables if CI function worked
        if(!is.na(first(confint.fe))){
            alpha1.ll=confint.fe["treat",1]
            alpha1.ul=confint.fe["treat",2]
        }else{alpha1.ll<-NA
        alpha1.ul<-NA}

        temp.lmer<-data.frame(Tr.Effect=as.numeric(coef(summary(reg))["treat","Estimate"]),
                              Tr.Effect.se=as.numeric(coef(summary(reg))["treat","Std. Error"]),
                              Tr.Effect.ll=alpha1.ll,
                              Tr.Effect.ul=alpha1.ul,
                              Sigma.Est=sigma.tr, #sd of LATE
                              Sigma.LL.Q=sigma.tr.Q["tau.est.ll"], #lower limit of 95% CI on sd of LATE
                              Sigma.UL.Q= sigma.tr.Q["tau.est.ul"], #upper limit of 95% CI on sd of LATE
                              Sigma.LL.Q80= sigma.tr.Q["tau.est.ll80"], #lower limit of 80% CI on sd of LATE,
                              Sigma.UL.Q80= sigma.tr.Q["tau.est.ul80"], #upper limit of 80% CI on sd of LATE
                              flag_novar.Q=sigma.tr.Q["flag_novar"], #no variance flag
                              flag_toovar.Q=sigma.tr.Q["flag_toovar"], #too much variance flag
                              sing=isSingular(reg), #flag for if regression was singular
                              AIC=AIC(reg), #regression AIC
                              sample.size=nrow(df),
                              total.clusters=length(unique(df.reg$cl.id)))
    }else{
        temp.lmer<-data.frame(model,sample.size=nrow(df),totsch.inbw=first(df$totsch.inbw),
                              total.clusters=length(unique(df.reg$cl.id)))}
    return(temp.lmer)
    # reg
}
