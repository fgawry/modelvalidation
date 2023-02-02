#' Extract model residuals
#'
#' Extract model specific results from `lm()`, `glm()`, `lmer()`,
#' `lme()` and `gls()` models.
#'
#' @param model A fitted model
#' @param bin.size For glm binomial family only. Binned residuals.
#' @return Model specific residuals.
extract.resid <-function(model, bin.size) {
if (class(model)[[1]] == "lm") {
  resid<-resid(model)
  ylab<-"Ord. residuals"
  return(list(resid=resid,ylab=ylab))
}
if (class(model)[[1]] == "gls") {
  resid<-resid(model, type="normalized")
  ylab<-"Stand. residuals"
  return(list(resid=resid,ylab=ylab))
}
if (class(model)[[1]] == "lme") {
  resid<-resid(model, type='pearson', level=1)
  nranef<-length(nlme::ranef(model))
  resid.random<-nlme::ranef(model)
  ylab<-"Pearson resid. Fixed (L1)"
  ylab.random<-"Pearson Resid.-Rand. Effects"
  return(list(resid=resid,
              nranef=nranef,
              resid.random=resid.random,
              ylab=ylab,
              ylab.random=ylab.random))

}
if (class(model)[[1]] == "lmerMod") {
  resid<-resid(model, type="pearson", level=1)
  ylab<-"Pearson resid. Fixed (L1)"
  nranef<-ncol(lme4::ranef(model)[[1]])
  resid.random<-lme4::ranef(model)[[1]]
  ylab.random<-"Pearson Resid.-Rand. Effects"
  return(list(resid=resid,
              nranef=nranef,
              resid.random=resid.random,
              ylab=ylab,
              ylab.random=ylab.random))
}
if (class(model)[[1]] == "glm") {

  family<-model$family[[1]]

  if(family=="gaussian" | family=="Gamma" |
     family=="inverse.gaussian" | family=="poisson" | family=="quasi") {
    resid<-resid(model, type="deviance")
    ylab<-"Deviance Residuals"
    return(list(resid=resid,ylab=ylab))
  }
  if(family=="quasipoisson" | family=="quasibinomial") {
    dispersion<-summary(model)$dispersion
    mu <- predict(model, type = "response")
    E <- model$y - mu
    resid <- E / sqrt(dispersion * mu)
    ylab<-"Scaled Pearson Residuals"
    return(list(resid=resid,ylab=ylab))
  }
  if(family=="binomial") {
    resid<-resid(model, type="pearson")
    ylab<-"Pearson Residuals"
    ylab.bin<-"Pearson Residuals (bins)"

    bin.size<-bin.size
    bins<-round(seq(from=1,to=length(resid), length.out=length(resid)/bin.size))
    bin.factor<-vector(length=length(resid))
    resid.mean<-vector(length = length(bins)-1)
    #resid.sd<-resid
    if (length(bins)<=1) {stop("less than two bins. Decrease bin size.")}
    for (i in 1:(length(bins)-1)){
      bin.factor[bins[[i]]:bins[[i+1]]]<-rep(i, length(bins[[i]]:bins[[i+1]]))
      resid.mean[[i]]<-mean(resid[bins[[i]]:bins[[i+1]]])
      #resid.sd[[i]]<-sd(resid[bins[[i]]:bins[[i+1]]])
    }
    bin.factor<-as.factor(bin.factor)
    return(list(resid=resid,
                ylab=ylab,
                ylab.bin=ylab.bin,
                bins=bins,
                resid.mean=resid.mean,
                bin.factor=bin.factor))

  }
}
}




#' Extract model data
#'
#' Extract model data from `lm()`, `glm()`, `lmer()`,
#' `lme()` and `gls()` models.
#'
#' @param model A fitted model
#' @return Model data.
extract.data <-function(model) {
  if (class(model)[[1]] == "lm") {
    data<-model$model
    return(data)
  }
  if (class(model)[[1]] == "gls") {
    data<-nlme::getData(model)
    return(data)
    if(is.null(data)==TRUE) {
      stop("Model has no data.")
    }
  }
  if (class(model)[[1]] == "lme") {
    data<-nlme::getData(model)
    return(data)
    if(is.null(data)==TRUE) {
      stop("Model has no data.")
    }
  }
  if (class(model)[[1]] == "lmerMod") {
    data<-model.frame(model)
    return(data)
  }
  if (class(model)[[1]] == "glm") {
    data<-model$model
    return(data)
  }
}





#' Model simulation
#'
#'Simulate new models based on `lm()`, `glm()` or `lmer()` models.
#'See also `simulate()` for further details.
#'
#' @param model A fitted model
#' @param nsim Number of simulations
#' @param data Original model data
#' @return Simulated models
insp.sim <-function(model,
                    nsim=19,
                    data=extract.data(model)) {
  lmer.sim <- list()
  if (class(model)[[1]] == "lmerMod" | class(model)[[1]] ==
      "glm" | class(model)[[1]] == "lm") {
    lmer.sim.data <- simulate(object = model, nsim = nsim)
    for (i in 1:ncol(lmer.sim.data)) {
      new_data<-data
      new_data$lmer.sim.data_y <- lmer.sim.data[, i]
      lmer.sim[[i]] <- update(model, lmer.sim.data_y ~
                                ., data = new_data)
    }
  }
  if (class(model)[[1]] == "gls") {
    print("These simulations may not work for correlation structures")
    lmer.sim.data <-vector("list", length = nsim)
    new_data<-data
    for (i in 1:nsim) {
      print(paste("Simulation #",i))
      lmer.sim.data[[i]] <- nlraa::simulate_gls(model, psim = 2)
      new_data$lmer.sim.data_y <- as.vector(lmer.sim.data[[i]])
      lmer.sim.data[[i]] <- new_data
      lmer.sim[[i]] <- update(model, lmer.sim.data_y ~
                                ., data = new_data)
    }
  }

  if (class(model)[[1]] == "lme") {
    print("These simulations may not work for correlation structures")
    lmer.sim.data <-vector("list", length = nsim)
    new_data<-data
    for (i in 1:nsim) {
      print(paste("Simulation #",i))
      lmer.sim.data[[i]] <- nlraa::simulate_lme(model, psim = 3)
      new_data$lmer.sim.data_y <- as.vector(lmer.sim.data[[i]])
      lmer.sim.data[[i]] <- new_data
      lmer.sim[[i]] <- update(model, lmer.sim.data_y ~
                                ., data = new_data)
    }
  }



  if (class(model)[[1]] == "gls" | class(model)[[1]] ==
      "lme") {
    lmer.sim[[length(lmer.sim) + 1]] <- model
    lmer.sim.data[[length(lmer.sim.data) + 1]] <- data
    r<-list(sim=lmer.sim, sim.data=lmer.sim.data)
    return(r)
  }

  if (class(model)[[1]] == "lmerMod" | class(model)[[1]] ==
      "glm" | class(model)[[1]] == "lm") {
    lmer.sim[[length(lmer.sim) + 1]] <- model
    return(lmer.sim)
  }


}

#' Plot simulated data
#'
#' Plot a series of validation plots based on results from
#' `insp.sim()` function.
#'
#' @param insp.sim A fitted model simulation list form `insp.sim`
#' @param resid.sim A list of model simulation residuals
#' @param fitted.sim A list of model simulation fitted data
#' @param data.sim A list of model simulation data
#' @param ask.sim Whether ask before showing real data plot
#' @examples
#' n<-100
#' pred1<-rnorm(n=n, mean=0, sd=1)
#' pred2<-rnorm(n=n, mean=1, sd=2)
#' pred3<-gl(n=5,k=n/5)
#' ERROR<-rnorm(n=n, mean=0, sd=2)
#' a<-1
#' b<-2
#' c<-rep(runif(n=length(levels(pred3)),min = -2, max = 2),
#'        each=n/length(levels(pred3)))
#' Y <- a*pred1 + b*pred2 + c + ERROR
#' lm<-lm(Y ~ pred1 + pred2 + pred3)
#' sim<-insp.sim(model=lm)
#' plots.insp.sim(insp.sim=sim,
#'               ask=FALSE)
plots.insp.sim<-function(insp.sim,
                         resid.sim=lapply(insp.sim,extract.resid),
                         fitted.sim=lapply(insp.sim,fitted),
                         data.sim=lapply(insp.sim,extract.data),
                         pred=NULL,
                         ask.sim=TRUE) {

  if(is.null(pred)==TRUE) {
    pred<-colnames(data.sim[[1]])[2:ncol(data.sim[[1]])]
  }

  ask.next<-function() {
    utils::askYesNo("Show next validation plot?")
  }

  if (class(insp.sim[[1]])[[1]]=="lm") {

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.hist.r(resid=resid.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.qqnorm.r(resid=resid.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.resid.x.fitted.r(resid=resid.sim,fitted=fitted.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.resid.x.pred.r(resid=resid.sim,
                           data=data.sim,
                           pred=pred,
                           ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.cook.r(models=insp.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.acf.r(models=insp.sim, ask=ask.sim)
    }
  }

  if (class(insp.sim[[1]])[[1]]=="glm") {

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.hist.r(resid=resid.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.qqnorm.r(resid=resid.sim,ask=ask.sim)
    }

    if(family(insp.sim[[1]])[[1]]=="binomial") {
      reply<-ask.next()
      if (is.na(reply)==TRUE) {stop("Cancelled.")}
      if (reply==TRUE) {
        plots.hist.bin.r(resid=resid.sim,ask=ask.sim)
      }
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.resid.x.fitted.r(resid=resid.sim,fitted=fitted.sim,ask=ask.sim)
    }

    if(family(insp.sim[[1]])[[1]]=="binomial") {
      reply<-ask.next()
      if (is.na(reply)==TRUE) {stop("Cancelled.")}
      if (reply==TRUE) {
        plots.resid.bin.r(resid=resid.sim,ask=ask.sim)
      }
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.resid.x.pred.r(resid=resid.sim,
                           data=data.sim,
                           pred=pred,
                           ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.cook.r(models=insp.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.acf.r(models=insp.sim,ask=ask.sim)
    }
  }

  if (class(insp.sim[[1]])[[1]]=="lmerMod") {

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.hist.r(resid=resid.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.qqnorm.r(resid=resid.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.resid.x.fitted.r(resid=resid.sim,fitted=fitted.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.resid.x.pred.r(resid=resid.sim,
                           data=data.sim,
                           pred=pred,
                           ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.qqnorm.random.r(resid=resid.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.acf.r(models=insp.sim,ask=ask.sim)
    }
  }

  if (class(insp.sim[[1]])[[1]]=="gls") {

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.hist.r(resid=resid.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.qqnorm.r(resid=resid.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.resid.x.fitted.r(resid=resid.sim,
                             fitted=fitted.sim,
                             ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.resid.x.pred.r(resid=resid.sim,
                           data=data.sim,
                           pred=pred,
                           ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.acf.r(models=insp.sim,
                  ask=ask.sim)
    }
  }

  if (class(insp.sim[[1]])[[1]]=="lme") {

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.hist.r(resid=resid.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.qqnorm.r(resid=resid.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.resid.x.fitted.r(resid=resid.sim,fitted=fitted.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.resid.x.pred.r(resid=resid.sim,
                           data=data.sim,
                           pred=pred,
                           ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.qqnorm.random.r(resid=resid.sim,ask=ask.sim)
    }

    reply<-ask.next()
    if (is.na(reply)==TRUE) {stop("Cancelled.")}
    if (reply==TRUE) {
      plots.acf.r(models=insp.sim,ask=ask.sim)
    }
  }

}

#' Asks whether actual datum plot sould be shown.
#'
#' This function is used internally.
#'
#' @param real.data plot position with real data
#' @param ask Whether ask before showing real data plot
ask.reveal.data<-function(real.data, ask) {
  if (ask==TRUE) {
    ask1<-utils::askYesNo(msg="Show real data plot?")
    if (ask1==TRUE) {
      screen(n=real.data, new=FALSE)
      box(col="blue")
    }
    if (is.na(ask1)==TRUE){stop("Cancelled")}
  } else {
    screen(n=real.data, new=FALSE)
    box(col="blue")
  }
}


#' Generates a plot random order sequence.
#'
#' This function is used internally.
#'
#' @param n number of plots
#' @param datum plot number with real data
random.plots<-function(n=20,datum=20) {
  random.plots<-sample(1:n)
  real.data<-match(datum,random.plots)
  return(list(random.plots=random.plots,real.data=real.data))
}

#' Plot residuals histograms from `insp.sim()` data.
#'
#'
#' @param resid list with models residuals
#' @param ask whether ask before showing real data plot
plots.hist.r<-function(resid, ask) {
  random.plots<-random.plots()
  dev.new()
  dev.off()
  #old.par<-par()
  close.screen()
  split.screen(figs=c(5,4))
  new.par<-list(mar = c(1, 1, 1, 1),
                mgp=c(1,0.1,0),
                tcl=0.5,
                ann=FALSE,
                col.axis="black",
                cex.axis=1,
                xpd=TRUE,
                cex.main=1,
                cex=0.75)
  for (i in 1:length(resid)) {
    screen(i)
    par(new.par)
    hist(resid[[ random.plots$random.plots[[i]] ]]$resid,
         main=NULL, xlab=NA,ylab=NA)
    if (i==1) {
      title(main=resid[[ random.plots$random.plots[[i]] ]]$ylab)
    }
  }
  ask.reveal.data(real.data=random.plots$real.data, ask=ask)
  #par(old.par)
}

#' Plot residuals histograms from `insp.sim()` binomial data.
#'
#'
#' @param resid list with models residuals
#' @param ask whether ask before showing real data plot
plots.hist.bin.r<-function(resid, ask) {
  random.plots<-random.plots()
  dev.new()
  dev.off()
  #old.par<-par()
  close.screen()
  split.screen(figs=c(5,4))
  new.par<-list(mar = c(1, 1, 1, 1),
                mgp=c(1,0.1,0),
                tcl=0.5,
                ann=FALSE,
                col.axis="black",
                cex.axis=1,
                xpd=TRUE,
                cex.main=1,
                cex=0.75)
  for (i in 1:length(resid)) {
    screen(i)
    par(new.par)
    hist(resid[[ random.plots$random.plots[[i]] ]]$resid.mean,
         main=NULL, xlab=NA,ylab=NA)
    if (i==1) {
      title(main=resid[[ random.plots$random.plots[[i]] ]]$ylab.bin)
    }
  }
  ask.reveal.data(real.data=random.plots$real.data, ask=ask)
  #par(old.par)
}

#' Plot residuals qqnorm from `insp.sim()` data.
#'
#'
#' @param resid list with models residuals
#' @param ask whether ask before showing real data plot
plots.qqnorm.r<-function(resid, ask) {
  random.plots<-random.plots()
  dev.new()
  dev.off()
  #old.par<-par()
  close.screen()
  split.screen(figs=c(5,4))
  new.par<-list(mar = c(1, 1, 1, 1),
                mgp=c(1,0.1,0),
                tcl=0.5,
                ann=FALSE,
                col.axis="black",
                cex.axis=1,
                xpd=FALSE,
                cex.main=1,
                cex=0.75)
  for (i in 1:length(resid)) {
    screen(i)
    par(new.par)
    qqnorm(resid[[ random.plots$random.plots[[i]] ]]$resid,
           main = NULL,
           ylab = NA,
           xlab = NA)
    qqline(resid[[ random.plots$random.plots[[i]] ]]$resid, col="red")
    if (i==1) {
      title(main=resid[[random.plots$random.plots[[i]] ]]$ylab)
    }
  }
  ask.reveal.data(real.data=random.plots$real.data, ask=ask)
}

#' Plot residuals x fitted from `insp.sim()` data.
#'
#'
#' @param resid list with models residuals
#' @param fitted list with models fitted values
#' @param ask whether ask before showing real data plot
plots.resid.x.fitted.r<-function(resid,fitted,ask) {
  random.plots<-random.plots()
  dev.new()
  dev.off()
  #old.par<-par()
  close.screen()
  split.screen(figs=c(5,4))
  new.par<-list(mar = c(1, 1, 1, 1),
                mgp=c(1,0.1,0),
                tcl=0.5,
                ann=FALSE,
                col.axis="black",
                cex.axis=1,
                xpd=FALSE,
                cex.main=1,
                cex=0.75)
  for (i in 1:length(resid)) {
    screen(i)
    par(new.par)
    plot(resid[[ random.plots$random.plots[[i]] ]]$resid~
           fitted[[ random.plots$random.plots[[i]] ]],
           main = NULL,
           ylab = NA,
           xlab = NA)
    lines(lowess(resid[[ random.plots$random.plots[[i]] ]]$resid~
                   fitted[[ random.plots$random.plots[[i]] ]]),col="red")
    if (i==1) {
      title(main="Resid. x Fitted")
    }
  }
  ask.reveal.data(real.data=random.plots$real.data,ask=ask)
}

#' Plot residuals from `insp.sim()` binomial data.
#'
#'
#' @param resid list with models residuals
#' @param ask whether ask before showing real data plot
plots.resid.bin.r<-function(resid,ask) {
  random.plots<-random.plots()
  dev.new()
  dev.off()
  #old.par<-par()
  close.screen()
  split.screen(figs=c(5,4))
  new.par<-list(mar = c(1, 1, 1, 1),
                mgp=c(1,0.1,0),
                tcl=0.5,
                ann=FALSE,
                col.axis="black",
                cex.axis=1,
                xpd=FALSE,
                cex.main=1,
                cex=0.75)
  for (i in 1:length(resid)) {
    screen(i)
    par(new.par)
    plot(resid[[ random.plots$random.plots[[i]] ]]$resid~
           resid[[ random.plots$random.plots[[i]] ]]$bin.factor,
         main = NULL,
         ylab = NA,
         xlab = NA)
    lines(lowess(resid[[ random.plots$random.plots[[i]] ]]$resid~
                   resid[[ random.plots$random.plots[[i]] ]]$bin.factor),col="red")
    if (i==1) {
      title(main="Resid. x Bins")
    }
  }
  ask.reveal.data(real.data=random.plots$real.data,ask=ask)
}

#' Plot residuals x all covariates from `insp.sim()` data.
#'
#'
#' @param resid list with model residuals
#' @param data list with model data
#' @param ask whether ask before showing real data plot
plots.resid.x.pred.r<-function (resid,
                               data,
                               pred,
                               ask) {
  n.pred<-length(pred)
  for (k in 1:n.pred) {
    random.plots<-random.plots()
    dev.new()
    dev.off()
    #old.par<-par()
    close.screen()
    split.screen(figs=c(5,4))
    new.par<-list(mar = c(1, 1, 1, 1),
                  mgp=c(1,0.1,0),
                  tcl=0.5,
                  ann=FALSE,
                  col.axis="black",
                  cex.axis=1,
                  xpd=FALSE,
                  cex.main=1,
                  cex=0.75)
    for (i in 1:length(resid)) {
      screen(i)
      par(new.par)
      temp.x.data<-data[[ random.plots$random.plots[[i]] ]][,pred[[k]]]
      if (is.character(temp.x.data)==TRUE) {
        temp.x.data<-as.factor(temp.x.data)
      }
      plot(resid[[ random.plots$random.plots[[i]] ]]$resid ~ temp.x.data,
           ylab=NA,
           xlab=NA)
      lines(lowess(resid[[ random.plots$random.plots[[i]] ]]$resid ~
                     temp.x.data), col="red")
      if (i==1) {
        title(main=paste("resid x", pred[[k]]))
      }
    }
    ask.reveal.data(real.data=random.plots$real.data,ask=ask)
    if (k<n.pred) {
      ask.pred<-utils::askYesNo(msg="Show next?")
      if (ask.pred==FALSE) {
        if (k<n.pred){k<-k+1} else {stop("There are no more plots to show")}
      }
      if (is.na(ask.pred)==TRUE){stop("Cancelled")}
    }
  }
}

#' Plot random residuals from `insp.sim()` mixed-model data.
#'
#'
#' @param resid list with model residuals
#' @param ask whether ask before showing real data plot
plots.qqnorm.random.r <- function(resid,
                               ask) {
  n.pred<-ncol(resid[[1]]$resid.random)
  for (k in 1:n.pred) {
    random.plots<-random.plots()
    dev.new()
    dev.off()
    #old.par<-par()
    close.screen()
    split.screen(figs=c(5,4))
    new.par<-list(mar = c(1, 1, 1, 1),
                  mgp=c(1,0.1,0),
                  tcl=0.5,
                  ann=FALSE,
                  col.axis="black",
                  cex.axis=1,
                  xpd=FALSE,
                  cex.main=1,
                  cex=0.75)
    for (i in 1:length(resid)) {
      screen(i)
      par(new.par)
      qqnorm(resid[[ random.plots$random.plots[[i]] ]]$resid.random[,k],
             main = NULL,
             ylab = NA,
             xlab = NA)
      qqline(resid[[ random.plots$random.plots[[i]] ]]$resid.random[,k], col="red")
      if (i==1) {
        title(main=paste(names(resid[[1]]$resid.random)[[k]],resid[[1]]$ylab.random))
      }
    }
    ask.reveal.data(real.data=random.plots$real.data,ask=ask)
    if (k<n.pred) {
      ask.pred<-utils::askYesNo(msg="Show next?")
      if (ask.pred==FALSE) {
        if (k<n.pred){k<-k+1} else {stop("There are no more plots to show")}
      }
      if (is.na(ask.pred)==TRUE){stop("Cancelled")}
    }
  }
}

#' Plot residuals auto-correlation estimates
#' from `insp.sim()` data.
#'
#' This function depends on the order of the data.
#'
#' @param models list with models from `insp.sim()` function.
#' @param ask whether ask before showing real data plot.
plots.acf.r<-function(models, ask) {
  random.plots<-random.plots()
  dev.new()
  dev.off()
  #old.par<-par()
  close.screen()
  split.screen(figs=c(5,4))
  new.par<-list(mar = c(1, 1, 1, 1),
                mgp=c(1,0.1,0),
                tcl=0.5,
                ann=FALSE,
                col.axis="black",
                cex.axis=1,
                xpd=FALSE,
                cex.main=1,
                cex=0.75)
  for (i in 1:length(models)) {
    screen(i)
    par(new.par)
    if(class(models[[1]])[[1]]=="gls"|class(models[[1]])[[1]]=="lme") {
      acf(resid(models[[ random.plots$random.plots[[i]] ]], type="normalized"),
          main=NULL)
    } else {
      acf(resid(models[[ random.plots$random.plots[[i]] ]]), main=NULL)
    }

    if (i==1) {
      if(class(models[[1]])[[1]]=="gls"|class(models[[1]])[[1]]=="lme") {
        title(main="ACF-Stand. Resididuals")
        } else {
          title(main="ACF-Ordinary Residuals")
        }
    }
  }
  ask.reveal.data(real.data=random.plots$real.data, ask=ask)
}

#' Plot residuals cook's distances from `insp.sim()` data.
#'
#' @param models list with models from `insp.sim()` function.
#' @param ask whether ask before showing real data plot.
plots.cook.r<-function(models, ask) {
  random.plots<-random.plots()
  dev.new()
  dev.off()
  #old.par<-par()
  close.screen()
  split.screen(figs=c(5,4))
  new.par<-list(mar = c(1, 1, 1, 1),
                mgp=c(1,0.1,0),
                tcl=0.5,
                ann=FALSE,
                col.axis="black",
                cex.axis=1,
                xpd=FALSE,
                cex.main=1,
                cex=0.75)
  for (i in 1:length(models)) {
    screen(i)
    par(new.par)
    plot(models[[ random.plots$random.plots[[i]] ]],
         caption=NA,
         main=NA,
         sub=NA,
         which=4)
    abline(h=1, col="red")
    abline(h=0.5, col="red", lty=2)
    if (i==1) {
      title(main="Cook's dist.")
    }
  }
  ask.reveal.data(real.data=random.plots$real.data, ask=ask)
}


#' Linear Model Validation Plots
#'
#' Linear Model Validation Plots for `lm()`, `glm()`, `lmer()`,
#' `lme()` and `gls()` models.
#'
#' @param model A fitted model. Accepts lm(), glm(), gls(), lme(), and lmer() models.
#' @param data Model data.
#' @param bin.size For glm binomial family only. Binned residuals.
#' @return Model specific residuals.
#' @references Zuur et al. 2009
#' @author Felipe M. Gawryszewski
insp.plot <- function (model, data=extract.data(model), bin.size=10) {

  dev.new()
  dev.off()

  data<-data.frame(data)
  ngraph<-6+ncol(data)
  if (class(model)[[1]]=="gls") {
    ngraph<-ngraph-1
  }
  if (class(model)[[1]]=="glm") {
    if (model$family[[1]]=="binomial") {
      ngraph<-ngraph+3
    }
  }
  if (class(model)[[1]] == "lme") {
    ngraph<-ngraph+1
  }

  ngraph1<-round(sqrt(ngraph),0)
  ngraph2<-ceiling(ngraph/ngraph1)

  par<-par()
  par(mfrow=c(ngraph1,ngraph2), mar = c(5, 4, 4, 2) + 0.1)

  if (class(model)[[1]] == "lm") {
    resid<-resid(model)
    ylab<-"Ordinary residuals"
  }
  if (class(model)[[1]] == "gls") {
    resid<-resid(model, type="normalized")
    ylab<-"Standardised residuals"
  }
  if (class(model)[[1]] == "lme") {
    resid<-resid(model, type='pearson', level=1)
    nranef<-length(nlme::ranef(model))
    resid.random<-nlme::ranef(model)
    ylab<-"Pearson resid. - Fixed Effects (L1)"
    ylab.random<-"Pearson resid.- Random Effects"
  }
  if (class(model)[[1]] == "lmerMod") {
    resid<-resid(model, type="pearson", level=1)
    ylab<-"Pearson resid. - Fixed Effects (L1)"
    nranef<-ncol(lme4::ranef(model)[[1]])
    resid.random<-lme4::ranef(model)[[1]]
    ylab.random<-"Pearson resid.- Random Effects"

  }
  if (class(model)[[1]] == "glm") {
    resid<-resid(model, type="deviance")
    ylab<-"Deviance Residuals"
    family<-model$family[[1]]

    if(family=="quasipoisson" | family=="quasibinomial") {
      dispersion<-summary(model)$dispersion
      mu <- predict(model, type = "response")
      E <- model$y - mu
      resid <- E / sqrt(dispersion * mu)
      ylab<-"Scaled Pearson Residuals"
    }
    if(family=="binomial") {
      resid<-resid(model, type="pearson")
      ylab<-"Pearson Residuals"
      ylab.bin<-"Pearson Residuals (bins)"

      bin.size<-bin.size
      bins<-round(seq(from=1,to=length(resid), length.out=length(resid)/bin.size))
      bin.factor<-vector(length=length(resid))
      resid.mean<-vector(length = length(bins)-1)
      #resid.sd<-resid
      if (length(bins)<=1) {stop("less than two bins. Decrease bin size.")}
        for (i in 1:(length(bins)-1)){
        bin.factor[bins[[i]]:bins[[i+1]]]<-rep(i, length(bins[[i]]:bins[[i+1]]))
        resid.mean[[i]]<-mean(resid[bins[[i]]:bins[[i+1]]])
        #resid.sd[[i]]<-sd(resid[bins[[i]]:bins[[i+1]]])
        }
      bin.factor<-as.factor(bin.factor)
      }
  }

  hist(resid, xlab=ylab, main=NULL)
  qqnorm(resid, main = "Main Effect Normal Q-Q Plot")
  qqline(resid, col="red")

  if (class(model)[[1]] == "glm") {
    if(family=="binomial") {
      boxplot(resid~bin.factor, ylab=ylab.bin, xlab="Bins")
      lines(lowess(resid ~ bin.factor), col="red")
      hist(resid.mean, xlab=ylab.bin, main=NULL)
      qqnorm(resid.mean, main = "Binned Residuals Q-Q Plot")
      qqline(resid.mean, col="red")
    }
  }

  if(class(model)[[1]]=='lme'){
    if(nranef==1) {
      qqnorm(resid.random[[1]], main = paste("Random Effect Q-Q Plot"))
      qqline(resid.random[[1]], col="red")
    }
    if(nranef>1) {
      for (i in 1:nranef) {
        qqnorm(unlist(resid.random[[i]]), main = paste("Random Effect Q-Q Plot", i))
        qqline(unlist(resid.random[[i]]), col="red")
      }
    }
  }
  if(class(model)[[1]]=="lmerMod"){
    if(nranef==1) {
      qqnorm(resid.random[[1]], main = paste("Random Effect Q-Q Plot"))
      qqline(resid.random[[1]], col="red")
    }
    if(nranef>1) {
      for (i in 1:nranef) {
        qqnorm(unlist(resid.random[[i]]), main = paste("Random Effect Q-Q Plot", i))
        qqline(unlist(resid.random[[i]]), col="red")
      }
    }
  }

  if(class(model)[[1]]=='lme') {

  plot(resid ~ fitted(model, level=0), ylab=ylab, xlab="Fitted")
  lines(lowess(resid ~ fitted(model, level=0)), col="red")

  plot(resid(model, type="pearson", level=0) ~ fitted(model, level=0), ylab="Stand. residuals - Fixed Effects - L0", xlab="Fitted")
  lines(lowess(resid(model, type="pearson", level=0) ~ fitted(model, level=0)), col="red")


  } else {

  plot(resid ~ fitted(model), ylab=ylab, xlab="Fitted")
  lines(lowess(resid ~ fitted(model)), col="red")

  }
  for (i in 1:ncol(data)) {
    plot(resid ~ data[,i], ylab=ylab, xlab=colnames(data)[[i]])
    lines(lowess(resid ~ data[,i]), col="red")
  }
  if(class(model)[[1]]=="lm" | class(model)[[1]]=="glm") {
    #plot(model, which=3)
    plot(model, which=4)
    abline(h=1, col="red")
    abline(h=0.5, col="red", lty=2)
  }


  if(class(model)[[1]]=="gls") {
    acf(resid(model, type="normalized"), main="Stand. Resid.")
  }

  if(class(model)[[1]]!="gls") {
    acf(resid(model), main="Ord. Resid.")
  }

  if(class(model)[[1]]=="glm") {
    print(paste("Residual Deviance/DF =", round(model$deviance/model$df.residual,3)))
  }
  par(mfrow=par$mfrow, mar = par$mar)
}
