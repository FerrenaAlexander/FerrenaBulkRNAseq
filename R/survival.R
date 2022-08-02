


# to do items
# finish adding saving parts of the functions to the end of this
# may need to split saving for continuous vs categorical

# also, in deplots.r, make sure to check heatmapplot:
# added suport for sig vs nonsig genes --> test it
# need to document some of the new params in that function


#' A function for survival analysis.
#'
#' An easy wrapper around survival analysis as implemented in the packages survival and survminer.
#'
#' @param testvardf - a data.frame of variables to test for association with survival, columns = variables, rows = observations. If you want to run a gene expression matrix, then genes need to be columns and samples need to be rows, ie you may need to transpose via the transpose function, t().
#' @param clinvardf - a data.frame with clinical variables, including time and status columns, columns = variables, rows = observations.
#' @param vartypes - a character vector of length = ncol(testvardf). Should say either "continuous" or "categorical", for each column of testvardf. If empty, will guess the vartype, continuous has >12 unique values.
#' @param multivarnames - a character vector. The colnames of clinvardf to include in multivariable modelling, if desired. If empty, does not perform multivar modelling.
#' @param timevarname - a string, the colname of the column in clinvardf with time to event information.
#' @param statusvarname - a string, the colname of the column in clinvardf with event information. The event column should be a character vector and should have the word "Censored" (capital C); ie alive/dead should be "Censored" and "Dead", event should be "Censored" and "Event".
#' @param outdir - a string, the path to write outputs to; defaults to './survival'
#'
#' @return A nested list object: highest level, categorical vs continuous; below that, lists containing plots, models, and data; plotlists will contain kaplan-meier plots; models will include cox regresion models; datalists contain the data DFs, which is important for reproducibility and model diagnostics including proportional hazards testing.
#' @export
#'
#' @examples
survival <- function(testvardf,
                     clinvardf,
                     vartypes,
                     multivarnames,
                     timevarname,
                     statusvarname,
                     outdir){




  message("\nInitiating analysis")

  if( missing(outdir) ) {outdir <- './survival'; message('\nNo outdir provided, defaulting to write outs to "./survival"')}

  ### try to get continuous vs cat test vars ###
  # if vartpyes is provided (not missing), use that
  if( !missing(vartypes) ) {

    message("\nUsing vartypes to separate continuous and categorical variables...")

    if(length(vartypes) != ncol(testvardf)){error("vartypes is not the same length as ncol(testvardf). make sure each column is accounter for in vartypes as either 'continuous' or 'categorical', or leave vartypes blank.")}

    # using the vartypes, get the continuous / categorical columns
    contvars <- testvardf[,vartypes=='continuous', drop=FALSE]
    catvars <- testvardf[,vartypes=='categorical', drop=FALSE]

    #get their names...
    contnames <- colnames(contvars)
    catnames <- colnames(catvars)

    message(' - ', ncol(contvars), ' continuous variables')
    message(' - ', ncol(catvars), ' categorical variables')

  }

  # if vartypes is missing, then try to guess based on continuous =  >12 unique values
  if( missing(vartypes) ) {

    message('\nGuessing vartypes based on continuous defined as > 12 unique values...')

    #get num unique values for each variable
    num_unique_vals <- sapply(testvardf, FUN = function(x){length(unique(x))}, simplify = T)

    #categorical
    catnames <- names(num_unique_vals[num_unique_vals<=12])
    catvars <- testvardf[,catnames, drop=FALSE]

    #continuous
    contnames <- names(num_unique_vals[num_unique_vals>12])
    contvars <- testvardf[,contnames, drop=FALSE]

    message(' - ', ncol(contvars), ' continuous variables')
    message(' - ', ncol(catvars), ' categorical variables')

  }








  #### continuous vars ####


  if( length(contnames) > 0){

    message("\nBeginning analysis for continuous variables...")

    #get varnames
    scorenames <- names(contvars)

    # prep progress bar
    total = ncol(contvars)
    pb <- txtProgressBar(min = 0, max = total, style = 3)

    #outs from loop: plots, models, data
    plotlist <- list()
    modellist <- list()
    datalist <- list()

    # this list will wrap the cat var lists, if there are cat vars; else, NULL.
    contvarout <- list()


    for(scoreidx in 1:length(scorenames) ) {


      #get the score name
      scorename <- colnames(contvars)[scoreidx]

      #get the score valeues, numeric continuopus
      scorevec <- contvars[,scorename]

      ## dichotomize the score ##
      highname <- paste0('High\nN=', table(na.omit(scorevec) >= median(na.omit(scorevec)) )[2] )
      lowname <- paste0('Low\nN=', table(na.omit(scorevec) >= median(na.omit(scorevec)) )[1] )

      dichscore <- ifelse(scorevec >= median(na.omit(scorevec)), yes = highname, no = lowname)
      dichscore <- factor(dichscore, levels = c(lowname, highname))

      ### split score to tertiles ###
      terciles <- factor(ntile(scorevec, 3))
      terciles <- plyr::mapvalues(x = terciles, from = levels(terciles),
                                  to = c( paste0('Low\nN=', table(terciles)[1]),
                                          paste0('Med\nN=', table(terciles)[2]),
                                          paste0('High\nN=', table(terciles)[2]) )
      )



      ### set up survival datadf ###
      datadf <- data.frame(continuous = scorevec,
                           dichotomized = dichscore,
                           terciles = terciles,
                           surv = clinvardf[,timevarname],
                           status = clinvardf[,statusvarname])




      # #check dist
      # #check dist over dead/alive
      # hist <- ggplot(datadf, aes(continuous, fill=status))+
      #   geom_histogram(col='black',position='identity', alpha=0.5) +
      #   facet_wrap(~status, nrow=2)
      #
      # #check dist cor with survival?
      #
      # corCensored <- cor.test(datadf[datadf$status=='Censored', "surv"], datadf[datadf$status=='Censored', "continuous"])
      # corNonCensored <- cor.test(datadf[datadf$status!='Censored', "surv"], datadf[datadf$status!='Censored', "continuous"])
      # corlab <- paste0('Non-Censored cor: ', round(corNonCensored$estimate, 3), ', Non-Censored P:', round(corNonCensored$p.value, 3),
      #                  '\nCensored cor: ', round(corCensored$estimate, 3), ', Censored P: ', round(corCensored$p.value, 3))
      # corp <- ggplot(datadf, aes(continuous, surv, col = status))+
      #   geom_point()+
      #   geom_smooth()+
      #   facet_wrap(~status, nrow=2)+
      #   labs(caption = corlab)
      #
      # distplots <- hist+corp

      ### code status as numeric ###
      # EVENT = 1, CENSORED = 0
      datadf$status_code <- ifelse(datadf$status == 'Censored', yes=0, no = 1)



      ### tertile analysis ###
      fit <- surv_fit(Surv(surv, status_code) ~ terciles, data = datadf)

      gterc <- ggsurvplot(fit,
                          conf.int = F,
                          pval = TRUE,
                          pval.method = T,
                          linetype = "strata",
                          test.for.trend = F,
                          title=scorename,
                          ggtheme = theme_classic2(),
                          palette =  c( "#8DA0CB", "#66C2A5", "#FC8D62" ) )


      ### dichotomized analysis ###
      fit <- surv_fit(Surv(surv, status_code) ~ dichotomized, data = datadf)

      gdich <- ggsurvplot(fit,
                          conf.int = F,
                          pval = TRUE, #pval.coord = c(0, 0.15),
                          pval.method = T, #pval.method.coord = c(0, 0.25),
                          linetype = "strata",
                          test.for.trend = F,
                          title=scorename,
                          ggtheme = theme_classic2(),
                          palette = c("#377EB8" , "#E41A1C") )

      # just get plots
      gterc <- gterc$plot
      gdich <- gdich$plot

      ### cox modelling ###
      #dich model, will match plot
      dich <- coxph(Surv(surv, status_code) ~ dichotomized, data = datadf)

      #terc
      terc <- coxph(Surv(surv, status_code) ~ terciles, data = datadf)

      #continuous univariate cox
      continuous <- coxph(Surv(surv, status_code) ~ continuous, data = datadf)



      #multivar cox test
      if( !missing(multivarnames) ){

        #prep for multivar
        datadf <- cbind(datadf, clinvardf)

        # run it
        multivar <- coxph(as.formula(
          paste0('Surv(surv, status_code) ~ ',
                 paste(multivarnames, collapse = ' + '))
        ),
        data=datadf)

      } else{multivar <- NULL}



      #bind all outputs
      # models
      models <- list(dich = dich,
                     terc = terc,
                     continuous = continuous,
                     multivar = multivar
      )

      # plots
      # plots <- list(gdich, gterc, distplots)
      plots <- list(gdich, gterc)

      #save outs
      plotlist[[scorename]] <- plots
      modellist[[scorename]] <- models
      datalist[[scorename]] <- datadf


      setTxtProgressBar(pb, scoreidx)

    } # end cont var loop.


    contvarout <- list(plotlist = plotlist,
                       modellist = modellist,
                       datalist = datalist)


  } else {contvarout <- NULL} # end cont var if statement.












  #### categorical vars ####

  if( length(catnames) > 0 ){

    message("\nBeginning analysis for categorical variables...")


    #get varnames
    scorenames <- names(catvars)

    # prep progress bar
    total = ncol(catvars)
    pb <- txtProgressBar(min = 0, max = total, style = 3)


    #prep output lists

    #outs from loop: plots, models, data
    plotlist <- list()
    modellist <- list()
    datalist <- list()

    # this list will wrap the cat var lists, if there are cat vars; else, NULL.
    catvarout <- list()

    for(scoreidx in 1:length(scorenames)) {

      #get the score name
      scorename <- colnames(catvars)[scoreidx]

      #get the score valeues, numeric continuopus
      scorevec <- catvars[,scorename]

      # change name to have N
      scorevec <- plyr::mapvalues(scorevec, from = levels(scorevec),
                                  to = paste0(levels(scorevec), '\nN=', table(scorevec)[levels(scorevec)] ) )



      datadf <- data.frame(var = scorevec,
                           surv = clinvardf[,timevarname],
                           status = clinvardf[,statusvarname])


      ### code status as numeric ###
      # EVENT = 1, CENSORED = 0
      datadf$status_code <- ifelse(datadf$status == 'Censored', yes=0, no = 1)



      fit <- surv_fit(Surv(surv, status_code) ~ var, data = datadf)

      gcat <- ggsurvplot(fit,
                         conf.int = F,
                         pval = TRUE,
                         pval.method = T,
                         linetype = "strata",
                         test.for.trend = F,
                         title=scorename,
                         ggtheme = theme_classic2())

      #test for trend, but only if >2 vars...
      if( length(levels(datadf$var)) > 2 )  {



        gtft <- ggsurvplot(fit,
                           conf.int = F,
                           pval = TRUE,
                           pval.method = T,
                           linetype = "strata",
                           test.for.trend = T,
                           title=scorename,
                           ggtheme = theme_classic2())

        gtft <- gtft$plot

      }

      gcat <- gcat$plot

      if( length(levels(datadf$var)) > 2 )  {
        plots <- list(gcat, gtft)

      } else{
        plots <- gcat
      }

      ### cox modelling ###
      #dich model, will match plot
      model <- coxph(Surv(surv, status_code) ~ var, data = datadf)


      #multivar cox test
      if( !missing(multivarnames) ){

        #prep for multivar
        datadf <- cbind(datadf, clinvardf)

        # run it
        multivar <- coxph(as.formula(
          paste0('Surv(surv, status_code) ~ ',
                 paste(multivarnames, collapse = ' + '))
        ),
        data=datadf)

      } else{multivar <- NULL}



      #bind all outputs
      # models
      models <- list(model = model,
                     multivar = multivar
      )

      # plots
      # plots <- list(gdich, gterc, distplots)
      plots <- plots

      #save outs
      plotlist[[scorename]] <- plots
      modellist[[scorename]] <- models
      datalist[[scorename]] <- datadf


      setTxtProgressBar(pb, scoreidx)

    } # end cat var loop.


    catvarout <- list(plotlist = plotlist,
                      modellist = modellist,
                      datalist = datalist)

  } else{catvarout <- NULL} #end cat var if statement.












  ##### WRITE OUTS, SAVE PLOTS ETC #####

  outlist <- list(continuous = contvarout,
                  categorical = catvarout)

  message('\nPreparing to save outputs to: "', outdir,'"')
  dir.create(outdir)



  ### save continuous vars ###

  if(is.null( outlist[[1]]) ){


    message('\nSaving continuous variable outputs...')

    #get outputs
    plotlist <- outlist[[1]][[1]]
    modellist <- outlist[[1]][[2]]
    datalist <- outlist[[1]][[3]]



    sigpways <- names(plotlist)

    sigresdir <- paste0(outdir, '/continuousvars/')
    dir.create(sigresdir)

    total = length(sigpways)
    pb <- txtProgressBar(min = 0, max = total, style = 3)


    #for significant results, make outputs
    pdfplotlist <- list()
    for(scoreidx in 1:length(sigpways) ){

      score <- sigpways[scoreidx]
      #message(score)

      #get models
      models <- modellist[[score]]
      datadf <- datalist[[score]] #load data into env or else, the diagnostics break

      #get plots
      plots <- plotlist[[score]]


      #get the coefficients, hazard ratios, confidence intervals, and P values
      summarytable <- bind_rows( lapply(models, function(x){as.data.frame( summary(x)$coefficients ) }) )
      summarytable <- bind_cols(summarytable,
                                bind_rows( lapply(models, function(x){as.data.frame( summary(x)$conf.int ) }) )[,3:4] )





      #only plot if significant...
      # ignore met
      # pvals <- summarytable$`Pr(>|z|)`[-(nrow(summarytable))]
      # if( !any(pvals < 0.10) ){
      #   next
      # }
      # this is defined by sigpways, not necessary

      #rename the continuous
      rownames <- rownames(summarytable)
      rownames[grepl('^continuous', rownames)] <- c('Univariable\nContinuous', 'Multivariable\nContinuous')
      rownames(summarytable) <- rownames

      #put conf int next to coef
      summarytable <- summarytable[,c(1,6,7,2,3,4,5)]

      #hr, lower/upperCI, P, and sig only...
      summarytable <- summarytable[,c(4,2,3,7)]


      #add significance marks
      summarytable$significance <- ''
      summarytable[summarytable$`Pr(>|z|)` < 0.1, "significance"] <- '.'
      summarytable[summarytable$`Pr(>|z|)` < 0.05, "significance"] <- '*'
      summarytable[summarytable$`Pr(>|z|)` < 0.01, "significance"] <- '**'
      summarytable[summarytable$`Pr(>|z|)` < 0.001, "significance"] <- '***'

      #round table
      summarytable[,-ncol(summarytable)] <- round( summarytable[,-ncol(summarytable)],3 )

      # plot it all together
      updatedplot <- (plots[[1]] + plots[[2]]) / ( gridExtra::tableGrob(summarytable)) + patchwork::plot_annotation(score) +
        plot_layout(heights = c(1,2))


      pdfplotlist[[score]] <- updatedplot



      #save all outputs...
      resultdir <- paste0(sigresdir, '/', score)
      dir.create(resultdir)


      # saveRDS( file = paste0(resultdir, '/models.rds'), models  )


      #check model assumptions...
      modeldiagnosticsfile <- paste0(resultdir, '/modeldiagnostics.pdf')
      pdf(modeldiagnosticsfile, height = 7, width = 7)
      for(model in models){

        print( ggcoxzph( cox.zph(model) ) )

      }
      dev.off()


      #save plots


      gdich <- plots[[1]]
      gterc <- plots[[2]]
      distplots <- plots[[3]]
      suppressMessages(ggsave( paste0(resultdir, '/dist.jpg'), distplots , height = 7, width = 7))

      ggsave( paste0(resultdir, '/gdich.jpg'), gdich , height = 5, width = 5)
      ggsave( paste0(resultdir, '/gterc.jpg'), gterc , height = 5, width = 5)


      pdf(paste0(resultdir, '/summaryplot.pdf'), height = 10, width = 9)
      updatedplot
      dev.off()

      summarytable <- cbind(rownames(summarytable), summarytable)

      write.csv(paste0(resultdir, '/modeltable.csv'), x=summarytable, row.names = F)


      ### saveas much as possible just in case ###
      list_to_save <- list(datadf = datadf, plots = plots, models = models)
      saveRDS(object = list_to_save, file = paste0(resultdir, '/objects.rds') )

      setTxtProgressBar(pb, scoreidx)


  }




  }




  return(outlist)

}
