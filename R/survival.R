


# to do items
# customizable forest plots
# add in a limit on how many plots to save in summary: if running thousands, only top 30 significant...?

# also, in deplots.r, make sure to check heatmapplot:
# added suport for sig vs nonsig genes --> test it
# need to document some of the new params in that function


#' A function for survival analysis.
#'
#' An easy wrapper around survival analysis as implemented in the packages `survival` and `survminer`.
#' This package is meant to easily facilitate survival analysis for test variables (given in `testvardf`), using time-to-event and event-status variables (given in `clinvardf`).
#' It also supports multivariable modelling using other variables in `clinvardf`, which are called according to their column names as provided by `multivarnames`.
#'
#' @param testvardf - a data.frame of variables to test for association with survival, columns = variables, rows = observations. If you want to run a gene expression matrix, then genes need to be columns and samples need to be rows, ie you may need to transpose via the transpose function, t().
#' @param clinvardf - a data.frame with clinical variables, including time and status columns, columns = variables, rows = observations.
#' @param vartypes - a character vector of length = ncol(testvardf). Should say either "continuous" or "categorical", for each column of testvardf. If empty, will guess the vartype, defining continuous variables as those with >12 unique values.
#' @param multivarnames - a character vector. The colnames of clinvardf to include in multivariable modelling, if desired. If empty, does not perform multivar modelling.
#' @param timevarname - a string, the colname of the column in clinvardf with time to event information.
#' @param statusvarname - a string, the colname of the column in clinvardf with event information. The event column should be a character vector and should have the word "Censored" (capital C); ie alive/dead should be "Censored" and "Dead", event should be "Censored" and "Event".
#' @param outdir - a string, the path to write outputs to; defaults to './survival'
#'
#' @return Does not return anything; instead, will save survival analysis plots, Cox regression models, and frozen data to outdir.
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
          paste0('Surv(surv, status_code) ~ continuous + ',
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

      #drop null list; this happens if no multivar was used.
      models <- models[lengths(models) != 0]

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
        plots <- list( gcat, plot_spacer() )

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
          paste0('Surv(surv, status_code) ~ var + ',
                 paste(multivarnames, collapse = ' + '))
        ),
        data=datadf)

      } else{multivar <- NULL}



      #bind all outputs
      # models
      models <- list(model = model,
                     multivar = multivar
      )

      #drop null list; this happens if no multivar was used.
      models <- models[lengths(models) != 0]

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

  message('\n\n\nPreparing to save outputs to: "', outdir,'"')
  dir.create(outdir, recursive = T)



  ### save continuous vars ###

  if( ! is.null( outlist[[1]]) ){


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




      #rename the continuous
      # make sure to
      rownames <- rownames(summarytable)

      rownames[1:4] <- paste0('Univariate_', rownames[1:4])
      rownames[4] <- 'Univariate_Continuous'

      if(length(rownames)>4){
        rownames[5:length(rownames)] <- paste0('Multivariate_', rownames[5:length(rownames)])
        rownames[5] <- 'Multivariate_Continuous'
      }

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
      updatedplot <- (plots[[1]] + plots[[2]]) / ( gridExtra::tableGrob(summarytable, theme=gridExtra::ttheme_default(base_size =8))) + patchwork::plot_annotation(score) +
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

      pdf(paste0(resultdir, '/KM_dichotomized.pdf'), height = 5, width = 5)
      print( gdich )
      dev.off()

      pdf(paste0(resultdir, '/KM_terciles.pdf'), height = 5, width = 5)
      print( gterc )
      dev.off()

      pdf(paste0(resultdir, '/summaryplot.pdf'), height = 10, width = 9)
      print( updatedplot )
      dev.off()

      summarytable <- cbind(rownames(summarytable), summarytable)
      colnames(summarytable)[1] <- 'variable'

      write.csv(paste0(resultdir, '/modeltable.csv'), x=summarytable, row.names = F)


      ### saveas much as possible just in case ###
      list_to_save <- list(datadf = datadf, plots = plots, models = models)
      saveRDS(object = list_to_save, file = paste0(resultdir, '/objects.rds') )

      setTxtProgressBar(pb, scoreidx)


    } # end cont saving loop


    pdfplotlist_cont <- pdfplotlist

  } else{pdfplotlist_cont <- NULL} # end cont saving if statement


  ### save categorical vars ###

  if( ! is.null( outlist[[2]]) ){


    message('\nSaving categorical variable outputs...')

    #get outputs
    plotlist <- outlist[[2]][[1]]
    modellist <- outlist[[2]][[2]]
    datalist <- outlist[[2]][[3]]



    sigpways <- names(plotlist)

    sigresdir <- paste0(outdir, '/categoricalvars/')
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






      #rename the univar
      # get the rows form then number of levels - 1....
      numtestlevs <- length(levels(datadf$var)) - 1
      rownames(summarytable)[1:numtestlevs] <- paste0('Univariate_', rownames(summarytable)[1:numtestlevs] )


      #rename the multivar
      # check if multivar even exists: nrow will be greater than numtestlevs
      # to get the rest of the rows, use numtestlevs + 1 --> all else is multivar...
      if(nrow(summarytable) > numtestlevs ){
        numtestlevs <- numtestlevs + 1
        rownames(summarytable)[numtestlevs:nrow(summarytable)] <- paste0('Multivariate_', rownames(summarytable)[numtestlevs:nrow(summarytable)] )
      }

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
      updatedplot <- (plots[[1]] + plots[[2]]) / ( gridExtra::tableGrob(summarytable, theme=gridExtra::ttheme_default(base_size =8))) + patchwork::plot_annotation(score) +
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


      gcat <- plots[[1]]
      gtft <- plots[[2]]

      pdf(paste0(resultdir, '/KM_categorical.pdf'), height = 5, width = 5)
      print( gcat )
      dev.off()

      pdf(paste0(resultdir, '/KM_categorical_test-for-trend.pdf'), height = 5, width = 5)
      print( gtft )
      dev.off()

      pdf(paste0(resultdir, '/summaryplot.pdf'), height = 10, width = 9)
      print( updatedplot )
      dev.off()

      summarytable <- cbind(rownames(summarytable), summarytable)
      colnames(summarytable)[1] <- 'variable'

      write.csv(paste0(resultdir, '/modeltable.csv'), x=summarytable, row.names = F)


      ### saveas much as possible just in case ###
      list_to_save <- list(datadf = datadf, plots = plots, models = models)
      saveRDS(object = list_to_save, file = paste0(resultdir, '/objects.rds') )

      setTxtProgressBar(pb, scoreidx)


    } # end cat saving loop

    pdfplotlist_cat <- pdfplotlist

  } else{pdfplotlist_cat <- NULL} # end cat saving if statement



  # summaryfiles:
  # forestplots;
  # summaryplots and tables

  summaryplotdir <- paste0(outdir, '/summaryplots')

  message('\n\n\nSaving summary plots to: "', summaryplotdir,'"')
  dir.create(summaryplotdir, recursive = T)


  if( ! is.null( outlist[[1]]) ){

    #forestplots for dich; cont; and multivar

    #get outputs
    plotlist <- outlist[[1]][[1]]
    modellist <- outlist[[1]][[2]]
    datalist <- outlist[[1]][[3]]



    forest_univar <- internal_forestplot_continuous_univar(modellist)
    forest_dich <- internal_forestplot_continuous_dich(modellist)


    if(length(modellist[[1]])==4){
      forest_multivar <- internal_forestplot_continuous_multivar(modellist)
    }



    pdf(  paste0(summaryplotdir, '/summary_continuousvars.pdf')  )
    print( forest_univar )
    print( forest_dich )
    if(length(modellist[[1]])==4){ print( forest_multivar ) }
    print( pdfplotlist_cont )

    suppressMessages(dev.off())


  }

  if( ! is.null( outlist[[2]]) ){

    #get outputs
    plotlist <- outlist[[2]][[1]]
    modellist <- outlist[[2]][[2]]
    datalist <- outlist[[2]][[3]]


    forest_cat_univar <- internal_forestplot_categorical_univar(modellist)

    if(length(modellist[[1]])==2){
      forest_cat_multivar <- internal_forestplot_categorical_multivar(modellist)
    }


    pdf(  paste0(summaryplotdir, '/summary_categoricalvar.pdf')  )
    print( forest_cat_univar )
    if(length(modellist[[1]])==2){ print ( forest_cat_multivar ) }
    print( pdfplotlist_cat )

    suppressMessages(dev.off())

  }

  message('\nAll done!')

}















##### internal functions #####

### fix long names function
#' Fix long strongs.
#'
#' Long strings may contain groups of words separated by underscores or other characters, such as pathway names. This function will try to cut the string in half, replacing the delimiter character closest to the middle of the string with a new-line. Good for plottng things.
#'
#' @param pnames character vector containing strings with long names. If longer than `threshold` characters, will split in half.
#' @param delimiter string. the character separating words in the strings in pnames. default = '_'
#' @param threshold integer. strings of length > threshold are cut in half with a new-line at the `delimiter` closest to the center of the string.
#'
#' @return
#' @export
#'
#' @examples
fixlongnames <- function(pnames, delimiter, threshold){

  if(missing(delimiter)){delimiter <- '_'}
  if(missing(threshold)){threshold <- 50}

  pnewnames <- c()
  for(p in pnames){


    if(str_length(p) > threshold){

      #try to find and replace underscore closest to halfway...
      halfway <- round(str_length(p)/2)
      underscorepos <- str_locate_all(p, delimiter)[[1]][,1]

      distances <- abs(halfway-underscorepos)

      which_und_is_closest <- which.min(distances)

      split <- str_split_fixed(p, delimiter,Inf)

      newp <- paste0(paste(split[1,1:which_und_is_closest], collapse = delimiter),
                     '        ', '\n',
                     paste(split[1,(which_und_is_closest+1):ncol(split)], collapse = delimiter)
      )



    } else{newp <- p}

    pnewnames <- c(pnewnames, newp)

  }

  pnewnames
}





internal_forestplot_continuous_univar <- function(modellist){
  ### check the hazard ratios ###


  #univar first
  univarlist <- list()
  for(score in names(modellist) ){
    summx <- summary(modellist[[score]]$continuous)
    univar <- as.data.frame(summx$conf.int)
    univar$p <- summx$coefficients[,5]
    colnames(univar) <- c('HR', 'hrneg', 'lower95', 'upper95', 'p')
    univar<- data.frame(score=score,univar, row.names = NULL)
    univarlist[[score]]  <- univar
  }
  resdf <- dplyr::bind_rows(univarlist)

  resdf$score_rename <- fixlongnames(resdf$score)


  #if more than 30, keep only 30 most significant on each side
  # no filtering for genes...
  resdf <- resdf[order(resdf$p),]

  #sort by HR
  resdf <- resdf[order(resdf$HR),]

  #significance...
  resdf$test <- ifelse(resdf$p < 0.05, 'significant', 'non-sig')
  resdf$test <- factor(resdf$test, levels = c('non-sig', 'significant'))

  # color...
  resdf$testcolor <- ifelse(resdf$test == 'significant', yes="#66C2A5", no = "#FC8D62")

  #factorize..
  resdf$score_rename <- factor(resdf$score_rename, levels=resdf$score_rename)

  forestunivar <- ggplot(resdf, aes(y = score_rename, x = HR, xmin=lower95, xmax=upper95, col = test))+
    geom_pointrange()+
    geom_vline(xintercept = 1, linetype = 'dashed')+
    theme_linedraw()+
    theme(axis.text.y = element_text(size = 7),
          axis.title.y = element_blank() )+
    # scale_color_manual(values=resdf$testcolor)+
    labs(title = 'HR Forest plot, univariate continuous')+
    xlab('Hazard Ratio\nLeft = better survival ; Right = worse survival')

  forestunivar

}


internal_forestplot_continuous_dich <- function(modellist){

  #hazard ratios for dich

  ### check the hazard ratios ###
  dichlist <- list()
  for(score in names(modellist) ){
    summx <- summary(modellist[[score]]$dich)
    dich <- as.data.frame(summx$conf.int)
    dich$p <- summx$coefficients[,5]
    colnames(dich) <- c('HR', 'hrneg', 'lower95', 'upper95', 'p')
    dich<- data.frame(score=score,dich, row.names = NULL)
    dichlist[[score]]  <- dich
  }
  resdf <- dplyr::bind_rows(dichlist)

  resdf$score_rename <- fixlongnames(resdf$score)


  #if more than 30, keep only 30 most significant on each side


  # if(nrow(resdf>30)){
  #   resdf <- resdf[order(resdf$p),]
  #   up <- resdf[resdf$HR>1,] ; down <- resdf[resdf$HR<1,]
  #   up <- up[up$p < 0.05,] ; down <- down[down$p < 0.05,]
  #
  #   resdf <- rbind(up,down)
  # }


  #sort by HR
  resdf <- resdf[order(resdf$HR),]

  #if more than 30, keep only 30 strongest...
  # if(nrow(resdf) > 30){
  #   resdf <- rbind( head(resdf, 15), tail(resdf, 15))
  # }


  #significance...
  resdf$test <- ifelse(resdf$p < 0.05, 'significant', 'non-sig')
  resdf$test <- factor(resdf$test, levels = c('non-sig', 'significant'))

  # color...
  resdf$testcolor <- ifelse(resdf$test == 'significant', yes="#66C2A5", no = "#FC8D62")

  #factorize..
  resdf$score_rename <- factor(resdf$score_rename, levels=resdf$score_rename)

  forestdich <- ggplot(resdf, aes(y = score_rename, x = HR, xmin=lower95, xmax=upper95, col=test))+
    geom_pointrange()+
    geom_vline(xintercept = 1, linetype = 'dashed')+
    theme_linedraw()+
    theme(axis.text.y = element_text(size = 7),
          axis.title.y = element_blank() )+
    # scale_color_manual(values=resdf$testcolor)+
    labs(title = 'HR Forest plot, dichotomized')+
    xlab('Hazard Ratio\nLeft = better survival ; Right = worse survival')

  forestdich

}


internal_forestplot_continuous_multivar <- function(modellist){

  #hazard ratios for res

  ### check the hazard ratios ###
  reslist <- list()
  for(score in names(modellist) ){
    summx <- summary(modellist[[score]]$multivar)
    res <- as.data.frame(summx$conf.int)['continuous',]
    res$p <- summx$coefficients[,5]['continuous']
    colnames(res) <- c('HR', 'hrneg', 'lower95', 'upper95', 'p')
    res<- data.frame(score=score,res, row.names = NULL)
    reslist[[score]]  <- res
  }
  resdf <- dplyr::bind_rows(reslist)

  resdf$score_rename <- fixlongnames(resdf$score)


  #if more than 30, keep only 30 most significant on each side


  # if(nrow(resdf>30)){
  #   resdf <- resdf[order(resdf$p),]
  #   up <- resdf[resdf$HR>1,] ; down <- resdf[resdf$HR<1,]
  #   up <- up[up$p < 0.05,] ; down <- down[down$p < 0.05,]
  #
  #   resdf <- rbind(up,down)
  # }


  #sort by HR
  resdf <- resdf[order(resdf$HR),]

  #if more than 30, keep only 30 strongest...
  # if(nrow(resdf) > 30){
  #   resdf <- rbind( head(resdf, 15), tail(resdf, 15))
  # }


  #significance...
  resdf$test <- ifelse(resdf$p < 0.05, 'significant', 'non-sig')
  resdf$test <- factor(resdf$test, levels = c('non-sig', 'significant'))

  # color...
  resdf$testcolor <- ifelse(resdf$test == 'significant', yes="#66C2A5", no = "#FC8D62")

  #factorize..
  resdf$score_rename <- factor(resdf$score_rename, levels=resdf$score_rename)

  forestres <- ggplot(resdf, aes(y = score_rename, x = HR, xmin=lower95, xmax=upper95, col=test))+
    geom_pointrange()+
    geom_vline(xintercept = 1, linetype = 'dashed')+
    theme_linedraw()+
    theme(axis.text.y = element_text(size = 7),
          axis.title.y = element_blank() )+
    # scale_color_manual(values=resdf$testcolor)+
    labs(title = 'HR Forest plot, multivariable continuous')+
    xlab('Hazard Ratio\nLeft = better survival ; Right = worse survival')

  forestres

}






internal_forestplot_categorical_univar <- function(modellist){

  #hazard ratios for res

  ### check the hazard ratios ###
  reslist <- list()
  for(score in names(modellist) ){
    summx <- summary(modellist[[score]]$model)
    res <- as.data.frame(summx$conf.int)
    res$p <- summx$coefficients[,5]

    colnames(res) <- c('HR', 'hrneg', 'lower95', 'upper95', 'p')
    res<- data.frame(score=score,level=rownames(res), res, row.names = NULL)
    reslist[[score]]  <- res
  }
  resdf <- dplyr::bind_rows(reslist)

  resdf$score_rename <- fixlongnames(resdf$score)


  #significance...
  resdf$test <- ifelse(resdf$p < 0.05, 'significant', 'non-sig')
  resdf$test <- factor(resdf$test, levels = c('non-sig', 'significant'))

  # color...
  resdf$testcolor <- ifelse(resdf$test == 'significant', yes="#66C2A5", no = "#FC8D62")

  #factorize..
  # resdf$score_rename <- factor(resdf$score_rename, levels=resdf$score_rename)

  forestres <- ggplot(resdf, aes(y = level, x = HR, xmin=lower95, xmax=upper95, col=test))+
    geom_pointrange()+
    geom_vline(xintercept = 1, linetype = 'dashed')+
    facet_wrap(~score_rename, ncol=1, scales='free_y')+
    theme_linedraw()+
    theme(axis.text.y = element_text(size = 7),
          axis.title.y = element_blank() )+
    # scale_color_manual(values=resdf$testcolor)+
    labs(title = 'HR Forest plot, categorical')+
    xlab('Hazard Ratio\nLeft = better survival ; Right = worse survival')


  forestres

}


internal_forestplot_categorical_multivar <- function(modellist){

  #hazard ratios for res

  ### check the hazard ratios ###
  reslist <- list()
  for(score in names(modellist) ){
    summx <- summary(modellist[[score]]$multivar)
    res <- as.data.frame(summx$conf.int)
    res$p <- summx$coefficients[,5]

    # make sure to select just the multivar ones of interest...
    # it will start with var, and alos will have a \n in it...
    res <- res[grepl('^var', rownames(res)),]
    res <- res[grepl('\nN', rownames(res)),]

    colnames(res) <- c('HR', 'hrneg', 'lower95', 'upper95', 'p')
    res<- data.frame(score=score,level=rownames(res), res, row.names = NULL)
    reslist[[score]]  <- res
  }
  resdf <- dplyr::bind_rows(reslist)

  resdf$score_rename <- fixlongnames(resdf$score)


  #significance...
  resdf$test <- ifelse(resdf$p < 0.05, 'significant', 'non-sig')
  resdf$test <- factor(resdf$test, levels = c('non-sig', 'significant'))

  # color...
  resdf$testcolor <- ifelse(resdf$test == 'significant', yes="#66C2A5", no = "#FC8D62")

  #factorize..
  # resdf$score_rename <- factor(resdf$score_rename, levels=resdf$score_rename)

  forestres <- ggplot(resdf, aes(y = level, x = HR, xmin=lower95, xmax=upper95, col=test))+
    geom_pointrange()+
    geom_vline(xintercept = 1, linetype = 'dashed')+
    facet_wrap(~score_rename, ncol=1, scales='free_y')+
    theme_linedraw()+
    theme(axis.text.y = element_text(size = 7),
          axis.title.y = element_blank() )+
    # scale_color_manual(values=resdf$testcolor)+
    labs(title = 'HR Forest plot, multivariable categorical')+
    xlab('Hazard Ratio\nLeft = better survival ; Right = worse survival')


  forestres

}

