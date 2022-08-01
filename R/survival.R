


# to do items
# consider saving things in the first survival function, just plotting in later funtions?
# finalize and check again categorical variable surv
# add in all later functions: summary plot, table, and forestplots

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
#' @param time - a string, the colname of the column in clinvardf with time to event information.
#' @param status - a string, the colname of the column in clinvardf with event information. The event column should be a character vector and should have the word "Censored" (capital C); ie alive/dead should be "Censored" and "Dead", event should be "Censored" and "Event".
#'
#' @return
#' @export
#'
#' @examples
survival <- function(testvardf,
                     clinvardf,
                     vartypes,
                     multivarnames,
                     time,
                     status){

  #if( missing(vartypes) ) {lfclimit <- 15}


  message("\nInitiating analysis")

  ### try to get continuous vs cat test vars ###
  # if vartpyes is provided (not missing), use that
  if( !missing(vartypes) ) {

    message("\nUsing vartypes to separate continuous and categorical variables...")

    if(length(vartypes) != ncol(testvardf)){error("vartypes is not the same length as ncol(testvardf). make sure each column is accounter for in vartypes as either 'continuous' or 'categorical', or leave vartypes blank.")}

    # using the vartypes, get the continuous / categorical columns
    contvars <- testvardf[,vartypes=='continuous']
    catvars <- testvardf[,vartypes=='categorical']

    message(' - ', ncol(contvars), ' continuous variables')
    message(' - ', ncol(catvars), ' categorical variables')

  }

  # if vartypes is missing, then try to guess based on continuous =  >12 unique values
  if( missing(vartypes) ) {

    message('\nGuessing vartypes based on continuous defined as > 12 unique values...')

    #get num unique values for each variable
    num_unique_vals <- sapply(md, FUN = function(x){length(unique(x))}, simplify = T)

    #categorical
    catnames <- names(num_unique_vals[num_unique_vals<=12])
    catvars <- testvardf[,catnames]

    #continuous
    contnames <- names(num_unique_vals[num_unique_vals>12])
    contvars <- testvardf[,contnames]

    message(' - ', ncol(contvars), ' continuous variables')
    message(' - ', ncol(catvars), ' categorical variables')

  }



  #prep output lists
  plotlist <- list()
  modellist <- list()
  datalist <- list()


  #### categorical vars ####

  if(ncol(catvars) > 0){

    message("\nBeginning analysis for categorical variables...")


    #get varnames
    scorenames <- names(catvars)

    # prep progress bar
    total = ncol(catvars)
    pb <- txtProgressBar(min = 0, max = total, style = 3)


    for(scoreidx in 1:length(scorenames)) {

      #get the score name
      scorename <- colnames(catvars)[scoreidx]

      #get the score valeues, numeric continuopus
      scorevec <- catvars[,scorename]

      # change name to have N
      scorevec <- plyr::mapvalues(scorevec, from = levels(scorevec),
                                  to = paste0(levels(scorevec), '\nN=', table(scorevec)[levels(scorevec)] ) )



      df <- data.frame(var = scorevec,
                       surv = md[,time],
                       status = md[,status])


      ### code status as numeric ###
      # EVENT = 1, CENSORED = 0
      df$status_code <- ifelse(df$status == 'Censored', yes=0, no = 1)



      fit <- survfit(Surv(surv, status_code) ~ var, data = df)

      gcat <- ggsurvplot(fit,
                         conf.int = F,
                         pval = TRUE,
                         pval.method = T,
                         linetype = "strata",
                         test.for.trend = F,
                         title=scorename,
                         ggtheme = theme_classic2())

      #test for trend, but only if >2 vars...
      if( length(levels(df$var)) > 2 )  {



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

      if( length(levels(df$var)) > 2 )  {
        plots <- list(gcat, gtft)

      } else{
        plots <- gcat
      }

      ### cox modelling ###
      #dich model, will match plot
      model <- coxph(Surv(surv, status_code) ~ var, data = df)


      #multivar cox test
      if( !missing(multivarnames) ){

        #prep for multivar
        df <- cbind(df, clinvardf)

        # run it
        multivar <- coxph(as.formula(
          paste0('Surv(surv, status_code) ~ ',
                 paste(multivarnames, collapse = ' + '))
        ),
        data=df)

      }



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
      datalist[[scorename]] <- df


      setTxtProgressBar(pb, scoreidx)

    } # end cont var loop.

  }





  #### continuous vars ####


  if(ncol(contvars) > 0){

    message("\nBeginning analysis for categorical variables...")

    #get varnames
    scorenames <- names(contvars)

    # prep progress bar
    total = ncol(contvars)
    pb <- txtProgressBar(min = 0, max = total, style = 3)


    for(scoreidx in 1:length(scorenames)) {

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



      ### set up survival df ###
      df <- data.frame(continuous = scorevec,
                       dichotomized = dichscore,
                       terciles = terciles,
                       surv = clinvardf[,time],
                       status = clinvardf[,status])



      # #check dist
      # #check dist over dead/alive
      # hist <- ggplot(df, aes(continuous, fill=status))+
      #   geom_histogram(col='black',position='identity', alpha=0.5) +
      #   facet_wrap(~status, nrow=2)
      #
      # #check dist cor with survival?
      #
      # corCensored <- cor.test(df[df$status=='Censored', "surv"], df[df$status=='Censored', "continuous"])
      # corNonCensored <- cor.test(df[df$status!='Censored', "surv"], df[df$status!='Censored', "continuous"])
      # corlab <- paste0('Non-Censored cor: ', round(corNonCensored$estimate, 3), ', Non-Censored P:', round(corNonCensored$p.value, 3),
      #                  '\nCensored cor: ', round(corCensored$estimate, 3), ', Censored P: ', round(corCensored$p.value, 3))
      # corp <- ggplot(df, aes(continuous, surv, col = status))+
      #   geom_point()+
      #   geom_smooth()+
      #   facet_wrap(~status, nrow=2)+
      #   labs(caption = corlab)
      #
      # distplots <- hist+corp

      ### code status as numeric ###
      # EVENT = 1, CENSORED = 0
      df$status_code <- ifelse(df$status == 'Censored', yes=0, no = 1)



      ### tertile analysis ###
      fit <- survfit(Surv(surv, status_code) ~ terciles, data = df)

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
      fit <- survfit(Surv(surv, status_code) ~ dichotomized, data = df)

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
      dich <- coxph(Surv(surv, status_code) ~ dichotomized, data = df)

      #terc
      terc <- coxph(Surv(surv, status_code) ~ terciles, data = df)

      #continuous univariate cox
      continuous <- coxph(Surv(surv, status_code) ~ continuous, data = df)



      #multivar cox test
      if( !missing(multivarnames) ){

        #prep for multivar
        df <- cbind(df, clinvardf)

        # run it
        multivar <- coxph(as.formula(
          paste0('Surv(surv, status_code) ~ ',
                 paste(multivarnames, collapse = ' + '))
        ),
        data=df)

      }



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
      datalist[[scorename]] <- df


      setTxtProgressBar(pb, scoreidx)

    } # end cont var loop.


  } # end cont var if statement.


}
