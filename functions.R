#########################################################
# EMA's 'Method A' (ANOVA)                              #
# fixed: sequence, subject(sequence), period, treatment #
#########################################################
method.A <- function(alpha = 0.05, path.in, path.out = tempdir(), file,
                     set = "", ext, na = ".", sep = ",", dec = ".",
                     logtrans = TRUE, regulator = "EMA", ola = FALSE,
                     print = TRUE, details = FALSE, adjust = FALSE,
                     verbose = FALSE, ask = FALSE, plot.bxp = FALSE,
                     fence = 2, data = NULL) {
  exec <- strftime(Sys.time(), usetz=TRUE)
  if (!missing(ext)) ext <- tolower(ext) # case-insensitive
  if (regulator == "HC")
    stop("For HC use method.B() option == 1 or option == 2 instead.")
  ret  <- CV.calc(alpha=alpha, path.in=path.in, path.out=path.out,
                  file=file, set=set, ext=ext, na=na, sep=sep,
                  dec=dec, logtrans=logtrans, regulator=regulator, ola=ola,
                  print=print, verbose=verbose, ask=ask,
                  plot.bxp=plot.bxp, fence=fence, data=data)
  results <- paste0(ret$res.file, "_MethodA.txt")
  # generate variables based on the attribute
  # 2nd condition: Otherwise, the header from a CSV file will be overwritten
  if (!is.null(data) & missing(ext)) {
    info  <- info.data(data)
    file  <- info$file
    set   <- info$set
    ref   <- info$ref
    # descr <- info$descr
    ext   <- ""
  }
  logtrans <- ret$logtrans
  os <- Sys.info()[[1]] # get OS for line-endings in output (Win: CRLF)
  ow <- options()       # save options
  options(digits=12)    # increase digits for anova()
  on.exit(ow)           # ensure that options are reset if an error occurs
  if (logtrans) {       # use the raw data and log-transform internally
    modA <- lm(log(PK) ~ sequence + subject%in%sequence + period + treatment,
               data = ret$data)
  } else {              # use the already log-transformed data
    modA <- lm(logPK ~ sequence + subject%in%sequence + period + treatment,
               data = ret$data)
  }

  PE  <- exp(coef(modA)[["treatmentT"]])
  CI  <- as.numeric(exp(confint(modA, "treatmentT", level=1-2*alpha)))
  DF  <- aov(modA)[[8]]
  res <- data.frame(ret$type, "A", ret$n, ret$nTT, ret$nRR,
                    paste0(ret$Sub.Seq, collapse="|"),
                    paste0(ret$Miss.seq, collapse="|"),
                    paste0(ret$Miss.per, collapse="|"), alpha,
                    DF, ret$CVwT, ret$CVwR, ret$swT, ret$swR, ret$sw.ratio,
                    ret$sw.ratio.upper, ret$BE1, ret$BE2, CI[1], CI[2],
                    PE, "fail", "fail", "fail", log(CI[2])-log(PE),
                    paste0(ret$ol, collapse="|"), ret$CVwR.rec,
                    ret$swR.rec, ret$sw.ratio.rec,
                    ret$sw.ratio.rec.upper, ret$BE.rec1,
                    ret$BE.rec2, "fail", "fail", "fail",
                    stringsAsFactors=FALSE)
  names(res)<- c("Design", "Method", "n", "nTT", "nRR", "Sub/seq",
                 "Miss/seq", "Miss/per", "alpha", "DF", "CVwT(%)",
                 "CVwR(%)", "swT", "swR", "sw.ratio", "sw.ratio.CL",
                 "L(%)", "U(%)", "CL.lo(%)", "CL.hi(%)", "PE(%)",
                 "CI", "GMR", "BE", "log.half-width", "outlier",
                 "CVwR.rec(%)", "swR.rec", "sw.ratio.rec",
                 "sw.ratio.rec.CL", "L.rec(%)", "U.rec(%)",
                 "CI.rec", "GMR.rec", "BE.rec")
  if (ret$BE2 == 1.25) { # change column names if necessary
    colnames(res)[which(names(res) == "L(%)")] <- "BE.lo(%)"
    colnames(res)[which(names(res) == "U(%)")] <- "BE.hi(%)"
  }
  if (!is.na(ret$BE.rec2) & ret$BE.rec2 == 1.25) { # change column names if necessary
    colnames(res)[which(names(res) == "L.rec(%)")] <- "BE.rec.lo(%)"
    colnames(res)[which(names(res) == "U.rec(%)")] <- "BE.rec.hi(%)"
  }
  # Convert CVs, limits, PE, and CI (till here as fractions) to percent
  res$"CVwT(%)" <- 100*res$"CVwT(%)"
  res$"CVwR(%)" <- 100*res$"CVwR(%)"
  if ("BE.lo(%)" %in% names(res)) { # conventional limits
    res$"BE.lo(%)" <- 100*res$"BE.lo(%)"
    res$"BE.hi(%)" <- 100*res$"BE.hi(%)"
  } else {                          # expanded limits
    res$"L(%)" <- 100*res$"L(%)"
    res$"U(%)" <- 100*res$"U(%)"
  }
  if (!is.na(res$"CVwR.rec(%)")) {
    res$"CVwR.rec(%)"    <- 100*res$"CVwR.rec(%)"
    if ("BE.rec.lo(%)" %in% names(res)) { # conventional limits
      res$"BE.rec.lo(%)" <- 100*res$"BE.rec.lo(%)"
      res$"BE.rec.hi(%)" <- 100*res$"BE.rec.hi(%)"
    } else {                              # expanded limits
      res$"L.rec(%)" <- 100*res$"L.rec(%)"
      res$"U.rec(%)" <- 100*res$"U.rec(%)"
    }
  }
  res$"PE(%)"    <- 100*res$"PE(%)"
  res$"CL.lo(%)" <- 100*res$"CL.lo(%)"
  res$"CL.hi(%)" <- 100*res$"CL.hi(%)"
  if (res$"CL.lo(%)" >= 100*ret$BE1 &
      res$"CL.hi(%)" <= 100*ret$BE2)
    res$CI <- "pass"  # CI within acceptance range
  if (res$"PE(%)" >= 80 & res[["PE(%)"]] <= 125)
    res$GMR <- "pass" # PE within 80.00-125.00%
  if (res$CI == "pass" & res$GMR == "pass")
    res$BE <- "pass"  # if passing both, conclude BE
  if (!is.na(res$"CVwR.rec(%)")) {
    if (res$"CL.lo(%)" >= 100*ret$BE.rec1 &
        res$"CL.hi(%)" <= 100*ret$BE.rec2)
      res$CI.rec <- "pass"  # CI within acceptance range
    res$GMR.rec <- res$GMR
    if (res$CI.rec == "pass" & res$GMR.rec == "pass")
      res$BE.rec <- "pass"  # if passing both, conclude BE
  }
  if (details) { # results in full precision
    ret <- res
    if (as.character(res$outlier) == "NA") {
      # remove superfluous columns if ola=FALSE or ola=TRUE
      # and no outlier(s) detected
      ret <- ret[, !names(ret) %in% c("outlier", "CVwR.rec(%)",
                                      "swR.rec", "sw.ratio.rec",
                                      "L.rec(%)", "U.rec(%)",
                                      "CI.rec", "GMR.rec", "BE.rec")]
    }
    #class(ret) <- "repBE"
    return(ret)
  }
  # Round percents to two decimals according to the GL
  res$"CVwT(%)" <- res$"CVwT(%)"
  res$"CVwR(%)" <- res$"CVwR(%)"
  if ("BE.lo(%)" %in% names(res)) { # conventional limits
    res$"BE.lo(%)" <- res$"BE.lo(%)"
    res$"BE.hi(%)" <- res$"BE.hi(%)"
  } else {                          # expanded limits
    res$"L(%)" <- res$"L(%)"
    res$"U(%)" <- res$"U(%)"
  }
  res$"PE(%)"    <- res$"PE(%)"
  res$"CL.lo(%)" <- res$"CL.lo(%)"
  res$"CL.hi(%)" <- res$"CL.hi(%)"
  if (!is.na(res$"CVwR.rec(%)")) {
    res$"CVwR.rec(%)" <- res$"CVwR.rec(%)"
    if ("BE.rec.lo(%)" %in% names(res)) { # conventional limits
      res$"BE.rec.lo(%)" <- res$"BE.rec.lo(%)"
      res$"BE.rec.hi(%)" <- res$"BE.rec.hi(%)"
    } else {                          # expanded limits
      res$"L.rec(%)" <- res$"L.rec(%)"
      res$"U.rec(%)" <- res$"U.rec(%)"
    }
  }
  overwrite <- TRUE # default
  if (print) { # to file in UTF-8
    if (ask & file.exists(results)) {
      answer <- tolower(readline("Results already exists. Overwrite the file [y|n]? "))
      if(answer != "y") overwrite <- FALSE
    }
    if (overwrite) { # either the file does not exist or should be overwritten
      # only binary mode supports UTF-8 and different line endings
      res.file <- file(description=results, open="wb")
      res.str  <- info.env(fun="method.A", option=NA, path.in, path.out,
                           file, set, ext, exec, data)
      if (os == "Windows") res.str <- gsub("\n", "\r\n", res.str) # CRLF
      if (os == "Darwin")  res.str <- gsub("\n", "\r", res.str)   # CR
      writeBin(charToRaw(res.str), res.file)
      close(res.file)
    }
  }
  # insert DF and alpha into the text generated by CV.calc()
  cut.pos    <- unlist(gregexpr(pattern="Regulator", ret$txt))
  left.str   <- substr(ret$txt, 1, cut.pos-1)
  right.str  <- substr(ret$txt, cut.pos, nchar(ret$txt))
  insert.str <- paste0("Degrees of freedom : ", sprintf("%3i", DF),
                       "\nalpha              :   ", alpha,
                       " (", 100*(1-2*alpha), "% CI)\n")
  txt <- paste0(left.str, insert.str, right.str)
  if (!is.na(res$"CVwR.rec(%)")) {
    txt1 <- paste0("\n\nAssessment based on original CVwR",
                   sprintf(" %.2f%%", res$"CVwR(%)"))
    txt <- paste0(txt, txt1, "\n", paste0(rep("\u2500", nchar(txt1)-2), collapse=""))
  }
  txt <- paste0(txt,
                "\nConfidence interval: ", sprintf("%6.2f%% ... %6.2f%%",
                                                   res$"CL.lo(%)", res$"CL.hi(%)"), "  ", res$CI,
                "\nPoint estimate     : ", sprintf("%6.2f%%", res$"PE(%)"),
                "              ", res$GMR,
                "\nMixed (CI & PE)    :                      ", res$BE, "\n")
  txt <- paste0(txt,
                repBE.draw.line(called.from="ABEL", L=ret$BE1, U=ret$BE2,
                                lo=CI[1], hi=CI[2], PE=PE), "\n")
  if (!is.na(res$"CVwR.rec(%)")) {
    txt1 <- paste0("\nAssessment based on recalculated CVwR",
                   sprintf(" %.2f%%", res$"CVwR.rec(%)"))
    txt <- paste0(txt, txt1, "\n", paste0(rep("\u2500", nchar(txt1)-1), collapse=""))
    txt <- paste0(txt, "\nConfidence interval: ", res$CI.rec,
                  "\nPoint estimate     : ", res$GMR.rec,
                  "\nMixed (CI & PE)    : ", res$BE.rec, "\n")
    txt <- paste0(txt,
                  repBE.draw.line(called.from="ABEL", L=ret$BE.rec1, U=ret$BE.rec2,
                                  lo=CI[1], hi=CI[2], PE=PE), "\n")
  }
  if (res$Design == "TRR|RTR")
    txt <- paste0(txt, "Note: The extra-reference design assumes lacking period effects. ",
                  "The treatment\ncomparison will be biased in the presence of a ",
                  "true period effect.\n")
  if (print & overwrite) {
    res.file <- file(description=results, open="ab")
    res.str  <- txt                                             # UNIXes LF
    if (os == "Windows") res.str <- gsub("\n", "\r\n", res.str) # CRLF
    if (os == "Darwin")  res.str <- gsub("\n", "\r", res.str)   # CR
    writeBin(charToRaw(res.str), res.file)
    close(res.file)
  }
  # Optional: Assess whether there is an inflation of the Type I Error
  # with the specified alpha. If yes, iteratively adjust alpha to
  # preserve the consumer risk.
  if (adjust) {
    if (!ret$type %in% c("TRTR|RTRT", "TRT|RTR", "TRR|RTR|RRT")) {
      message("Assessment of the Type I Error for this design is not implemented.")
    } else {
      if (PE > 0.80 & PE < 1.25) { # PE must be within limits!
        if (ret$type == "TRTR|RTRT")   des <- "2x2x4"
        if (ret$type == "TRT|RTR")     des <- "2x2x3"
        if (ret$type == "TRR|RTR|RRT") des <- "2x3x3"
        if (print) {
          txt <- paste0("\n", paste0(rep("\u2500", 74), collapse=""))
          adj <- scABEL.ad(theta0=PE, CV=ret$CVwR, design=des,
                           n=ret$Sub.Seq, alpha.pre=alpha, print=FALSE)
          if (!is.na(ret$CVwR.rec)) { # 2nd run for recalculated CVwR
            adj1 <- scABEL.ad(theta0=PE, CV=ret$CVwR.rec, design=des,
                              n=ret$Sub.Seq, alpha.pre=alpha, print=FALSE)
          }
          no.infl <- "  TIE not > nominal 0.05; consumer risk is controlled.\n"
          ifelse (is.na(ret$CVwR.rec),
                  txt <- paste0(txt, "\nAssessment of the empiric Type I Error (TIE); "),
                  txt <- paste0(txt, "\nAssessment of the empiric Type I Error (TIE) based on original CVwR;\n"))
          iter <- (adj$sims - adj$sims %% 1e6) / 1e6
          ifelse (iter == 1,
                  txt <- paste0(txt, "1,000,000 studies simulated.\n"),
                  txt <- paste0(txt, "1,000,000 studies in each of the ", iter,
                                " iterations simulated.\n"))
          if (is.na(adj$alpha.adj)) {
            txt <- paste0(txt, no.infl)
          } else {
            txt <- paste0(txt, "  TIE for alpha",
                          sprintf(" %1.6f         : %1.5f",
                                  alpha, adj$TIE.unadj), "\n")
            txt <- paste0(txt, "  TIE for adjusted alpha",
                          sprintf(" %1.6f: %1.5f",
                                  adj$alpha.adj, adj$TIE.adj), "\n")
          } # EO orginal CVwR
          if (!is.na(ret$CVwR.rec)) {
            txt <- paste0(txt, "Assessment of the empiric Type I Error (TIE) based on recalculated CVwR;")
            iter <- (adj1$sims - adj1$sims %% 1e6) / 1e6
            ifelse (iter == 1,
                    txt <- paste0(txt, "\n1,000,000 studies simulated.\n"),
                    txt <- paste0(txt, "\n1,000,000 studies in each of the ", iter,
                                  " iterations simulated.\n"))
            if (is.na(adj1$alpha.adj)) {
              txt <- paste0(txt, no.infl)
            } else {
              if (alpha)
                txt <- paste0(txt, "  TIE for alpha",
                              sprintf(" %1.6f         : %1.5f",
                                      alpha, adj1$TIE.unadj), "\n")
              if (round(adj1$TIE.adj, 5) != round(adj1$TIE.unadj, 5)) { # this should not happen...
                txt <- paste0(txt, "  TIE for adjusted alpha",
                              sprintf(" %1.6f: %1.5f", adj1$alpha.adj, adj1$TIE.adj), "\n")
              } else {
                txt <- paste0(txt, no.infl)
              }
            }
          } # EO recalculated CVwR
          if (overwrite) {
            res.file <- file(results, open="ab")                        # line endings
            res.str  <- txt                                             # LF (UNIXes, Solaris)
            if (os == "Windows") res.str <- gsub("\n", "\r\n", res.str) # CRLF (Windows)
            if (os == "Darwin")  res.str <- gsub("\n", "\r", res.str)   # CR (OSX)
            writeBin(charToRaw(res.str), res.file)
            close(res.file)
          }
        } else { # to console
          scABEL.ad(theta0=PE, CV=ret$CVwR, design=des, n=ret$Sub.Seq,
                    alpha.pre=alpha, details=TRUE)
          if (!is.na(ret$CVwR.rec)) { # 2nd run for recalculated CVwR
            scABEL.ad(theta0=PE, CV=ret$CVwR.rec, design=des, n=ret$Sub.Seq,
                      alpha.pre=alpha, details=TRUE)
          }
        }
      } else {
        message("PE must be within 80.00\u2013125.00% for assessment of the Type I Error.")
      }
    }
  } # end of iteratively adjusting alpha
  
  res_and_anova <- list()
  res_and_anova$res <- res
  if (verbose) {
    name <-  paste0(file, set)
    len  <- max(27+nchar(name), 35)
    # cat(paste0("\nData set ", name, ": Method A by lm()"),
        # paste0("\n", paste0(rep("\u2500", len), collapse = "")), "\n")
    # change from type I (default as in versions up to 1.0.17)
    # to type III to get the correct carryover test
    typeIII <- stats::anova(modA) # otherwise summary of lmerTest is used
    attr(typeIII, "heading")[1] <- "Type III Analysis of Variance Table\n"
    MSdenom <- typeIII["sequence:subject", "Mean Sq"]
    df2     <- typeIII["sequence:subject", "Df"]
    fvalue  <- typeIII["sequence", "Mean Sq"] / MSdenom
    df1     <- typeIII["sequence", "Df"]
    typeIII["sequence", 4] <- fvalue
    typeIII["sequence", 5] <- pf(fvalue, df1, df2, lower.tail = FALSE)
    res_and_anova$anova <- typeIII
  }
  return(res_and_anova)
} # end of function method.A()

####################################################
# Calculate CV according to the EMA's Q&A-document #
####################################################
CV.calc <- function(alpha = 0.05, path.in, path.out, file, set = "",
                    ext, na, sep = ",", dec = ".", logtrans = TRUE,
                    regulator, ola = FALSE, details = FALSE,
                    adjust = FALSE, print, verbose = FALSE, ask = FALSE,
                    theta1 = theta1, theta2 = theta2, plot.bxp = FALSE,
                    fence = 2, data) {
  if (missing(path.in)) path.in <- NULL
  if (missing(data)) data <- NULL
  called.from <- as.character(sys.call(-1))[1]
  ret   <- get.data(path.in = path.in, path.out = path.out, file = file,
                    set = set, ext = ext, na = na, sep = sep, dec = dec,
                    logtrans = logtrans, print = print,
                    plot.bxp = plot.bxp, data = data)
  logtrans <- ret$logtrans
  ow <- options("digits")
  options(digits=12) # more digits for anova
  on.exit(ow)        # ensure that options are reset if an error occurs
  if (logtrans) {    # use raw data and log-transform internally
    modCVR <- lm(log(PK) ~ sequence + subject%in%sequence + period,
                 data = ret$ref)
  } else {           # use already log-transformed data
    modCVR <- lm(logPK ~ sequence + subject%in%sequence + period,
                 data = ret$ref)
  }
  aovCVR <- anova(modCVR)
  msewR  <- aovCVR["Residuals", "Mean Sq"]
  DfCVR  <- aovCVR["Residuals", "Df"]
  CVwR   <- mse2CV(msewR)
  if (ret$design == "full") { # not for the EMA but the WHO
    if (logtrans) {
      modCVT <- lm(log(PK) ~ sequence + subject%in%sequence + period,
                   data = ret$test)
    } else {
      modCVT <- lm(logPK ~ sequence + subject%in%sequence + period,
                   data = ret$test)
    }
    aovCVT <- anova(modCVT)
    msewT  <- aovCVT["Residuals", "Mean Sq"]
    DfCVT  <- aovCVT["Residuals", "Df"]
    CVwT  <- mse2CV(msewT)
    sw.ratio <- sqrt(msewT)/sqrt(msewR)
    sw.ratio.CI <- c(sw.ratio/sqrt(qf(0.1/2, df1 = DfCVT, df2 = DfCVR,
                                      lower.tail = FALSE)),
                     sw.ratio/sqrt(qf(1-0.1/2, df1 = DfCVT, df2 = DfCVR,
                                      lower.tail = FALSE)))
    names(sw.ratio.CI) <- c("lower", "upper")
  }
  outlier <- FALSE
  BE.rec  <- rep(NA, 2)
  if (ola & !called.from == "ABE") { # check for outliers
    stud.res  <- rstudent(modCVR)  # studentized (SAS)
    stud.res  <- stud.res[!is.na(stud.res)] # get rid of NAs and zeros
    stud.res  <- stud.res[which(stud.res != 0)]
    stud.res  <- stud.res[c(TRUE, FALSE)] # need only the 1st occasions
    bp1       <- boxplot(stud.res, range=fence, plot=FALSE)
    names.ol1 <- names(bp1$out)
    stand.res <- rstandard(modCVR) # standardized (SAS, PHX/WNL)
    stand.res <- stand.res[!is.na(stand.res)]
    stand.res <- stand.res[which(stand.res != 0)]
    stand.res <- stand.res[c(TRUE, FALSE)] # need only the 1st occasions
    bp2       <- boxplot(stand.res, range=fence, plot=FALSE)
    names.ol2 <- names(bp2$out)
    if (length(names.ol1) > 0) { # >=1 detected
      outlier   <- TRUE
      ol.value1 <- as.numeric(bp1$out)
      ol.seq1   <- as.character(ret$ref[names(bp1$out), "sequence"])
      ol.subj1  <- as.character(ret$ref[names(bp1$out), "subject"])
      ol.value2 <- as.numeric(bp2$out)
      ol.seq2   <- as.character(ret$ref[names(bp2$out), "sequence"])
      ol.subj2  <- as.character(ret$ref[names(bp2$out), "subject"])
      pars <- list(boxwex = 0.5, boxfill = "lightblue", medcol = "blue",
                   outpch = 21, outcex = 1.35, outcol = "red", outbg = "#FFCCCC")
      overwrite <- TRUE
      if (as.logical(capabilities("png"))) {
        if (plot.bxp) {   # save in PNG format to path.out
          if (ask & file.exists(ret$png.path)) {
            answer <- tolower(readline("Boxplot already exists. Overwrite the PNG  [y|n]? "))
            if(answer != "y") overwrite <- FALSE
          }
          if (overwrite) { # either the file does not exist or should be overwritten
            png(ret$png.path, width = 720, height = 720, pointsize = 18)
          }
        }
      } else {
        message("png-device is not available; changed to plot.bxp = FALSE")
      }
      bxp(bp1, xlim = c(0, 3), ylim = c(-1, 1)*max(abs(c(ol.value1, ol.value2))),
          las = 1, ylab = "residual", pars = pars, main = "")
      title(main = expression(paste("EMA\u2019s model for ", CV[wR], ":")),
            line = 3, cex.main = 1.1)
      title(main = expression(paste("log(response) ~ sequence + subject(sequence) + period; data = R")),
            line = 2, cex.main = 1.05)
      title(main = bquote(paste("Outlier fence ", .(fence), "\u00D7IQR")),
            line = 1, cex.main = 1.05)
      if (length(names.ol1) == 0) {
        lab.txt <- "no outlier"
      } else {
        ifelse (length(ol.value1) == 1,
                lab.txt <- "1 outlier",
                lab.txt <- paste(length(ol.value1), "outliers"))
      } # Note: not /exactly/ equal. SAS uses 'type=2'
      mtext(paste0("studentized\n(R, ~SAS)\n\n", lab.txt), 1, line = 4, at = 1)
      text(rep(1.1, 2), bp1$stats[c(1, 5)], adj = c(0, 0.25), cex = 0.8,
           sprintf("%+.3f", bp1$stats[c(1, 5)]))
      if (!identical(ol.value1, numeric(0))) { # only if stud. outlier
        text(rep(0.9, length(ol.value1)), ol.value1, adj = c(1, 0.25), cex = 0.8,
             paste0("# ", ol.subj1, " (", ol.seq1, ")"))
        text(rep(1.1, length(ol.value1)), ol.value1, adj = c(0, 0.25), cex  =0.8,
             sprintf("%+.3f", ol.value1))
      }
      bxp(bp2, axes = FALSE, at = 2, add = TRUE, pars = pars)
      if (length(names.ol2) == 0) {
        lab.txt <- "no outlier"
      } else {
        ifelse (length(ol.value2) == 1,
                lab.txt <- "1 outlier",
                lab.txt <- paste(length(ol.value2), "outliers"))
      } # Note: not /exactly/ equal. SAS uses 'type=2' and PHX/WNL 'type=6'
      mtext(paste0("standardized\n(R, ~SAS,\n~Phoenix WinNonlin)\n", lab.txt), 1,
            line = 4, at = 2)
      text(rep(2.1, 2), bp2$stats[c(1, 5)], adj = c(0, 0.25), cex = 0.8,
           sprintf("%+.3f", bp2$stats[c(1, 5)]))
      if (!identical(ol.value2, numeric(0))) { # only if stand. outlier
        text(rep(1.9, length(ol.value2)), ol.value2, adj = c(1, 0.25), cex = 0.8,
             paste0("# ", ol.subj2, " (", ol.seq2, ")"))
        text(rep(2.1, length(ol.value2)), ol.value2, adj = c(0, 0.25), cex = 0.8,
             sprintf("%+.3f", ol.value2))
      }
      abline(h = 0, lty = "dotted")
      if (plot.bxp && file.exists(ret$png.path)) {
        if (!is.null(dev.list()["png"])) {
          invisible(dev.off(dev.list()["png"]))
        }
      }
      # recalculate CVwR for data without outlier(s)
      ol   <- ret$ref[names(bp1$out), "subject"]
      excl <- ret$ref[!ret$ref$subject %in% ol, ]
      if (logtrans) { # use the raw data and log-transform internally
        modCVR.rec <- lm(log(PK) ~ sequence + subject%in%sequence + period,
                         data = excl)
      } else {
        modCVR.rec <- lm(logPK ~ sequence + subject%in%sequence + period,
                         data = excl)
      }
      aovCVR.rec <- anova(modCVR.rec)
      msewR.rec  <- aovCVR.rec["Residuals", "Mean Sq"]
      DfCVR.rec  <- aovCVR.rec["Residuals", "Df"]
      CVwR.rec   <- mse2CV(msewR.rec)
      if (ret$design == "full") {
        sw.ratio.rec <- sqrt(msewT)/sqrt(msewR.rec)
        sw.ratio.rec.CI <- c(sw.ratio.rec/sqrt(qf(0.1/2, df1 = DfCVT,
                                                  df2 = DfCVR.rec,
                                                  lower.tail = FALSE)),
                             sw.ratio.rec/sqrt(qf(1-0.1/2, df1 = DfCVT,
                                                  df2 = DfCVR.rec,
                                                  lower.tail = FALSE)))
        names(sw.ratio.rec.CI) <- c("lower", "upper")
      }
      if (verbose) {
        stud.res.whiskers <- signif(range(bp1$stats[, 1]), 7)
        stud.res.outliers <- data.frame(ol.subj1, ol.seq1, signif(ol.value1, 7))
        names(stud.res.outliers) <- c("subject", "sequence", "stud.res")
        stand.res.whiskers <- signif(range(bp2$stats[, 1]), 7)
        stand.res.outliers <- data.frame(ol.subj2, ol.seq2, signif(ol.value2, 7))
        names(stand.res.outliers) <- c("subject", "sequence", "stand.res")
        cat(paste0("\nOutlier analysis\n (externally) studentized residuals",
                   "\n Limits (", fence, "\u00D7IQR whiskers): ",
                   stud.res.whiskers[1], ", ", stud.res.whiskers[2],
                   "\n Outliers:\n")); print(stud.res.outliers, row.names = FALSE)
        cat(paste0("\n standarized (internally studentized) residuals\n Limits (",
                   fence, "\u00D7IQR whiskers): ", stand.res.whiskers[1], ", ",
                   stand.res.whiskers[2], "\n Outliers:\n"))
        # since standardized residuals are more liberal,
        # we have to deal with such a special case
        if (nrow(stand.res.outliers) == 0) {
          cat(" none detected\n")
        } else {
          print(stand.res.outliers, row.names = FALSE)
        }
      } # EO verbose
    } # EO >= 1 outlier
  } # EO outlier analysis (only if called from method.A()/method.B() & ola=TRUE)
  if (!called.from == "ABE") {
    if (regulator == "EMA") {
      reg_set <- reg_const("EMA")
      BE <- as.numeric(scABEL(CV = CVwR, regulator = "EMA"))
    }
    if (regulator == "GCC") {
      reg_set <- reg_const("GCC")
      BE <- as.numeric(scABEL(CV = CVwR, regulator = "GCC"))
    }
    if (regulator == "HC") {
      reg_set <- reg_const("HC")
      if (!alpha == 0.5) {
        BE <- as.numeric(scABEL(CV = CVwR, regulator = "HC"))
      } else { # Cmax
        BE <- c(0.80, 1.25)
      }
    }
  } else {
    BE <- c(theta1, theta2)
  }
  txt <- ret$txt
  if (called.from != "ABE") { # only for scaling
    if (regulator == "EMA") {
      txt <- paste0(ret$txt,
                    "\nRegulator          : EMA",
                    "\nSwitching CV       : ",
                    sprintf("%6.2f%%", 100*reg_set$CVswitch),
                    "\nScaling cap        : ",
                    sprintf("%6.2f%%", 100*reg_set$CVcap),
                    "\nRegul. constant (k):  ",
                    sprintf("%.3f", reg_set$r_const))
    }
    if (regulator == "GCC") {
      txt <- paste0(ret$txt,
                    "\nRegulator          : GCC",
                    "\nSwitching CV       : ",
                    sprintf("%6.2f%%", 100*reg_set$CVswitch))
    }
    if (regulator == "HC") {
      if (!alpha == 0.5) {
        txt <- paste0(ret$txt,
                      "\nRegulator          : HC",
                      "\nSwitching CV       : ",
                      sprintf("%6.2f%%", 100*reg_set$CVswitch),
                      "\nScaling cap        : ",
                      sprintf("%6.2f%%", 100*reg_set$CVcap),
                      "\nRegul. constant (k):  ",
                      sprintf("%.3f", reg_set$r_const))
      } else { # Cmax
        txt <- paste0(ret$txt,
                      "\nRegulator          : HC")
        
      }
    }
    if (CVwR > 0.3 & !regulator == "HC" & !alpha == 0.5) {
      txt <- paste0(txt, "\nGMR restriction    :  80.00% ... 125.00%")
    }
  }
  if (ret$design == "full") {
    # sw.ratio <- CV2se(CVwT)/CV2se(CVwR) # we should already have it, right?
    txt <- paste0(txt, "\nCVwT               : ",
                  sprintf("%6.2f%%", 100*CVwT))
    if (called.from != "ABE") { # not needed for ABE
      txt <- paste0(txt, "\nswT                :   ",
                    sprintf("%.5f", CV2se(CVwT)))
    }
  }
  txt <- paste0(txt, "\nCVwR               : ", sprintf("%6.2f%%", 100*CVwR))
  if (called.from != "ABE") { # only for scaling
    txt <- paste0(txt, " (reference-scaling ")
    if (CVwR <= 0.3) txt <- paste0(txt, "not ")
    txt <- paste0(txt, "applicable)")
    if (CVwR <= 0.3) {
      txt <- paste0(txt, "\nBE-limits          :  80.00% ... 125.00%")
    } else {
      txt <- paste0(txt, "\nswR                :   ",
                    sprintf("%.5f", CV2se(CVwR)))
      if (regulator == "EMA") {
        txt <- paste0(txt, "\nExpanded limits    : ",
                      sprintf("%6.2f%% ... %.2f%%",
                              100*BE[1], 100*BE[2]), " [100exp(\u00B1",
                      sprintf("%.3f", reg_set$r_const), "\u00B7swR)]")
      }
      if (regulator == "GCC") {
        txt <- paste0(txt, "\nWidened limits      :",
                      sprintf("%6.2f%% ... %.2f%%",
                              100*BE[1], 100*BE[2]))
      }
    }
    if (ret$design == "full") {
      txt <- paste0(txt, "\nswT / swR          :   ",
                    sprintf("%.4f", sw.ratio))
      if (sw.ratio >= 2/3 & sw.ratio <= 3/2) { # like in PBE/IBE
        txt <- paste0(txt, " (similar variabilities of T and R)")
      } else {
        ifelse (sw.ratio < 2/3,
                txt <- paste0(txt, " (T lower variability than R)"),
                txt <- paste0(txt, " (T higher variability than R)"))
      }
      txt <- paste0(txt, "\nsw-ratio (upper CL):   ",
                    sprintf("%.4f", sw.ratio.CI[["upper"]]))
      if (sw.ratio.CI[["upper"]] <= 2.5) { # like in the FDA's warfarin guidance
        txt <- paste0(txt, " (comparable variabilities of T and R)")
      } else {
        txt <- paste0(txt, " (T higher variability than R)")
      }
    }
    if (ola) {
      if (outlier) {
        if (ret$design == "full") sw.ratio.rec <- sqrt(msewT)/sqrt(msewR.rec)
        BE.rec <- as.numeric(scABEL(CV=CVwR.rec, regulator="EMA"))
        txt <- paste0(txt, "\n\nOutlier fence      :  ", fence,
                      "\u00D7IQR of studentized residuals.")
        txt1 <- paste0("\nRecalculation due to presence of ",
                       length(ol))
        ifelse (length(ol) == 1,
                txt1 <- paste0(txt1, " outlier (subj. "),
                txt1 <- paste0(txt1, " outliers (subj. "))
        txt1 <- paste0(txt1, paste0(ol, collapse="|"), ")")
        txt <- paste0(txt, txt1, "\n",
                      paste0(rep("\u2500", nchar(txt1)-1), collapse=""))
        txt <- paste0(txt,
                      "\nCVwR (outl. excl.) : ", sprintf("%6.2f%%", 100*CVwR.rec),
                      " (reference-scaling ")
        if (CVwR.rec > 0.3) {
          txt <- paste0(txt, "applicable)")
          txt <- paste0(txt, "\nswR (recalculated) :   ",
                        sprintf("%.5f", CV2se(CVwR.rec)))
          txt <- paste0(txt, "\nExpanded limits    : ",
                        sprintf("%6.2f%% ... %.2f%%",
                                100*BE.rec[1], 100*BE.rec[2]), " [100exp(\u00B1",
                        sprintf("%.3f", reg_set$r_const), "\u00B7swR)]")
        } else {
          txt <- paste0(txt, "not applicable)")
          txt <- paste0(txt, "\nUnscaled BE-limits :  80.00% ... 125.00%")
        }
        if (ret$design == "full") {
          txt <- paste0(txt, "\nswT / swR (recalc.):   ",
                        sprintf("%.4f", sw.ratio.rec))
          if (sw.ratio.rec >= 2/3 & sw.ratio.rec <= 3/2) { # like in PBE/IBE
            txt <- paste0(txt, " (similar variabilities of T and R)")
          } else {
            ifelse (sw.ratio.rec < 2/3,
                    txt <- paste0(txt, " (T lower variability than R)"),
                    txt <- paste0(txt, " (T higher variability than R)"))
          }
          txt <- paste0(txt, "\nsw-ratio (upper CL):   ",
                        sprintf("%.4f", sw.ratio.rec.CI[["upper"]]))
          if (sw.ratio.rec.CI[["upper"]] <= 2.5) { # like in the FDA's warfarin guidance
            txt <- paste0(txt, " (comparable variabilities of T and R)")
          } else {
            txt <- paste0(txt, " (T higher variability than R)")
          }
        }
      } else {
        txt <- paste0(txt, "\n\nOutlier fence      :  ", fence,
                      "\u00D7IQR of studentized residuals.")
        txt <- paste0(txt, "\nNo outlier detected.",
                      "\n", paste0(rep("\u2500", 49), collapse=""))
      }
    } # EO ola
  } else { # called from ABE()
    txt <- paste0(txt, "\nBE-limits          : ",
                  sprintf("%6.2f%% ... %.2f%%", 100*BE[1], 100*BE[2]))
  }
  ret$txt <- txt
  ret <- c(ret, CVswitch=ifelse(called.from != "ABE", reg_set$CVswitch, NA),
           CVcap=ifelse(called.from != "ABE", reg_set$CVcap, NA),
           r_const=ifelse(called.from != "ABE", reg_set$r_const, NA), BE=BE,
           CVwT=ifelse(ret$design == "full", CVwT, NA), CVwR=CVwR,
           swT=ifelse(ret$design == "full", sqrt(msewT), NA), swR=sqrt(msewR),
           sw.ratio=ifelse(called.from != "ABE" & ret$design == "full",
                           sw.ratio, NA),
           sw.ratio.upper=ifelse(called.from != "ABE" & ret$design == "full",
                                 sw.ratio.CI[["upper"]], NA),
           ol=ifelse(called.from != "ABE" & outlier, list(ol.subj1), NA),
           CVwR.rec=ifelse(called.from != "ABE" & outlier, CVwR.rec, NA),
           swR.rec=ifelse(called.from != "ABE" & outlier, sqrt(msewR.rec), NA),
           sw.ratio.rec=ifelse(called.from != "ABE" & outlier &
                                 ret$design == "full", sw.ratio.rec, NA),
           sw.ratio.rec.upper=ifelse(called.from != "ABE" & outlier &
                                       ret$design == "full",
                                     sw.ratio.rec.CI[["upper"]], NA),
           BE.rec=BE.rec)
  return(ret)
} # end of function CV.calc()

################################################
# Information about the computing environment, #
# packages, data set, and method.              #
################################################
info.env <- function(fun, option=NA, path.in, path.out,
                     file, set, ext, exec, data) {
  if (!is.null(data) & missing(ext)) {
    info  <- info.data(data)
    file  <- info$file
    set   <- info$set
    ref   <- info$ref
    # descr <- info$descr
    ext   <- ""
  }
  ext.xls <- c("XLS", "xls", "XLSX", "xlsx")
  system <- Sys.info()
  node   <- system[["nodename"]]
  user   <- system[["user"]]
  OS     <- system[["sysname"]]
  OSrel  <- system[["release"]]
  OSver  <- system[["version"]]
  rver   <- sessionInfo()$R.version$version.string
  rver   <- substr(rver, which(strsplit(rver, "")[[1]]=="n")+2, nchar(rver))
  ryear  <- substr(rver, which(strsplit(rver, "")[[1]]=="(")-1, nchar(rver))
  rver   <- substr(rver, 1, which(strsplit(rver, "")[[1]]=="(")-2)
  cit    <- citation("replicateBE")
  year   <- paste0(" (", substr(cit, regexpr("year", cit)+8,
                                regexpr("year", cit)+11), ")")
  cit    <- citation("readxl")
  year1  <- paste0(" (", substr(cit, regexpr("year", cit)+8,
                                regexpr("year", cit)+11), ")")
  cit    <- citation("PowerTOST")
  year2  <- paste0(" (", substr(cit, regexpr("year", cit)+8,
                                regexpr("year", cit)+11), ")")
  cit    <- citation("nlme")
  year3  <- paste0(" (", substr(cit, regexpr("year", cit)+8,
                                regexpr("year", cit)+11), ")")
  cit    <- citation("lmerTest")
  year4  <- paste0(" (", substr(cit, regexpr("year", cit)+8,
                                regexpr("year", cit)+11), ")")
  lic    <- paste0("This code is copyright \u00A9 by Helmut Sch\u00FCtz, Michael Tomashevskiy, Detlew Labes.\n",
                   "This code is open source; you can redistribute it and/or modify it under the\n",
                   "terms of the GNU General Public License as published by the Free Software Foun-\n",
                   "dation; either version 3, or (at your option) any later version. See the GNU\n",
                   "GPL for more details. Copies of the GPL-3 versions are available at:\n",
                   "https://www.gnu.org/licenses/gpl-3.0.html")
  discl  <- paste0("\n\u2554", paste0(rep("\u2550", 76), collapse=""), "\u2557\n",
                   "\u2551 Program offered for Use without any Guarantees and Absolutely No Warranty. \u2551\n",
                   "\u2551 No Liability is accepted for any Loss and Risk to Public Health Resulting  \u2551\n",
                   "\u2551 from Use of this R-Code.                                                   \u2551\n",
                   "\u255A", paste0(rep("\u2550", 76), collapse=""), "\u255D")
  if (fun == "model.B") {
    ifelse (option == 1, hr.len <- 62+nchar(exec), hr.len <- 57+nchar(exec))
  } else {
    hr.len <- 79
  }
  hr     <- paste0(rep("\u2500", hr.len), collapse="")
  if (!is.null(data)) { # internal data
    info <- paste(lic, discl, "\nReference data set :", set, "(internal data)")
  } else {              # CSV or XLS(X)
    if (missing(path.in) |
        regexpr("/library/replicateBE/extdata/", path.in)[1] >= 1) { # internal CSV
      info <- paste(lic, discl, "\nReference data set :", set, "(internal CSV)")
    } else {                                                         # external CSV or XLS(X)
      info <- paste(lic, discl, "\nInput from         : ")
      if (is.null(path.in)) {
        info <- paste0(info, normalizePath(getwd()), winslash = "/")
      } else {
        info <- paste0(info, normalizePath(path.in), winslash = "/")
      }
      if (ext %in% ext.xls) {
        info <- paste0(info, "\nFile [sheet]       : ", file, ".", ext,
                       " [", set, "]")
      } else {
        info <- paste0(info, "\nFile               : ", file, set, ".", ext)
      }
    }
  }
  info <- paste(info, "\nOutput to          : ")
  if (is.null(path.out)) {
    #info <- paste0(info, getwd(), "/")
    info <- paste0(info, normalizePath(getwd(), winslash = "/"))
  } else {
    info <- paste0(info, normalizePath(path.out, winslash = "/"))
  }
  info <- paste0(info, "\nSystem             : ", node,
                 "\nUser               : ", user,
                 "\nOperating System   : ", OS, " ", OSrel)
  if (OS == "Darwin") { # special treatment (long system[["version"]])
    tmp <- strwrap(OSver, width = 79, prefix="\n                     ")
    for (j in 1:length(tmp)) {
      info <- paste0(info, tmp[[j]])
    }
  } else {
    info <- paste0(info, " ", OSver)
  }
  info <- paste0(info, "\nR version          : ",
                 sprintf("%-10s", rver), ryear)
  if (ext %in% ext.xls) {
    info <- paste0(info, "\nreadxl version     : ",
                   sprintf("%-10s", packageVersion("readxl")), year1)
  }
  info <- paste0(info, "\nPowerTOST version  : ",
                 sprintf("%-10s", packageVersion("PowerTOST")), year2)
  if (fun == "method.B") {
    if (option == 2) {
      info <- paste0(info, "\nnlme version       : ",
                     sprintf("%-10s", packageVersion("nlme")), year3)
    } else {
      info <- paste0(info, "\nlmerTest version   : ",
                     sprintf("%-10s", packageVersion("lmerTest")), year4)
    }
  }
  info <- paste0(info, "\nreplicateBE version: ",
                 sprintf("%-10s", packageVersion("replicateBE")), year)
  info <- paste0(info, "\n", hr,
                 "\nFunction           : CV.calc(): stats:lm() executed ", exec,
                 "\n  Fixed effects    : sequence, subject(sequence), period",
                 "\n  Data             : treatment = R")
  info <- paste0(info, "\nFunction           : ", fun, "(")
  if (is.na(option)) {
    info <- paste0(info, "): stats:lm()")
    info <- paste0(info, " executed ", exec)
  } else {
    info <- paste0(info, "option=", option, "): ")
    ifelse (option == 2,
            info <- paste0(info, "nlme:lme()\n"),
            info <- paste0(info, "lmerTest:lmer()\n"))
    info <- paste(info, "                    executed", exec)
  }
  if (fun %in% c("ABE", "method.A")) {
    info <- paste0(info, "\n  Fixed effects    : sequence, ",
                   "subject(sequence), period, treatment",
                   "\n  Data             : all")
  } else {
    info <- paste0(info, "\n  Fixed effects    : sequence, period, treatment",
                   "\n  Random effect    : subject(sequence)",
                   "\n  Data             : all")
  }
  info <- paste0(info, "\n", hr, "\n")
  return(info)
} # end of function env.info()

#########################
# function to show the  #
# BE-limits, CI, and PE #
#########################
repBE.draw.line <- function(called.from, L, U, lo, hi, PE, theta1, theta2) {
  # unicode symbols:
  # confidence interval:   filled black square
  #                        ]max.range[ left and/or right triangle
  # point estimate:        white rhombus
  # expanded limits:       double vertical line
  # BE-limits, GMR-restr.: single vertical line
  # spaghetti Viennese to catch all possible combinations
  # the 'resolution' is ca. 0.5%
  # the CI and PE have presedence over the limits
  if(missing(theta1)) theta1 <- 0.8
  if(missing(theta2)) theta2 <- 1/theta1
  s         <- c("\u256B", "\u255F", "\u2562",
                 "\u253C", "\u251C", "\u2524",
                 "\u25CA", "\u25A0", "\u25C4", "\u25BA",
                 "\u2500", "\u00A0")
  names(s) <- c("EX", "EX1", "EX2",
                "BE", "BE1", "BE2",
                "PE", "CI", "CL.lo", "CL.hi",
                "li", "sp")
  sf   <- 107.6     # scaling factor to get a 79 character string
  L.0  <- 0.6983678 # max. lower expansion
  U.0  <- 1.4319102 # max. upper expansion
  repl <- function(l, sf, L.0, loc, sym) {
    substr(l, sf*(loc-L.0)+1, sf*(loc-L.0)+1) <- sym
    return(l)
  }
  l <- paste0(rep(s["sp"], sf*(U.0-L.0)+1), collapse="")
  l <- repl(l, sf, L.0, 1, s["BE"])
  if (!called.from == "ABE") {
    ifelse (L == 0.8  & lo > 0.8, sym <- s["BE1"], sym <- s["BE"])
    l <- repl(l, sf, L.0, 0.80, sym)
    ifelse (U == 1.25 & hi < 1.25, sym <- s["BE2"], sym <- s["BE"])
    l <- repl(l, sf, L.0, 1.25, sym)
    if (L != 0.8) { # scaled
      ifelse (lo >= L, sym <- s["EX1"], sym <- s["EX"])
      l <- repl(l, sf, L.0, L, sym)
    } else {        # unscaled
      ifelse (lo >= 0.8, sym <- s["BE1"], sym <- s["BE"])
      l <- repl(l, sf, L.0, 0.8, sym)
    }
    if (U != 1.25) {
      ifelse (hi <= U, sym <- s["EX2"], sym <- s["EX"])
      l <- repl(l, sf, L.0, U, sym)
    } else {
      ifelse (hi >= 1.25, sym <- s["BE"], sym <- s["BE2"])
      l <- repl(l, sf, L.0, 1.25, sym)
    }
  } else {
    ifelse (lo < theta1, sym <- s["BE"], sym <- s["BE1"])
    l <- repl(l, sf, L.0, theta1, sym)
    ifelse (hi > theta2, sym <- s["BE"], sym <- s["BE2"])
    l <- repl(l, sf, L.0, theta2, sym)
  }
  l <- repl(l, sf, L.0, PE, s["PE"])
  ifelse (lo < L.0, l <- repl(l, sf, L.0, L.0, s["CL.lo"]),
          l <- repl(l, sf, L.0, lo, s["CI"]))
  ifelse (hi > U.0, l <- repl(l, sf, L.0, U.0, s["CL.hi"]),
          l <- repl(l, sf, L.0, hi, s["CI"]))
  last <- sf*(U.0-L.0)+1
  while (last <= sf*(U.0-L.0)+1) { # last non-space character
    last <- last - 1
    if (substr(l, last, last) != s["sp"]) break
  }
  if (substr(l, sf*(U.0-L.0)+1, sf*(U.0-L.0)+1) != s["sp"]) # special case
    last <- sf*(U.0-L.0)+1
  first <- 0
  while (first < last) {          # first non-space character
    first <- first + 1
    if (substr(l, first, first) != s["sp"]) break
  }
  while (first < last) { # replace space with line
    if (substr(l, first, first) == s["sp"])
      substr(l, first, first) <- s["li"]
    first <- first + 1
  }
  return(l) # TODO: trim trailing whitespace. How for unicode string?
}

###########################################
# return information of the internal data #
# containg variables as required by other #
# functions. Don't forget to update for   #
# new reference data sets!                #
############################################
info.data <- function(data = NULL) {
  if (missing(data) | is.null(data)) stop()
  sets     <- 30
#    descr    <- c("Dataset I given by the EMA available at
#  https://www.ema.europa.eu/en/documents/other/31-annex-ii-statistical-analysis-bioequivalence-study-example-data-set_en.pdf",
#                 "Dataset II given by the EMA (Q&A document) available at https://www.ema.europa.eu/en/documents/other/statistical-method-equivalence-studies-annex-iii_en.pdf",
#                 "Modified dataset I given by the EMA: Period\u00A03 removed.",
#                 "Cmax data of Table II from Patterson SD, Jones B. Viewpoint: observations on scaled average bioequivalence. Pharm Stat. 2012:11(1):1\u20137. doi:10.1002/pst.498",
#                 "Cmax data of the Appendix from Metzler CM, Shumaker RC. The Phenytoin Trial is a Case Study of \u2018Individual\u2019 Bioequivalence. Drug Inf J. 1998:32:1063\u201372.",
#                 "Modified dataset I given by the EMA: T and R switched.",
#                 "Dataset simulated with CVwT\u00A0=\u00A0CVwR\u00A0=\u00A035%, GMR\u00A0=\u00A00.90.",
#                 "Dataset simulated with CVwT\u00A0=\u00A070%, CVwR\u00A0=\u00A080%, CVbT\u00A0=\u00A0CVbR\u00A0=\u00A0150%, GMR\u00A0=\u00A00.85.",
#                 "Dataset with wide numeric range (based of rds08: Data of last 37 subjects multiplied by 1,000,000).",
#                 "Table 9.3.3 (AUC) from: Chow SC, Liu JP. Design and Analysis of Bioavailability and Bioequivalence Studies. Boca Raton: CRC Press; 3rd edition 2009. p275.",
#                 "Table 9.6 (Cmax) from: Hauschke D, Steinijans VW, Pigeot I. Bioequivalence Studies in Drug Development. Chichester: John Wiley: 2007. p216. (Drug\u00A017a of the FDA\u2019s bioequivalence study files: available at https://www.fda.gov/downloads/Drugs/ScienceResearch/UCM301481.zip).",
#                 "Dataset simulated with extreme intra- and intersubject variability, GMR\u00A0= 1.6487.",
#                 "Highly incomplete dataset (based of rds08: Approx. 50% of period\u00A04 data deleted).",
#                 "Dataset simulated with extreme intra- and intersubject variability, GMR\u00A0= 1. Dropouts as a hazard function growing with period.",
#                 "Highly incomplete data set (based of rds08: Approx. 50% of period\u00A04 data are coded as missing '.').",
#                 "Drug 14a, Cmax data: MAO inhibitor\u00A0- IR of the FDA\u2019s bioequivalence study files: available at https://wayback.archive-it.org/7993/20170723175533/https://www.fda.gov/downloads/Drugs/ScienceResearch/UCM301481.zip.",
#                 "Highly unbalanced dataset (based on rds03: 12 subjects in RTR and 7 in TRT).",
#                 "Highly incomplete dataset (based on rds14: T data of subjects 63\u201378 removed).",
#                 "Highly incomplete dataset (based on rds18: Data of subjects 63\u201378 removed).",
#                 "Highly incomplete dataset (based on rds19: Outlier of R (subject\u00A01) introduced: original value \u00D7100).",
#                 "Modified dataset I given by the EMA: One extreme result of subjects 45
# & 52 set to NA.",
#                 "Dataset simulated with CVwT\u00A0= CVwR\u00A0=\u00A045%, CVbT\u00A0= CVbR\u00A0=\u00A0100%, GMR\u00A0=\u00A00.90.",
#                 "Drug 7a, Cmax data: Beta-adrenergic blocking agent\u00A0- IR of the FDA\u2019s bioequivalence study files: available at https://www.fda.gov/downloads/Drugs/ScienceResearch/UCM301481.zip.",
#                 "Drug 1, Cmax data: Antianxiety agent\u00A0- IR of the FDA\u2019s bioequivalence study files: available at https://www.fda.gov/downloads/Drugs/ScienceResearch/UCM301481.zip.",
#                 "Dataset simulated with CVwT\u00A0=\u00A050%, CVwR\u00A0=\u00A080%, CVbT\u00A0=\u00A0CVbR\u00A0=\u00A0130%, GMR\u00A0=\u00A00.90.",
#                 "Example 4.4 (Cmax) from: Patterson SD, Jones B. Bioequivalence and Statistics in Clinical Pharmacology. Boca Raton: CRC Press; 2nd edition 2016. p105\u20136.",
#                 "Dataset simulated with CVwT\u00A0= CVwR\u00A0=\u00A035%, CVbT = CVbR\u00A0=\u00A075%, GMR\u00A0=\u00A00.90.",
#                 "Dataset simulated with CVwT\u00A0= CVwR\u00A0=\u00A035%, CVbT = CVbR\u00A0=\u00A075%, GMR\u00A0=\u00A00.90.",
#                 "Highly imbalanced & incomplete dataset simulated with CVwT\u00A0=\u00A014%, CVwR\u00A0=\u00A028%, CVbT\u00A0=\u00A028% CVbR\u00A0=\u00A056%, GMR\u00A0=\u00A00.90.",
#                 "Highly imbalanced & incomplete dataset simulated with CVwT\u00A0=\u00A014%, CVwR\u00A0=\u00A028%, CVbT\u00A0=\u00A028% CVbR\u00A0=\u00A056%, GMR\u00A0=\u00A00.90.")
  file     <- rep("DS", sets)
  set      <- sprintf("%02i", 1:sets)
  ref      <- paste0("rds", set)
  id       <- data.frame(file, set, ref, 
                         # descr, 
                         stringsAsFactors = FALSE)
  act      <- attr(data, "rset")
  if (!(act %in% ref)) {
    info <- NULL
  } else {
    info <- id[as.integer(substr(act, 4, 5)), ]
  }
  return(info)
}

#########################################
# Sequences of known and tested designs #
# according to the preferred order.     #
# Sort sequences of unknown design (T   #
# first) and throw a message that the   #
# design is untested.                   #
#########################################
info.design <- function(seqs = NA) {
  sequences <- length(seqs)
  if (sequences < 2)
    stop("At least 2 sequences required.")
  if (!is.character(seqs))
    stop("Sequences must be given as strings, not numbers.")
  if (sequences != length(unique(seqs)))
    stop(paste("The", sequences,"sequences must be unique."))
  periods <- unique(nchar(seqs))
  if (periods < 2)
    stop("Not a crossover design.")
  if (periods == 2 & sequences == 2) {
    stop("Not a replicate design.")
  }
  if (length(periods) > 1)
    stop("Each sequence must have the same number of periods.")
  reordered <- NA
  
  # 4-period 4-sequence full replicates
  if (periods == 4 & sequences == 4 & is.na(reordered[1])) {
    if (sum(seqs %in% c("RTRT", "RTTR", "TRRT", "TRTR")) == 4) {
      reordered <- seqs[order(match(seqs, c("TRTR", "RTRT", "TRRT", "RTTR")))]
    }
    if (sum(seqs %in% c("RRTT", "RTTR", "TRRT", "TTRR")) == 4 & is.na(reordered[1])) {
      reordered <- seqs[order(match(seqs, c("TRRT", "RTTR", "TTRR", "RRTT")))]
    }
    if (is.na(reordered[1])) {
      message("Untested design.")
      reordered <- rev(seqs) # at least T first
    }
    design <- "full"
  }
  
  # 4-period 2-sequence full replicates
  if (periods == 4 & sequences == 2 & is.na(reordered[1])) {
    if (sum(seqs %in% c("RTRT", "TRTR")) == 2 & is.na(reordered[1])) {
      reordered <- seqs[order(match(seqs, c("TRTR", "RTRT")))]
    }
    if (sum(seqs %in% c("RTTR", "TRRT")) == 2 & is.na(reordered[1])) {
      reordered <- seqs[order(match(seqs, c("TRRT", "RTTR")))]
    }
    if (sum(seqs %in% c("TTRR", "RRTT")) == 2 & is.na(reordered[1])) {
      reordered <- seqs[order(match(seqs, c("TTRR", "RRTT")))]
    }
    if (is.na(reordered[1])) {
      message("Untested design.")
      reordered <- rev(seqs) # at least T first
    }
    design <- "full"
  }
  
  # 3-period 2-sequence full replicates or partial replicate
  if (periods == 3 & sequences == 2 & is.na(reordered[1])) {
    if (sum(seqs %in% c("RTR", "TRT")) == 2 & is.na(reordered[1])) { # full
      reordered <- seqs[order(match(seqs, c("TRT", "RTR")))]
      design <- "full"
    }
    if (sum(seqs %in% c("RTT", "TRR")) == 2 & is.na(reordered[1])) { # full
      reordered <- seqs[order(match(seqs, c("TRR", "RTT")))]
      design <- "full"
    }
    if (sum(seqs %in% c("RTR", "TRR")) == 2 & is.na(reordered[1])) { # extra-reference
      reordered <- seqs[order(match(seqs, c("TRR", "RTR")))]
      design <- "partial"
    }
    if (is.na(reordered[1])) {
      message("Untested design.")
      reordered <- rev(seqs) # at least T first
      design <- "partial" # likely...
    }
  }
  
  # 3-period 3-sequence partial replicate
  if (periods == 3 & sequences == 3 & is.na(reordered[1])) {
    if (sum(seqs %in% c("RRT", "RTR", "TRR")) == 3 & is.na(reordered[1])) {
      reordered <- seqs[order(match(seqs, c("TRR", "RTR", "RRT")))]
    }
    if (is.na(reordered[1])) {
      message("Untested design.")
      reordered <- rev(seqs) # at least T first
    }
    design <- "partial"
  }
  
  # 2-period 4-sequence (Balaam's)
  if (periods == 2 & sequences == 4 & is.na(reordered[1])) {
    if (sum(seqs %in% c("RR", "RT", "TR", "TT")) == 4 & is.na(reordered[1])) {
      reordered <- seqs[order(match(seqs, c("TR", "RT", "TT", "RR")))]
    }
    if (is.na(reordered[1])) {
      message("Untested design.")
      reordered <- rev(seqs) # at least T first
    }
    design <- "full"
  }
  design <- list(reordered, paste0(reordered, collapse="|"),
                 sequences, periods, design)
  names(design) <- c("reordered", "type", "sequences", "periods", "design")
  return(design)
}

#################################
# Get the data from a file or   #
# internal data and generate    #
# output common to all methods. #
#################################
get.data <- function(path.in, path.out, file, set = "",
                     ext, na = ".", sep = ",", dec = ".",
                     logtrans = TRUE, print, plot.bxp, data) {
  transf <- logtrans # default
  if (is.null(data)) { # checking external data
    if (!missing(ext)) ext <- tolower(ext) # case-insensitive
    ext.csv <- "csv"
    ext.xls <- c("xls", "xlsx")
    if (missing(file))
      stop("Argument 'file' must be given.")
    if (is.numeric(file))
      stop("Argument 'file' must be a string (i.e., enclosed in single or double quotes).")
    if (is.numeric(set))
      stop("Argument 'set' must be a string (i.e., enclosed in single or double quotes).")
    if (missing(ext))
      stop("Argument 'ext' (file-extension) must be given.")
    if (!ext %in% c(ext.csv, ext.xls))
      stop("Data format not supported (must be Character Separated Variables or Excel.")
    if (ext %in% ext.csv) {
      if (!sep %in% c(";", ",", "\t"))
        stop(paste0("Reading CSV-file\n       Argument 'sep' (variable separator) must be any of",
                    "\n       ',' (comma = default), ';' (semicolon), or '\\t' (tab)."))
      if (!dec %in% c(".", ","))
        stop("Reading CSV-file: Argument 'dec' (decimal separator) must be\n'.' (period = default) or ',' (comma).")
      if (sep == dec)
        stop("Reading CSV-file\n       Arguments 'sep' and decimal 'dec' must be different.")
    }
    if (ext %in% ext.xls & (set == ""))
      stop("Reading Excel\n       Argument 'set' (name of worksheet) must be given.")
    if (is.null(path.in) | missing(path.in)) {
      stop("Argument 'path.in' not given. Please specify one or use '~/' for your home folder.")
    }
    if (!dir.exists(path.in)) {
      stop("Folder given in 'path.in' does not exist; please specify an existing one.")
    }
    path.in <- normalizePath(path.in, winslash = "/")
    # Adds trailing '/' to path if missing
    path.in <- ifelse(regmatches(path.in, regexpr(".$", path.in)) == "/",
                      path.in, paste0(path.in, "/"))
  } # EO checking external data
  if (print | plot.bxp) { # check only if necessary
    if (missing(path.out)) {
      stop("Argument 'path.out' not given. Please specify one or use '~/' for your home folder.")
    }
    if (is.null(path.out)) {
      stop("Argument 'path.out' not given. Please specify one or use '~/' for your home folder.")
    }
    if (!dir.exists(path.out)) {
      stop("Folder given in 'path.out' does not exist; please specify an existing one.")
    }
    path.out <- normalizePath(path.out, winslash = "/")
    # Adds trailing '/' to path if missing
    path.out <- ifelse(regmatches(path.out, regexpr(".$", path.out)) == "/",
                       path.out, paste0(path.out, "/"))
  } # EO print/plot checks
  if (is.null(data)) {      # read data from file
    if (ext %in% ext.csv) full.name <- paste0(path.in, file, set, ".", ext)
    if (ext %in% ext.xls) full.name <- paste0(path.in, file, ".", ext)
    if (!file.exists(full.name)) {
      setwd(dirname(file.choose()))
      path.in <- paste0(getwd(), "/")
      full.name <- paste0(path.in, file, ".", ext)
    }
    # Read the entire content
    if (ext %in% ext.xls) { # read from Excel to the data frame
      datawithdescr <- as.data.frame(read_excel(path=full.name, sheet=set,
                                                na=c("NA", "ND", ".", "", "Missing"),
                                                skip=0, col_names=FALSE, .name_repair="minimal"))
    } else {
      datawithdescr <- read.csv(file=full.name, sep=sep, dec=dec, quote="'\"",
                                header=FALSE, strip.white=TRUE,
                                na.strings=c("NA", "ND", ".", "", "Missing"),
                                stringsAsFactors=FALSE)
    }
    namesvector <- c("subject", "period", "sequence", "treatment")
    # Looking for a row with namesvector, summing all its members and
    # if all are there, mark as TRUE
    Nnamesdf <- c(t(apply(datawithdescr, 1, function(row, table) {
      sum(match(tolower(row), table=table), na.rm=TRUE)}, table=namesvector)) == 10)
    if (sum(Nnamesdf) == 0)
      stop("Column names must be given as 'subject', 'period', 'sequence', 'treatment'.")
    if (sum(Nnamesdf) > 1) {
      err.msg <- paste("More than 1 row with column names 'subject', 'period'",
                       "\n       'sequence', 'treatment' detected.")
      stop(err.msg)
    }
    # If there are some # comments in datafiles, collapse them
    # if (which(Nnamesdf == TRUE)-1) {
    #   # Selecting rows before names of dataset
    #   if (!ext %in% ext.xls) {
    #     descr <- scan(file=full.name, what=character(), quiet = TRUE, sep = "\n",
    #                   nlines = (which(Nnamesdf == TRUE)-1))
    #     descr <- paste0(descr[startsWith(descr, "#")], collapse="\n")
    #     descr <- gsub("#", "", descr)
    #   } else {
    #     descrdf <- datawithdescr[1:(which(Nnamesdf == TRUE)-1), ]
    #     descr <- unname(apply(descrdf, 1, function(row){paste0(row[!is.na(row)], collapse=" ")}))
    #   }
    #   descr <- paste0(strwrap(descr, width = 78), collapse="\n")
    # } else {
    #   descr <- ""
    # }
    data <- datawithdescr[(which(Nnamesdf == TRUE)+1):nrow(datawithdescr), ]
    names(data) <- lapply(datawithdescr[which(Nnamesdf == TRUE), ], as.character)
    # Convert eventual mixed or upper case variable names to lower case
    facs  <- which(!names(data) %in% c("PK", "logPK")) # will be factors later
    names(data)[facs] <- tolower(names(data)[facs])
    # if the file contains a first column named NA (imported row.name):
    if (typeof(data[[1]]) == "integer") data <- data[, -1] # remove it
    # from demo(error.catching)
    tryCatch.W.E <- function(expr) {
      W <- NULL
      w.handler <- function(w) { # warning handler
        W <<- w
        invokeRestart("muffleWarning")
      }
      list(value=withCallingHandlers(tryCatch(expr, error=function(e) e),
                                     warning=w.handler), warning=W)
    }
    PKcols <- which(names(data) %in% c("PK", "logPK")) # PK columns
    for (j in seq_along(PKcols)) {                     # transform to numeric
      msg1 <- tryCatch.W.E(data[, PKcols[j]] <- as.numeric(data[, PKcols[j]]))$warning
      invisible(grepl("NAs introduced by coercion", msg1))
      if (!is.null(msg1)) { # NAs as specified in "na" or complete data
        if (!grepl("NAs introduced by coercion", msg1)) {
          msg2 <- paste0("Import: Missing data according to your specifier na='", na, "'")
          msg2 <- paste0(msg2, " not found\n  in column ", names(data)[PKcols[j]], ".")
          msg2 <- paste0(msg2, " Any other non-numeric value was converted to NA.")
          warning(msg2)
        } else {
          print(msg1) # Problem
        }
      }
    }
    # data[, PKcols] <- lapply(data[, PKcols], as.numeric) # throws warnings
    # EO reading file
    if (print) res.file <- paste0(path.out, file, set, "_ABEL")
    if (plot.bxp) png.path <- paste0(path.out, file, set, "_boxplot.png")
    # If the user erroneously asks for analysis of logPK - which does not
    # exist in the data set - change to internal log-transformation.
    if (logtrans == FALSE & !"logPK" %in% colnames(data)) {
      warn.msg <- paste0("Requested analysis of already transformed data ('logtrans = FALSE')\n",
                         "  not possible since column 'logPK' does not exist in the dataset.\n",
                         "  Analysis of log-transformed column 'PK' instead.")
      warning(warn.msg)
      logtrans <- TRUE
      transf <- TRUE
    } # factorize variables except response(s)
    cols       <- c("subject", "period", "sequence", "treatment")
    data[cols] <- lapply(data[cols], factor)
    if (sum(!unique(data$treatment) %in% c("R", "T")) !=0)
      stop("treatments must be given as 'R' and 'T'.")
  } else { # EO reading external data
    if (missing(ext)) { # get information of internal data set
      info  <- info.data(data)
      file  <- info$file
      set   <- info$set
      ref   <- info$ref
      # descr <- info$descr
      if (logtrans == FALSE & !"logPK" %in% colnames(data)) {
        warn.msg <- paste0("Requested analysis of already transformed data ('logtrans = FALSE')\n",
                           "  not possible since column 'logPK' does not exist in the dataset.\n",
                           "  Analysis of log-transformed column 'PK' instead.")
        warning(warn.msg)
        logtrans <- TRUE
        transf <- TRUE
      }
    }
  }
  if (print) res.file <- paste0(path.out, file, set, "_ABEL")
  if (plot.bxp) png.path <- paste0(path.out, file, set, "_boxplot.png")
  subjs  <- unique(data$subject)          # Subjects
  seqs   <- levels(unique(data$sequence)) # Sequences
  design <- info.design(seqs=seqs)        # fetch info
  seqs   <- design$reordered              # preferred reordered sequences (T first)
  Npers  <- design$periods                # Number of periods
  Nseqs  <- design$sequences              # Number of sequences
  type   <- design$type                   # Nice identifier string
  design <- design$design                 # "full" or "partial"
  if (nchar(type) == 19) {  # 4-period 4-sequence full replicate designs
    if (Npers != 4) stop("4 periods required in this full replicate design.")
    if (Nseqs != 4) stop("4 sequences required in this full replicate design.")
  }
  if (nchar(type) == 9) {  # 4-period full replicate designs
    if (Npers != 4) stop("4 periods required in this full replicate design.")
    if (Nseqs != 2) stop("2 sequences required in this full replicate design.")
  }
  if (nchar(type) == 7) {  # 3-period replicates
    if (type %in% c("TRT|RTR", "TRR|RTT")) {
      if (Npers != 3) stop("3 periods required in this full replicate design.")
      if (Nseqs != 2) stop("2 sequences required in this full replicate design.")
    }
    if (type == "TRR|RTR") {
      if (Npers != 3) stop("3 periods required in the extra-reference design.")
      if (Nseqs != 2) stop("2 sequences required in the extra-reference design.")
    }
  }
  if (nchar(type) == 11) { # 3-sequence partial replicate or Balaam's design
    if (!type == "TR|RT|TT|RR") {
      if (Npers != 3) stop("3 periods required in this partial replicate design.")
      if (Nseqs != 3) stop("3 sequences required in this partial replicate design.")
    } else {
      if (Npers != 2) stop("2 periods required in Balaam's design.")
      if (Nseqs != 4) stop("4 sequences required in Balaam's design.")
    }
  }
  # next line introduced for DS24 where all data of subject 16 are NA
  # Given that: Do we need na.action=na.omit in lme() any more?
  data <- na.omit(data)
  Nsub.seq <- table(data$sequence[!duplicated(data$subject)])
  # adapt to reordered sequences
  Nsub.seq <- Nsub.seq[order(match(names(Nsub.seq), seqs))]
  ref   <- data[data$treatment == "R", ]
  test  <- data[data$treatment == "T", ]
  if (type == "TR|RT|TT|RR") {
    T.subj <- test[(test$sequence == "RT" | test$sequence == "TR"), "subject"]
    R.subj <- ref[(ref$sequence == "RT" | ref$sequence == "TR"), "subject"]
    NTR <- length(unique(T.subj) %in% unique(R.subj))
  } else {
    NTR <- length(unique(test$subject) %in% unique(ref$subject)) #  >= 1 T & >= 1 R
  }
  n     <- length(unique(data$subject))
  # Get the wide data frame subject+sequence\period
  compl <- reshape(data[c("subject", "sequence", "period", "PK")],
                   idvar=c("subject", "sequence"), timevar=c("period"),
                   direction="wide")
  uncompletedata <- compl[!complete.cases(compl), ] # exclude complete
  if (nrow(uncompletedata)) {
    # Select NAs in period columns
    Miss.per <- data.frame(sapply(compl, function(y) sum(is.na(y))))
    Miss.per <- t(Miss.per[!rownames(Miss.per) %in% c("subject", "sequence"), ])
    colnames(Miss.per) <- paste0("PK.", 1:ncol(Miss.per))
    # Reformat uncomplete
    uncompletedeshaped <- reshape(uncompletedata, direction="long")
    # Exclude periods with data
    deshapeduncomplete <- uncompletedeshaped[!complete.cases(uncompletedeshaped), ]
    Miss.seq <- table(deshapeduncomplete$sequence)
  } else {
    Miss.seq <- rep(0, Nseqs)
    names(Miss.seq) <- seqs
    Miss.per <- rep(0, Npers)
  }
  Miss.seq <- Miss.seq[order(match(names(Miss.seq), seqs))]
  # Data of subjects with two R treatments
  RR   <- ref[duplicated(ref$subject, fromLast=TRUE)|
                duplicated(ref$subject, fromLast=FALSE), ]
  RR   <- RR[!is.na(RR$PK), ]        # exclude NAs
  nRR  <- length(unique(RR$subject)) # number of subjects
  nTT  <- NA
  if (design == "full") { # only full replicates
    # Data of subjects with two T treatments
    TT  <- test[duplicated(test$subject, fromLast=TRUE)|
                  duplicated(test$subject, fromLast=FALSE), ]
    TT  <- TT[!is.na(TT$PK), ]        # exclude NAs
    nTT <- length(unique(TT$subject)) # number of subjects
  }
  # if (!any(is.na(descr)) && !is.null(descr) && length(descr) >= 1)
  #   txt <- paste0(strwrap(descr, width = 78), collapse="\n")
  txt = ''
  if (logtrans) {
    txt <- paste0(txt,
                  "\nAnalysis performed on column \u2018PK\u2019 ",
                  "(data internally log-transformed)")
  } else {
    txt <- paste0(txt,
                  "\nAnalysis performed on column \u2018logPK\u2019 ",
                  "(data already log-transformed)")
  }
  txt <- paste0(txt, "\nSequences (design) : ", type)
  if (nchar(type) == 19) txt <- paste(txt, "(4-period 4-sequence full replicate)")
  if (nchar(type) == 9) txt <- paste(txt, "(4-period full replicate)")
  if (type %in% c("TRT|RTR", "TRR|RTT")) txt <- paste(txt, "(3-period full replicate)")
  if (type == "TR|RT|TT|RR") txt <- paste(txt, "(Balaam\u2019s 2-period 4-sequence replicate)")
  if (type == "TRR|RTR|RRT") txt <- paste(txt, "(partial replicate)")
  if (type == "TRR|RTR") txt <- paste(txt, "(partial replicate; extra-reference)")
  x <- paste0(Nsub.seq, collapse = "|")
  if (sum(Miss.seq) > 0) x <- c(x, paste0(Miss.seq, collapse = "|"))
  if (sum(Miss.per) > 0) x <- c(x, paste0(Miss.per, collapse = "|"))
  x.len <- nchar(x)
  x.max <- max(x.len)
  txt <- paste0(txt, "\nSubjects / sequence: ", x[1])
  if (length(unique(Nsub.seq)) == 1) {
    txt <- paste0(txt, paste0(rep(" ", x.max-x.len[1]+1), collapse=""),
                  "(balanced)")
  } else {
    txt <- paste0(txt, paste0(rep(" ", x.max-x.len[1]+1), collapse=""),
                  "(unbalanced)")
  }
  if (sum(Miss.seq) > 0) {
    txt <- paste0(txt, "\nMissings / sequence: ", paste0(x[2],
                                                         paste0(rep(" ", x.max-x.len[2]+1), collapse=""),
                                                         "(incomplete)"))
  }
  if (sum(Miss.per) > 0) {
    txt <- paste0(txt, "\nMissings / period  : ", paste0(x[3],
                                                         paste0(rep(" ", x.max-x.len[3]+1), collapse=""),
                                                         "(incomplete)"))
  }
  txt <- paste0(txt,
                "\nSubjects (total)   : ", sprintf("%3i", n),
                "\nSubj\u2019s with T and R: ", sprintf("%3i", NTR),
                " (calculation of the CI)")
  if (NTR < 12) {
    txt <- paste0(txt, "\n                     ",
                  "Less than 12 as required acc. to the BE-GL.")
  }
  if (design == "full") {
    txt <- paste0(txt, "\nSubj\u2019s with two Ts : ", sprintf("%3i", nTT))
  }
  txt <- paste0(txt, "\nSubj\u2019s with two Rs : ", sprintf("%3i", nRR))
  if ((type == "TRT|RTR" | type == "TRR|RTT") & nRR < 12) {
    txt <- paste(txt, "(uncertain CVwR acc. to Q&A Rev. 12)")
  }
  ret <- list(data=data, ref=ref, RR=RR, test=test, type=type, n=n,
              nTT=ifelse(design == "full", nTT, NA), nRR=nRR,
              design=design, 
              txt=txt, 
              Sub.Seq=Nsub.seq,
              Miss.seq=Miss.seq, Miss.per=Miss.per, logtrans=transf,
              res.file=ifelse(print, res.file, NA),
              png.path=ifelse(plot.bxp, png.path, NA))
  return(ret)
} # end of function get.data()

