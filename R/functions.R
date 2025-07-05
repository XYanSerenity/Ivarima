#' Title
#'
#' @param key  The key values of the evaluated models
#' @param xdata Intervention variables
#' @param xdatar Intervention impact model parameter r
#' @param xdatab  Intervention impact model parameter b
#' @param xdatah  Intervention impact model parameter h
#' @param UnEvaledlist  List of unassessed models
#' @return Updated list of unassessed models.
#' @export
AddunEvaled <- function(key, xdata, xdatar, xdatab, xdatah, UnEvaledlist) {
  if (!key %in% c(names(UnEvaledlist))) {
    UnEvaledlist[[key]] <- list(
      Xdata = xdata,
      Xdatar = xdatar,
      Xdatab = xdatab,
      Xdatah = xdatah
    )
  }
  return(UnEvaledlist)
}


#' Find the intermediate position corresponding to the maximum positive or negative change in the v-weights of the LTF model
#' @param x  LTF model v-weights
#' @return  the intermediate position corresponding to the maximum positive or negative change
#' @export
find_middle_position <- function(x) {
  max_pos <- which(x == max(x[x > 0]))
  min_neg <- which(x == min(x[x < 0]))
  if (length(max_pos) > 0 && length(min_neg) > 0) {
    if (abs(max_pos - min_neg) == 2) {
      return(seq(min(max_pos, min_neg) + 1, max(max_pos, min_neg) - 1))
    } else {
      return(0)
    }
  } else {
    return(0)
  }
}

#' Find  the position of  the maximum positive or negative change in v-weights of the LTF model
#' @param coefs  LTF model v-weights
#' @return the position of the maximum positive or negative change
#' @export
#'
find_max_sign_diff <- function(coefs) {
  max_diff <- 0
  max_diff_index <- NULL
  for (i in 1:(length(coefs) - 1)) {
    if (coefs[i] * coefs[i + 1] < 0) {
      diff <- abs(coefs[i] - coefs[i + 1])
      if (diff > max_diff) {
        max_diff <- diff
        max_diff_index <- i
      }
    }
  }
  if (is.null(max_diff_index)) {
    return(0)
  } else {
    return(max_diff_index)
  }
}

#' Find the peak value of the v-weights in the LTF model
#' @param model the fitted LTF model
#' @param num Number of LTF model v-weights
#' @return  the position of the peak value
#' @export
find_max_v <- function(model, num) {
  coef <- lmtest::coeftest(model)
  vcoef <- tail(coef, num)
  df_v <- data.frame(vcoef[!is.na(vcoef[, 1]), ])
  colnames(df_v) <- c("Coefficient", "Std_Error", "Z_Value", "P_Value")
  peak <- which.max(df_v$Coefficient)
  Ptpeak <- df_v$Coefficient[peak]
  return(Ptpeak)
}
#' Building the LTF model
#' @param y  a univariate time series
#' @param x_p a pulse function interrupted by an intervention event
#' @param rlist  The ARIMA model order (p, d, q) of noise component
#' @param LTFmaxk The longest time-lagged response k in the LTF model
#' @return a fitted LTF model
#' @export
LTFxpFun <- function(y, x_p, rlist, LTFmaxk) {
  lagged_xp <- embed(x_p, LTFmaxk)
  xpvn <- LTFmaxk - 1
  xp_names <- paste0("xpv", 0:xpvn)
  colnames(lagged_xp) <- xp_names
  y_clean <- tail(y, nrow(lagged_xp))
  TFpt <- try(arima(y_clean, order = rlist, xreg = lagged_xp, method = "ML"))
  if (inherits(TFpt, "try-error")) {
    return(NULL)
  } else {
    return(TFpt)
  }
}

#' Using the identification rules of the LTF model, automatically identify possible sets of b
#' @param modelxp  a fitted LTF model
#' @return  the possible sets of b
#' @export
##AAA
IdenPatternb <-function(modelxp)
{
  coefpt <- lmtest::coeftest(modelxp)
  ptvcoef <- tail(coefpt, 11)  # 提取最后11个系数
  df_ptv <- data.frame(ptvcoef[!is.na(ptvcoef[, 1]), ])
  colnames(df_ptv) <- c("Coefficient", "Std_Error", "Z_Value", "P_Value")
  # print(df_ptv)
  ##调试时的画图
  barplot(tail(coef(modelxp), 11),las=2)
  # step1,找p最小点
  Ptminp <- which.max(abs(df_ptv$Coefficient))
  # print("Ptminp")
  #print(Ptminp)
  peak <- which.max(abs(df_ptv$Coefficient))
  Ptpeak <-df_ptv$Coefficient[peak]
  dmaxvsign <- find_max_sign_diff(df_ptv$Coefficient)
  if(Ptminp==1)
  {
    return(list(b = 0))
  }else{
    condition <- all(
      (abs(df_ptv$Coefficient[1:(Ptminp - 1)] / Ptpeak) < 0.3)
    )
    if (condition)
    {
      return(list(b = Ptminp - 1))
    }else{
      numb <- 0
      if(Ptminp>1){
        for (pi in 1:(Ptminp- 1)) {
          if (df_ptv$P_Value[pi] > 0.05 | abs(df_ptv$Coefficient[pi] / Ptpeak) < 0.3) {
            numb <- numb + 1
          } else {
            break
          }
        }}
      print(numb)
      Ptminp1 <- Ptminp - 1
      if(!is.null(dmaxvsign))
      {
        if(numb == dmaxvsign)
        {
          numb <- numb
        }else
        {
          numb <- c(numb,dmaxvsign)
        }
      }
      print(numb)
      ##diffmav
      dmaxv <-find_middle_position(df_ptv$Coefficient)
      if(dmaxv!=0){
        if(!dmaxv%in%numb)
        {
          numb <-c(numb,dmaxv)
        }
      }
      if(Ptminp1>0){
        if(df_ptv$Coefficient[Ptminp1] / Ptpeak< 0.3)
        {
          if (!Ptminp1%in%numb){
            return(list(b = c(numb,Ptminp1)))
          }else {
            return(list(b = numb))
          }
        }else{
          if (!Ptminp1%in%numb){
            numb = c(numb,Ptminp1)
          }
        }}
      if(Ptminp > 2 )
      {
        if(abs(df_ptv$Coefficient[Ptminp - 2] / Ptpeak)< 0.3)
        {
          Ptminp2<-Ptminp - 2
          if(!Ptminp2 %in%numb){
            numb <-c(numb,Ptminp2)
          }
        }
      }
      return(list(b = numb))
    }
  }
}
#' Find the optimal model with the smallest AIC from the evaluated models
#' @param Evaled the evaluated models
#' @return The optimal model's name and its AIC value.
#' @export
#'
FindMinAic <- function(Evaled) {
  min_aic <- Inf
  min_ele <- NULL
  for (name in names(Evaled)) {
    if (!is.null(Evaled[[name]]$aic)) {
      cur_aic <- Evaled[[name]]$aic
      if (cur_aic < min_aic) {
        coefz <- Evaled[[name]]$result
        colnames(coefz) <- c("Coefficient", "Std_Error", "Z_Value", "P_Value")
        if ("Pt-AR1" %in% row.names(coefz) && "Pt-MA1" %in% row.names(coefz)) {
          PtAR1p <- coefz["Pt-AR1", "P_Value"]
          PtAR1f <- coefz["Pt-AR1", "Coefficient"]
          PtMA1p <- coefz["Pt-MA1", "P_Value"]
          PtMA1f <- coefz["Pt-MA1", "Coefficient"]
          PtMA0f <- coefz["Pt-MA0", "Coefficient"]
          if (is.nan(PtAR1p) || is.nan(PtMA1p)) { } else if (PtAR1p < 0.05 && PtAR1f > 0 && PtAR1f < 1 && PtMA1p < 0.05 && (PtMA1f * PtMA0f > 0)) {
            min_aic <- cur_aic
            min_ele <- name
          } else {
            min_aic <- min_aic
            min_ele <- min_ele
          }
        } else if ("Pt-AR1" %in% row.names(coefz)) {
          PtAR1p <- coefz["Pt-AR1", "P_Value"]
          PtAR1f <- coefz["Pt-AR1", "Coefficient"]
          if (is.nan(PtAR1p)) { } else if (PtAR1p < 0.05 && PtAR1f > 0 && PtAR1f < 1) {
            min_aic <- cur_aic
            min_ele <- name
          } else {
            min_aic <- min_aic
            min_ele <- min_ele
          }
        } else if ("Pt-MA1" %in% row.names(coefz)) {
          PtMA1p <- coefz["Pt-MA1", "P_Value"]
          PtMA1f <- coefz["Pt-MA1", "Coefficient"]
          PtMA0f <- coefz["Pt-MA0", "Coefficient"]
          if (is.nan(PtMA1p)) { } else if (PtMA1p < 0.05 && (PtMA1f * PtMA0f > 0)) {
            min_aic <- cur_aic
            min_ele <- name
          } else {
            min_aic <- min_aic
            min_ele <- min_ele
          }
        } else if ("St-MA1" %in% row.names(coefz)) {
          StMA1p <- coefz["St-MA1", "P_Value"]
          StMA1f <- coefz["St-MA1", "Coefficient"]
          StMA0f <- coefz["St-MA0", "Coefficient"]
          if (is.nan(StMA1p)) { } else if (StMA1p < 0.05 && (StMA0f * StMA1f > 0)) {
            min_aic <- cur_aic
            min_ele <- name
          } else {
            min_aic <- min_aic
            min_ele <- min_ele
          }
        } else {
          min_aic <- cur_aic
          min_ele <- name
        }
      }
    }
  }

  return(list(min_ele = min_ele, aic = min_aic))
}
#' Extract the intervention component and intercept coefficients of the ARIMA model.
#' @param mode  a fitted LTF model
#' @return the coefficients of the intervention component and intercept
Tracoef <- function(mode) {
  coef <- coef(mode)
  coefx <- coef[!grepl("ar|ma", names(coef))]
  return(coefx)
}
#' Fitting an ARIMAx model with known initial values for the intervention component and the intercept term
#'
#' @param y a univariate time series
#' @param order regular ARIMA order
#' @param xtransf intervention variables
#' @param transfer a list consisting of the ARMA orders for each transfer (distributed lag) covariate
#' @param inintra initial values for the intervention component and the intercept term
#'
#' @return The Optarimax function returns the ARIMAX model fitted with known initial values for the intervention component and the intercept term
#' @export
Optarimax <- function(y, order, xtransf, transfer, inintra) {
  pqn <- order[1] + order[3]
  inincoef <- c(rep(NA, pqn), inintra)
  curmod <- TSA::arimax(y, order = order, xtransf = xtransf, transfer = transfer, method = "ML", init = inincoef)
  while (1) {
    cur_init <- c(rep(NA, pqn), Tracoef(curmod))
    temod2 <- TSA::arimax(y,
      order = order, xtransf = xtransf, init = cur_init,
      transfer = transfer, method = "ML"
    )
    if (abs(temod2$loglik - curmod$loglik) < 0.1 || temod2$loglik < curmod$loglik) {
      break
    } else {
      curmod <- temod2
    }
  }
  return(curmod)
}
#'
#' A improved parameters estimation method for TSA::arimax
#' @param y  a univariate time series
#' @param order regular ARIMA order
#' @param xtransf intervention variables
#' @param transfer a list consisting of the ARMA orders for each transfer (distributed lag) covariate.
#' @return The ArimaOptExtr function returns a fitted model with robust parameter estimates.
#' @export
#'
OptimInarimax <- function(y, order, xtransf, transfer) {
  # step1
  armalist <- order
  inittra <- c(0)
  mod <- try(TSA::arimax(y, order = order, xtransf = xtransf, transfer = transfer, method = "ML"))
  if (inherits(mod, "try-error")) {
    curList <- order
    Loop <- TRUE
    while (Loop && curList[1] < 5 && curList[3] < 5) {
      if (curList[1] > 0) {
        curList <- c(curList[1] + 1, curList[-1])
      } else if (curList[3] > 0) {
        curList <- c(curList[-3], curList[3] + 1)
      }
      newmod <- try(TSA::arimax(y, order = curList, xtransf = xtransf, transfer = transfer, method = "ML"))
      if (inherits(newmod, "try-error")) {

      } else {
        Loop <- FALSE
      }
      inittra <- Tracoef(newmod)
    }
  } else {
    inittra <- Tracoef(mod)
  }
  optmod <- Optarimax(y, armalist, xtransf, transfer, inittra)
  if (!inherits(mod, "try-error") && inherits(optmod, "try-error")) {
    return(mod)
  } else if (!inherits(mod, "try-error") && !inherits(optmod, "try-error")) {
    if (optmod$aic < mod$aic) {
      return(optmod)
    } else {
      return(mod)
    }
  } else if (inherits(mod, "try-error")) {
    return(optmod)
  }
}

#' Extract the p, d, and q parameters of the ARIMA model.
#'
#' @param modelAR  A fitted ARIMA model.
#'
#' @return The ArimaOptExtr function returns a vector which is the p,d,q parameters
#' @export
ArimaOptExtr <- function(modelAR) {

  non_seasonal_order <- modelAR$arma[c(1, 6, 2)]
  return(non_seasonal_order)
}

#' An automated, rule-based approach for identifying the optimal intervention impact model in interrupted time series,
#' constrained to cases where the intervention variable is an impulse function, the impact model follows
#' a first-order intervention transfer function, and the noise component is a non-seasonal ARIMA model
#' @param  y is a time series.
#' @param its_start a number between 1 and length(y)-1 stating the time point of the start of the intervention.
#' @param LTFmaxk The LTFmaxk accounts for the longest time-lagged response in Linear Transfer Function  model.
#' @return The function returns the optimal model fitting results found.
#' @export
#' @examples
#' w0 <- 0.5
#' n <- 100
#' xt6 <- c(rep(0, 56), rep(1, 44))
#' e <- arima.sim(model = list(order = c(1, 0, 0), ar = 0.4, ma = NULL), n = 100, sd = 0.1)
#' xxt <- w0 * xt6
#' mu <- 4 + xxt
#' y <- mu + e
#' model <- auto.ivarima(y = y, its_start = 51, LTFmaxk = 10)
#' summary(model)
auto.ivarima <- function(y, its_start, LTFmaxk) {
  LTFmaxk <- LTFmaxk + 1
  UnEvaledlist <- list()
  Evaledlist <- list()
  Evaledmodel <- list()
  rlistp0 <- c(1, 0, 0)
  n <- length(y)
  pren <- its_start - 1
  postn <- n - pren
  # pulse function generation
  x_p <- c(rep(0, pren), 1, rep(0, postn - 1))
  # step  function generation
  x_t <- c(rep(0, pren), rep(1, postn))
  TFptMod <- LTFxpFun(y, x_p, rlistp0, LTFmaxk)
  resloop <- TRUE
  loopi <- 1
  if (!is.null(TFptMod)) {
    while (resloop && loopi < 6) {
      # Perform the Ljung-Box test
      TFactest <- Box.test(TFptMod$residuals, lag = 20, type = "Ljung-Box")
      if (TFactest$p.value < 0.05) {
        autores <- forecast::auto.arima(TFptMod$residuals)
        # Extract p, d, and q parameters
        respdq <- ArimaOptExtr(autores)
        respdq <- as.numeric(respdq)
        if (identical(respdq, c(0, 0, 0))) {
          resloop <- FALSE
        } else {
          TFptMod <- LTFxpFun(y, x_p, respdq, LTFmaxk)
        }
        loopi <- loopi + 1
      } else {
        resloop <- FALSE
      }
    }
    xpresult <- IdenPatternb(TFptMod)
    # possible combinations of b and h
    combhs1 <- expand.grid(b = xpresult$b, h = c(1, 0))
    # possible combinations of b,h and r
    for (i in seq_len(nrow(combhs1))) {
      pb <- combhs1[i, "b"]
      ph <- combhs1[i, "h"]
      pname <- paste0("1", pb, ph, "PulseD")
      UnEvaledlist <- AddunEvaled(pname, data.frame(Pt = x_p), 1, pb, ph, UnEvaledlist) # 更新外部变量
    }
    while (length(UnEvaledlist) > 0) {
      first <- UnEvaledlist[[1]]
      # Extract key
      first_keys <- names(UnEvaledlist)[1]
      # If b > 0, lag the intervention variable.
      if (first$Xdatab > 0) {
        lag <- first$Xdatab
        first$Xdata <- lag(zoo::zoo(first$Xdata), -lag, na.pad = TRUE)
        first$Xdata <- ifelse(is.na(first$Xdata), 0, first$Xdata)
      }
      bestmod <- try(OptimInarimax(y,
        order = c(1, 0, 0), xtransf = first$Xdata,
        transfer = list(c(first$Xdatar, first$Xdatah))
      ))
      if (inherits(bestmod, "try-error")) {
        UnEvaledlist[[1]] <- NULL
      } else {
        coefz <- lmtest::coeftest(bestmod)
        # Delete rows named "intercept" and those starting with "ar" or "ma".
        coefz <- coefz[!rownames(coefz) %in% "intercept" & !grepl("^(ar|ma)", rownames(coefz)), , drop = FALSE]
        if ("Pt-AR1" %in% row.names(coefz)) {
          PtAR1f <- coefz["Pt-AR1", "Estimate"]
          if (PtAR1f > 0.9 || PtAR1f < 0) {
            addnameS <- paste0("0", first$Xdatab, first$Xdatah, "Step")
            UnEvaledlist <- AddunEvaled(addnameS, data.frame(St = x_t), 0, first$Xdatab, first$Xdatah, UnEvaledlist)
          }
          PtAR1p <- coefz["Pt-AR1", "Pr(>|z|)"]
          if (is.nan(PtAR1p) || PtAR1p >= 0.05 || PtAR1f < 0) {
            addnameP <- paste0("0", first$Xdatab, first$Xdatah, "Pulse")
            UnEvaledlist <- AddunEvaled(addnameP, data.frame(Pt = x_p), 0, first$Xdatab, first$Xdatah, UnEvaledlist)
          }
        }
        if ("Pt-MA1" %in% row.names(coefz)) {
          PtMA1p <- coefz["Pt-MA1", "Pr(>|z|)"]
          if (is.nan(PtMA1p) || PtMA1p >= 0.05) {
            strPt <- sub("^[0-9]+", "", first_keys)
            addnamePM <- paste0(first$Xdatar, first$Xdatab, 0, strPt)
            UnEvaledlist <- AddunEvaled(addnamePM, data.frame(Pt = x_p), first$Xdatar, first$Xdatab, 0, UnEvaledlist)
          }
        }
        if ("St-MA1" %in% row.names(coefz)) {
          StMA1p <- coefz["St-MA1", "Pr(>|z|)"]
          if (is.nan(StMA1p) || StMA1p >= 0.05) {
            strSt <- sub("^[0-9]+", "", first_keys)
            addnameSM <- paste0(first$Xdatar, first$Xdatab, 0, strSt)
            UnEvaledlist <- AddunEvaled(addnameSM, data.frame(St = x_t), first$Xdatar, first$Xdatab, 0, UnEvaledlist)
          }
        }
        Evaledmodel[[first_keys]] <- list(model = bestmod)
        Evaledlist[[first_keys]] <- list(
          aic = bestmod$aic,
          result = data.frame(coefz)
        )
        UnEvaledlist[[1]] <- NULL
      }
    }
  }
  if (length(Evaledlist) > 0) {
    findresult <- FindMinAic(Evaledlist)
    optname <- findresult$min_ele
    return(list(optname=optname,model= Evaledmodel[[optname]]$model))
  } else {
    return(NULL)
  }
}
