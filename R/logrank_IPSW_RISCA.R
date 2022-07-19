#' Calculate median survival time
#'
#' The user enters individual survival data and the weights previously calculated (by using logistic regression for instance).
#' The usual log-rank test is adapted to the corresponding adjusted survival curves.
#' @param time A numeric vector with the follow up time.
#' @param status A numeric vector with the event indicator (0=right censored, 1=event).
#' @param var A numeric vector with the binary var under interest (only two groups).
#' @param weights The weights for correcting the contribution of each individual.
#' By default, the weights are all equaled to 1 and the survival curves correspond to the usual Kaplan-Meier estimator.
#' @export


logrank_IPSW_RISCA <-function (time, status, var, weights = NULL)
{
  if (sum(time < 0) > 0) {
    print("Error time must be positive")
  }  else {
    if (sum(weights <= 0) > 0) {
      print("Error weights must be superior to 0")
    }   else {
      if (sum(status != 0 & status != 1) > 0) {
        print("Error status must be must be a vector of 0 or 1")
      }    else {
        if (is.null(weights)) {
          .w <- rep(1, length(time))
        }   else {
          .w <- weights
        }
        crosssum <- function(seq, point, value) {
          crosssum <- numeric(length(seq))
          for (i in 1:length(point)) {
            loc <- sum(seq <= point[i])
            crosssum[loc] <- crosssum[loc] + value[i]
          }
          crosssum
        }
        time <- round(time, 4)
        n.unit <- rep(1, length(time))
        d.time <- survfit(Surv(time, n.unit) ~ 1)$time
        n.risk0 <- crosssum(d.time, time[var ==
                                            0], n.unit[var == 0])
        n.risk1 <- crosssum(d.time, time[var ==
                                            1], n.unit[var == 1])
        n.risk0 <- rev(cumsum(rev(n.risk0)))
        n.risk1 <- rev(cumsum(rev(n.risk1)))
        n.risk <- n.risk0 + n.risk1
        n.event <- crosssum(d.time, time[status ==
                                            1], status[status == 1])
        mod.rate <- ifelse(n.risk == 1, 0, n.event *
                             (n.risk - n.event)/(n.risk * (n.risk - 1)))
        mod.rate[is.na(mod.rate)] <- 0
        w.risk0 <- crosssum(d.time, time[var ==
                                            0], .w[var == 0])
        w.risk1 <- crosssum(d.time, time[var ==
                                            1], .w[var == 1])
        w.risk0 <- rev(cumsum(rev(w.risk0)))
        w.risk1 <- rev(cumsum(rev(w.risk1)))
        w.risk <- w.risk0 + w.risk1
        w.risk0[w.risk0 == 0] <- 1e-04
        w.risk1[w.risk1 == 0] <- 1e-04
        subset0 <- (status == 1) * (var == 0)
        w.event0 <- crosssum(d.time, time[subset0 ==
                                             1], .w[subset0 == 1])
        subset1 <- (status == 1) * (var == 1)
        w.event1 <- crosssum(d.time, time[subset1 ==
                                             1], .w[subset1 == 1])
        w.event <- w.event0 + w.event1
        ww.risk0 <- crosssum(d.time, time[var ==
                                             0], (.w * .w)[var == 0])
        ww.risk1 <- crosssum(d.time, time[var ==
                                             1], (.w * .w)[var == 1])
        ww.risk0 <- rev(cumsum(rev(ww.risk0)))
        ww.risk1 <- rev(cumsum(rev(ww.risk1)))
        log.mean <- sum(w.event1 - w.risk1 * w.event/w.risk)
        log.var <- sum(mod.rate * (w.risk0 * w.risk0 *
                                     ww.risk1 + w.risk1 * w.risk1 * ww.risk0)/(w.risk *
                                                                                 w.risk))
        v.akm <- log.mean/sqrt(log.var)
        return(list(statistic = v.akm, p.value = 2 *
                      (1 - stats::pnorm(abs(v.akm)))))
      }
    }
  }
}
