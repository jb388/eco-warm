library(SoilR)
library(dplyr)
library(readxl)
library(forecast)
library(ISRaD)
library(ggplot2)

# Atm 14C
#####
# use Hua2021 (SoilR package)
NHZone2 <- read_excel(
  "/Users/jeff/sra-frc/data/external/Hua_2021/S0033822221000953sup002.xls", sheet = 2, skip = 5,
  col_names = c("Year.AD", "mean.Delta14C", "sd.Delta14C", "mean.F14C", "sd.F14C")) %>%
  data.frame

# forecast data through 2025
yrs <- seq(2000, 2019.25, by = 1/4) # final sequence of years by quarters
nz2 <- spline(NHZone2[ , c(1, 4)], xout = yrs) # quarterly spline interpolation of fm data
nhz2 <- ts(nz2$y, start = 2000, freq = 4) # Transformation into a time-series object
m <- ets(nhz2) # Fits an exponential smoothing state space model to the time series
f2 <- forecast(m, h = 5.5 * 4) # Uses the fitted model to forecast  years into the future

# bind pre and post-bomb curves; add forecasted data
atm14c <- rbind(
  bind.C14curves(IntCal20, NHZone2, "AD"),
  data.frame(Year.AD = seq(tsp(f2$mean)[1], tsp(f2$mean)[2],
                           by = 1 / tsp(f2$mean)[3]),
             Delta14C = suppressWarnings(
               convert_fm_d14c(
                 fm = as.numeric(f2$mean), 
                 obs_date_y = seq(2019.75, 2022, by = .25), 
                 verbose = FALSE)),
             Sigma = NA)
)

# filter to 1900-2022 and calc annual averages
Datm <- data.frame(Date = seq(1900.5, 2024.5), d14c = NA)
for (i in seq_along(Datm$Date)) {
  ix <- which(atm14c$Year.AD >= Datm[i, "Date"] & atm14c$Year.AD < Datm[i, "Date"] + 1)
  Datm[i, "d14c"] <- mean(atm14c[ix, "Delta14C"], na.rm = TRUE)
}
atm14c.1991 <- Datm[Datm$Date > 1991, ] # HF
atm14c.2006 <- Datm[Datm$Date > 2006, ] # B4W
atm14c.2013 <- Datm[Datm$Date > 2013, ] # Blgt
atm14c.2015 <- Datm[Datm$Date > 2015, ] # SPRUCE
atm14c.2016 <- Datm[Datm$Date > 2016, ] # SWELTR

# steady-state stock calc fx
calc.soc <- function(A, in_vector) {
  (-1 * solve(A) %*% in_vector)
}

# lambda (true half-life of 14C)
lambda <- 1 / 8267

# get initial fm value from k (pre-bomb)
fm <- function (k) {
  k / (k + lambda)
}

# mod fun 2pp
mod.fun.2pp <- function(atm14c, in_vector, pars, cstock, F0_Delta14C, scenario) {
  
  mod <- TwopParallelModel14(
    t = atm14c[,1],
    ks = pars[1:2],
    C0 = cstock[,1],
    F0_Delta14C = F0_Delta14C,
    In = In,
    gam = gam,
    inputFc = atm14c,
    lambda = lambda,
    lag = 0,
    pass = FALSE) 

  # ctl scenario
  C14m <- getF14C(mod)
  C14p <- getF14(mod) 
  C14r <- getF14R(mod)
  Ctot <- getC(mod)
  
  c14 <- data.frame(
    years = rep(atm14c$Date, (ncol(C14p) + 3)),
    d14C = c(c(C14p),
             C14m,
             C14r,
             atm14c$d14c),
    pool = rep(c("fast", "slow", "bulkC", "respiration", "atm"), 
               each = nrow(C14p)),
    scenario = scenario) %>% 
    distinct
  list(c14, Ctot)
}


# ctl
## pars
In <- 1
k1.2pp <- .2
k2.2pp <- .01
gam <- .9
pars <- c(k1.2pp, k2.2pp, gam)

# calculate initial 14C
F0_Delta14C <- Delta14C_from_AbsoluteFractionModern(fm(pars[1:2]))

# 2pp steady-state C stocks
in_vector <- c(In * pars[3], In * (1 - pars[3]))

# 2pp mod matrix
A <- -diag(pars[1:2])

# define initial cstock
ctl.ss.cstock <- calc.soc(A, in_vector)

## run mod
ctl.ls <- mod.fun.2pp(Datm, in_vector, pars, ctl.ss.cstock, F0_Delta14C, "ctl")

# warm, Q10 = 2
## pars
In <- 1
k1.2pp <- .2 * 2^(4/10)
k2.2pp <- .01 * 2^(4/10)
gam <- .9
pars <- c(k1.2pp, k2.2pp, gam)

# input vector
in_vector <- c(In * pars[3], In * (1 - pars[3]))

# t0 14C, 1991
F0_Delta14C <- ctl.ls[[1]] %>%
  filter(years == 1991.5) %>%
  filter(pool == "fast" | pool == "slow") %>%
  select(d14C) %>%
  unlist

## run mod
wrm.s20.f20.df <- mod.fun.2pp(atm14c.1991, in_vector, pars, ctl.ss.cstock, F0_Delta14C, "wrm, s = 2.0, f = 2.0")

# warm, slow Q10 = 1; fast Q10 = 2
## pars
In <- 1
k1.2pp <- .2 * 2^(4/10)
k2.2pp <- .01
gam <- .9
pars <- c(k1.2pp, k2.2pp, gam)

## run mod
wrm.s10.f20.df <- mod.fun.2pp(atm14c.1991, in_vector, pars, ctl.ss.cstock, F0_Delta14C, "wrm, s = 1.0, f = 2.0")

# warm, Q10 = 2, increased inputs
## pars
In <- 2
k1.2pp <- .2 * 2^(4/10)
k2.2pp <- .01 * 2^(4/10)
gam <- .9
pars <- c(k1.2pp, k2.2pp, gam)

# input vector
in_vector <- c(In * pars[3], In * (1 - pars[3]))

## run mod
wrm.s20.f20.2xI.df <- mod.fun.2pp(atm14c.1991, in_vector, pars, ctl.ss.cstock, F0_Delta14C, "wrm, s = 2.0, f = 2.0")

# warm, slow Q10 = 1; fast Q10 = 2
## pars
In <- 2
k1.2pp <- .2 * 2^(4/10)
k2.2pp <- .01
gam <- .9
pars <- c(k1.2pp, k2.2pp, gam)

# input vector
in_vector <- c(In * pars[3], In * (1 - pars[3]))

## run mod
wrm.s10.f20.2xI.df <- mod.fun.2pp(atm14c.1991, in_vector, pars, ctl.ss.cstock, F0_Delta14C, "wrm, s = 1.0, f = 2.0")


# plot
## control
ctl.p <- ctl.ls[[1]] %>% 
  mutate(scenario = factor(scenario, 
                           levels = c("ctl", 
                                      "wrm, s = 1.0, f = 2.0",
                                      "wrm, s = 2.0, f = 2.0"))) %>%
  ggplot(., aes(years, d14C, color = pool, linetype = scenario)) +
  geom_line() +
  scale_color_manual(
    values = c("atm" = "gray", 
               "bulkC" = "black",
               "respiration" = "#8549c3", #6c36a2, #8549c3
               "fast" = "#d8006c",
               "slow" = "#006cd8")) +
  scale_linetype(drop = FALSE) +
  scale_x_continuous(limits = c(1991, 2025)) +
  scale_y_continuous(limits = c(-20, 200)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16))
ggsave("/Users/jeff/eco-warm/src/plots/1991/ctl.2pp.png", plot = ctl.p)

wrm.s1.f2.p <- ctl.ls[[1]] %>% 
  rbind(wrm.s10.f20.df[[1]]) %>%
  mutate(scenario = factor(scenario, 
                           levels = c("ctl", 
                                      "wrm, s = 1.0, f = 2.0",
                                      "wrm, s = 2.0, f = 2.0"))) %>%
  ggplot(., aes(years, d14C, color = pool, linetype = scenario)) +
  geom_line() +
  scale_color_manual(
    values = c("atm" = "gray", 
               "bulkC" = "black",
               "respiration" = "#8549c3", #6c36a2, #8549c3
               "fast" = "#d8006c",
               "slow" = "#006cd8")) +
  scale_linetype(drop = FALSE) +
  scale_x_continuous(limits = c(1991, 2025)) +
  scale_y_continuous(limits = c(-20, 200)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16))
ggsave("/Users/jeff/eco-warm/src/plots/1991/wrm.s1.f2.2pp.png", plot = wrm.s1.f2.p)

wrm.s2.f2.p <- ctl.ls[[1]] %>% 
  rbind(wrm.s10.f20.df[[1]], wrm.s20.f20.df[[1]]) %>%
  mutate(scenario = factor(scenario, 
                           levels = c("ctl", 
                                      "wrm, s = 1.0, f = 2.0",
                                      "wrm, s = 2.0, f = 2.0"))) %>%
  ggplot(., aes(years, d14C, color = pool, linetype = scenario)) +
  geom_line() +
  scale_color_manual(
    values = c("atm" = "gray", 
               "bulkC" = "black",
               "respiration" = "#8549c3", #6c36a2, #8549c3
               "fast" = "#d8006c",
               "slow" = "#006cd8")) +
  scale_linetype(drop = FALSE) +
  scale_x_continuous(limits = c(1991, 2025)) +
  scale_y_continuous(limits = c(-20, 200)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16))
ggsave("/Users/jeff/eco-warm/src/plots/1991/wrm.s2.f2.2pp.png", plot = wrm.s2.f2.p)


# t0 14C, 2016
F0_Delta14C <- ctl.ls[[1]] %>%
  filter(years == 2016.5) %>%
  filter(pool == "fast" | pool == "slow") %>%
  select(d14C) %>%
  unlist

## pars
In <- 1
k1.2pp <- .2 * 2^(4/10)
k2.2pp <- .01 * 2^(4/10)
gam <- .9
pars <- c(k1.2pp, k2.2pp, gam)

## run mod
wrm.s20.f20.16.df <- mod.fun.2pp(atm14c.2016, in_vector, pars, ctl.ss.cstock, F0_Delta14C, "wrm, s = 2.0, f = 2.0")

# warm, slow Q10 = 1; fast Q10 = 2
## pars
In <- 1
k1.2pp <-  .2 * 2^(4/10) # modified Van't Hoff: k_ref * Q10^((Temp - T_ref)/10)
k2.2pp <- .01
gam <- .9
pars <- c(k1.2pp, k2.2pp, gam)

## run mod
wrm.s10.f20.16.df <- mod.fun.2pp(atm14c.2016, in_vector, pars, ctl.ss.cstock, F0_Delta14C, "wrm, s = 1.0, f = 2.0")

wrm.s20.f20.16.p <- ctl.ls[[1]] %>% 
  rbind(wrm.s10.f20.16.df[[1]], wrm.s20.f20.16.df[[1]]) %>%
  mutate(scenario = factor(scenario, 
                           levels = c("ctl", 
                                      "wrm, s = 1.0, f = 2.0",
                                      "wrm, s = 2.0, f = 2.0"))) %>%
  ggplot(., aes(years, d14C, color = pool, linetype = scenario)) +
  geom_line() +
  scale_color_manual(
    values = c("atm" = "gray", 
               "bulkC" = "black",
               "respiration" = "#8549c3", #6c36a2, #8549c3
               "fast" = "#d8006c",
               "slow" = "#006cd8")) +
  scale_linetype(drop = FALSE) +
  scale_x_continuous(limits = c(2016, 2025)) +
  scale_y_continuous(limits = c(-20, 100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16))
ggsave("/Users/jeff/eco-warm/src/plots/2016/wrm.s2.f2.16.2pp.png", plot = wrm.s20.f20.16.p)

wrm.s20.f20.2xI.p <- ctl.ls[[1]] %>% 
  rbind(wrm.s10.f20.2xI.df[[1]], wrm.s20.f20.2xI.df[[1]]) %>%
  mutate(scenario = factor(scenario, 
                           levels = c("ctl", 
                                      "wrm, s = 1.0, f = 2.0",
                                      "wrm, s = 2.0, f = 2.0"))) %>%
  ggplot(., aes(years, d14C, color = pool, linetype = scenario)) +
  geom_line() +
  scale_color_manual(
    values = c("atm" = "gray", 
               "bulkC" = "black",
               "respiration" = "#8549c3", #6c36a2, #8549c3
               "fast" = "#d8006c",
               "slow" = "#006cd8")) +
  scale_linetype(drop = FALSE) +
  scale_x_continuous(limits = c(1991, 2025)) +
  scale_y_continuous(limits = c(-20, 200)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16))
ggsave("/Users/jeff/eco-warm/src/plots/1991/wrm.s2.f2.2xI.2pp.png", plot = wrm.s20.f20.2xI.p)
