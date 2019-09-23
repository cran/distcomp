## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    message = FALSE,
    warning = FALSE,
    error = FALSE,
    tidy = FALSE,
    cache = FALSE
)
suppressMessages(suppressWarnings(library(distcomp)))
dcWorkspace <- tempdir()
ignore <- distcompSetup(workspace = dcWorkspace,
                        ssl_verifyhost = 0L, ssl_verifypeer = 0L)

## ------------------------------------------------------------------------
available <- availableComputations()
( computationNames <- names(available) )

## ------------------------------------------------------------------------
svdDef <- data.frame(compType = computationNames[2],
                     rank = 20L,
                     ncol = 20L,
                     id = "SVD",
                     stringsAsFactors = FALSE)

## ------------------------------------------------------------------------
set.seed(12345)
## Three sites
nSites <- 3
siteData <- lapply(seq.int(nSites), function(i) matrix(rnorm(10000), ncol=20))

## ------------------------------------------------------------------------
sites <- lapply(seq.int(nSites),
                function(x) list(name = paste0("site", x),
                                 worker = makeWorker(defn = svdDef, data = siteData[[x]])
                                 ))

## ------------------------------------------------------------------------
master <- makeMaster(svdDef)
for (site in sites) {
  master$addSite(name = site$name, worker = site$worker)
}

## ------------------------------------------------------------------------
system.time(result <- master$run(max.iter = 10000))

## ------------------------------------------------------------------------
full_data <- do.call(rbind, siteData)
full_svd <- svd(full_data)

## ------------------------------------------------------------------------
d_table <- data.frame(truth = full_svd$d, distcomp = result$d)
knitr::kable(d_table)

## ------------------------------------------------------------------------
norm(abs(result$v) - abs(full_svd$v), "F")

## ------------------------------------------------------------------------
coxDef <- data.frame(compType = computationNames[1],
                     formula = "Surv(time, censor) ~ age + becktota + ndrugfp1 + ndrugfp2 + ivhx3 + race + treat",
                     projectName = "STCoxTest",
                     projectDesc = "STCox Project Desc",
                     stringsAsFactors = FALSE)

## ------------------------------------------------------------------------
## Two sites
siteDataFiles <- file.path(system.file("ex", package="distcomp"), c("uis-site1.csv", "uis-site2.csv"))
siteData <- lapply(siteDataFiles, read.csv)

## ------------------------------------------------------------------------
sites <- lapply(seq_along(siteData),
                function(x) list(name = paste0("site", x),
                                 worker = makeWorker(defn = coxDef, data = siteData[[x]])
                                 ))

## ------------------------------------------------------------------------
master <- makeMaster(coxDef)
for (site in sites) {
  master$addSite(name = site$name, worker = site$worker)
}

## ------------------------------------------------------------------------
result <- master$run()
print(master$summary(), digits = 4)

## ------------------------------------------------------------------------
coxFit <- survival::coxph(formula=Surv(time, censor) ~ age + becktota + ndrugfp1 + ndrugfp2 +
                              ivhx3 + race + treat + strata(site),
                          data = rbind(siteData[[1]], siteData[[2]]))
summary(coxFit)$coefficients

## ------------------------------------------------------------------------
sessionInfo()

