# Read in the data
# install required R libraries if not installed already
list.of.packages <- c("car", "tidyverse", "psych", "metafor", "esc", "lme4", "ggplot2", "knitr", "puniform", "kableExtra", "lmerTest", "pwr", "Amelia", "multcomp", "magrittr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load required libraries
lapply(list.of.packages, require, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)
select <- dplyr::select

datNature <- read_delim("dataNature.csv", delim = ",", trim_ws = T)
datSocial <- read_delim("dataSocialSupport.csv", delim = ",", trim_ws = T)
dat <- bind_rows(datNature, datSocial)

dat <- dat %>% modify_at(., .at = c("F", "beta", "t", "r", "chiSq", "df1", "df2", "sd1", "sd2", "se1", "se2", "n1", "n2", "mean1", "mean2", "published", "researchDesign", "pubYear", "nMale", "nFemale"), .f = ~as.numeric(as.character(.)))
dat$result <- 1:nrow(dat)

# Some data wrangling to get the right type of data (formatting of the raw dataset in Excel introduces a lot of junk otherwise)
dat$pReported <- as.numeric(as.character(gsub("[^0-9.]", "", dat$pReported)))

# Which designs are present?
table(dat$researchDesign, useNA="ifany")

# Compute gender ratio (% of female)
dat$percFemale <- ifelse(is.na(dat$percFemale), (dat$nFemale/(dat$nFemale + dat$nMale))*100, dat$percFemale)

# Compute SDs from SEs
dat <- dat %>% mutate(sd1 = ifelse(is.na(sd1) & !is.na(se1), se1*sqrt(n1), sd1),
                      sd2 = ifelse(is.na(sd2) & !is.na(se2), se2*sqrt(n2), sd2))

# Assume balanced design when cell sizes not reported
dat <- dat %>% mutate(n1 = ifelse(!is.na(mean1) & !is.na(nTotal) & is.na(n1), nTotal/2, n1),
                      n2 = ifelse(!is.na(mean2) & !is.na(nTotal) & is.na(n2), nTotal/2, n2))

dat <- escalc(measure = "SMD", m1i = mean1, m2i = mean2, sd1i = sd1, sd2i = sd2, n1i = n1, n2i = n2, data = dat, include = researchDesign %in% c(1, 2))
dat[dat$researchDesign == 3,]$yi <- dat %>% filter(researchDesign == 3) %$% escalc(measure = "SMCC", m1i = mean1, m2i = mean2, sd1i = sd1, sd2i = sd2, ni = n1, ri = c(rep(rmCor, nrow(.))))$yi
dat[dat$researchDesign == 3,]$vi <- dat %>% filter(researchDesign == 3) %$% escalc(measure = "SMCC", m1i = mean1, m2i = mean2, sd1i = sd1, sd2i = sd2, ni = n1, ri = c(rep(rmCor, nrow(.))))$vi
dat <- dat %>% mutate(yi = abs(yi) * directionEffect)
dat <- dat %>% rowwise %>% mutate(p = case_when(researchDesign %in% c(1, 2) ~  as.numeric(round(tTestSummary(mean1, mean2, sd1, sd2, n1, n2, withinSS = FALSE)["p-value"], 5)),
                                                researchDesign == 3 ~ as.numeric(round(tTestSummary(mean1, mean2, sd1, sd2, n1, n2, withinSS = TRUE)["p-value"], 5)))) %>% data.frame()

# Initialize new variables
dat$gConv <- NA
dat$gVarConv <- NA
dat$useCellN <- NA

# Create result and study ID
dat$paperID
dat$study <- paste(dat$paperID, "/", dat$studyID, sep = "")

# # F-Test between with df1 == 1 ---------------------------------------------------------------------
# # Specify the design, compute ni and p
# dat <- dat %>% mutate(finalDesign = case_when(!is.na(F) & !is.na(df1) & !is.na(df2) & df1 == 1 ~ "F1"))
# dat <- dat %>% mutate(ni = ifelse(finalDesign == "F1", df2 + 2, NA))
# 
# # Decide whether n1+n2 approximately corresponds to the reported df (in order for n1 and n2 be used in equations).
# dat <- dat %>% mutate(useCellN = ifelse((dat$n1 + dat$n2) >= (dat$ni - 2) & (dat$n1 + dat$n2) <= (dat$ni + 2), 1, 0),
#                       useCellN = ifelse(is.na(useCellN), 0, useCellN))
# 
# # Compute ES and var based on total ni. Note to self: mutate not working due to a glitch in handling NAs by esc.
# dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign),]$gConv <- dat %>% filter(dat$finalDesign == "F1") %$% esc_f(f = F, totaln = df2 + 2, es.type = "g")$es
# dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign),]$gVarConv <- dat %>% filter(dat$finalDesign == "F1") %$% esc_f(f = F, totaln = df2 + 2, es.type = "g")$var
# 
# # Compute ES and var based on n1 and n2 if available. Note to self: mutate not working due to a glitch in handling NAs by esc.
# dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$gConv <- dat %>% filter(dat$finalDesign == "F1" & dat$useCellN == 1) %$% esc_f(f = F, grp1n = n1, grp2n = n2, es.type = "g")$es
# dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$gVarConv <- dat %>% filter(dat$finalDesign == "F1" & dat$useCellN == 1) %$% esc_f(f = F, grp1n = n1, grp2n = n2, es.type = "g")$var
# 
# # Show the converted ESs
# dat %>% filter(finalDesign == "F1") %>% select(gConv, gVarConv, researchDesign, df2, n1, n2, ni, useCellN)

# t-tests between or t for B in continuous designs ---------------------------------------------------------
# Specify the design, compute ni and p
dat <- dat %>% mutate(finalDesign = ifelse(!is.na(t) & !is.na(df2), "tBtw", NA),
                      ni = ifelse(finalDesign == "tBtw", df2 + 2, NA),
                      p = ifelse(finalDesign == "tBtw" & is.na(p), 2*pt(abs(t), df2, lower.tail = FALSE), p))

# Decide whether n1+n2 approximately corresponds to the reported df (in order for n1 and n2 be used in equations).
dat <- dat %>% mutate(useCellN = ifelse((dat$n1 + dat$n2) >= (dat$ni - 2) & (dat$n1 + dat$n2) <= (dat$ni + 2), 1, useCellN),
                      useCellN = ifelse(is.na(useCellN), 0, useCellN))

# Compute ES and var based on total ni. Note to self: mutate not working due to a glitch in handling NAs by esc.
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign),]$gConv <- dat %>% filter(dat$finalDesign == "tBtw") %$% esc_t(t = abs(t), totaln = df2 + 2, es.type = "g")$es
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign),]$gVarConv <- dat %>% filter(dat$finalDesign == "tBtw") %$% esc_t(t = abs(t), totaln = df2 + 2, es.type = "g")$var

# Compute ES and var based on n1 and n2 if available. Note to self: mutate not working due to a glitch in handling NAs by esc.
# dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$gConv <- dat %>% filter(dat$finalDesign == "tBtw" & dat$useCellN == 1) %$% esc_t(t = abs(t), grp1n = n1, grp2n = n2, es.type = "g")$es
# dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$gVarConv <- dat %>% filter(dat$finalDesign == "tBtw" & dat$useCellN == 1) %$% esc_t(t = abs(t), grp1n = n1, grp2n = n2, es.type = "g")$var

# Show the converted ESs
dat %>% filter(finalDesign == "tBtw") %>% select(t, gConv, gVarConv, researchDesign, df2, n1, n2, ni, useCellN)

# Correlation -------------------------------------------------------------
# Specify the design, convert rho to r, and compute ni and p
dat <- dat %>% mutate(finalDesign = ifelse((!is.na(r) | !is.na(r_s)) & (!is.na(df2) | !is.na(nTotal)), "cor", finalDesign),
                      r = ifelse(is.na(r) & !is.na(r_s), abs(2*sin(r_s*pi/6)), abs(r)),
                      ni = ifelse(finalDesign == "cor" & is.na(ni), ifelse(!is.na(df2), df2 + 2, nTotal), ni),
                      p = ifelse(finalDesign == "cor" & is.na(p), 2*pt(abs(r*sqrt(ni - 2 / (1 - r^2))), ni - 2, lower.tail = FALSE), p))

dat <- dat %>% mutate(yi = ifelse(finalDesign == "cor" & is.na(yi), abs(escalc(measure = "COR", ri = r, ni = ni, data = dat)$yi) * directionEffect, yi),
                      vi = ifelse(finalDesign == "cor" & is.na(vi), escalc(measure = "COR", ri = r, ni = ni, data = dat)$vi, vi))

# Show the converted ESs
dat %>% filter(finalDesign == "cor") %>% select(yi, vi, r, r_s, directionEffect, ni, p)

# Betas in between-subjects designs ---------------------------------------

#Betas in between-subjects designs (in ML, beta considered as covariate/confounding adjusted r, then using r to d conversion)
dat <- dat %>% mutate(finalDesign = ifelse(!is.na(dat$beta) & (!is.na(dat$df2) | !is.na(dat$nTotal)), "regression", finalDesign))
dat <- dat %>% mutate(ni = ifelse(finalDesign == "regression", df2 + 2, ni))

# dat[dat$finalDesign == "regression", ] %<>% mutate(
#   t.from.beta = abs(beta)*sqrt(df2 / (1 - beta^2)),
#   gConv = ((2*abs(beta))/sqrt(1-beta^2))*(1 - (3/(4*df2 - 1))),
#   ni = df2 + 2,
#   beta.var = escalc(measure = "COR", ri = beta, ni = df2 + 2, data = dat[dat$finalDesign == "regression", ])$vi,
#   gVarConv = (1 - (3/(4*df2 - 1))) * (4 * beta.var)/((1 - beta^2)^3),
#   p = 2*pt(abs(t.from.beta), df2, lower.tail=FALSE)
# )

# Show the converted ESs
dat %>% filter(finalDesign == "regression") %>% select(gConv, beta, gVarConv, finalDesign, ni)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$Use.for.Bias.Test == "Yes" & (dat$Design == "Between" | dat$Design == "Continuous") & !is.na(dat$beta) & !is.na(dat$p.reported), "biasTest"] <- "Beta.between"
# dat$label[dat$biasTest == "Beta.between"] <- paste(dat[dat$biasTest == "Beta.between",]$paperID, "/", dat[dat$biasTest == "Beta.between",]$Study.Indicator, "/", dat[dat$biasTest == "Beta.between",]$Variable.Indicator, ": ",
#                                                       "Z=", qnorm(1-dat[dat$biasTest == "Beta.between",]$p.reported/2), sep = "")


dat <- dat %>% mutate(ni = ifelse(is.na(ni) & !is.na(yi), n1 + n2, ni),
                      nTerm = 2/ni)

# # Multiply the ES by -1 if in the opposite direction
dat <- dat %>% mutate(yi = ifelse(is.na(yi) & !is.na(gConv) & !is.na(directionEffect), directionEffect * gConv, yi),
                      vi = ifelse(is.na(vi) & !is.na(gVarConv) & !is.na(directionEffect), gVarConv, vi),
                      label = paste(paperID, "/", studyID, "/", effectID, sep = ""),
                      # subgroup analysis A: if !is.na(affect) | !is.na(stressComponentType)) = 1; if !is.na(affectiveConsequencesStress) = 2
                      # subgroup analysis B: if (stressComponentType %in% c(1:4) | !is.na(affect)) ~ 1; stressComponentType = 5 ~ 2; stressComponentType = 6 ~ 3
                      stressAffective = as.factor(ifelse(!is.na(affect) | !is.na(stressComponentType), 1, ifelse(!is.na(affectiveConsequencesStress), 2, NA))),
                      stressCompRecoded = as.factor(case_when(stressComponentType %in% c(1:4) | !is.na(affect) ~ 1,
                                                              stressComponentType == 5 ~ 2,
                                                              stressComponentType == 6 ~ 3)))

# GRIM & GRIMMER Test
grim(dat)
grimmer(dat)

# Subset and create data objects
dat$useMeta <- 1
dat <- dat %>% filter_all(any_vars(!is.na(.)))

dataNature <- dat %>% filter(strategy == 1 & !is.na(yi))
dataSocial <- dat %>% filter(strategy == 2 & !is.na(yi))
dataNatureRoB <- dataNature %>% 
  mutate_at(vars(robDomain1:robOverall), funs(dplyr::recode(., '1' = 'Low', '2' = 'Some concerns', '3' = 'High'))) %>% 
  filter(!is.na(yi)) %>% 
  select(study, robDomain1, robDomain2, robDomain3, robDomain4, robDomain5, robOverall, useMeta)
dataSocialRoB <- dataSocial %>% 
  mutate_at(vars(robDomain1:robOverall), funs(dplyr::recode(., '1' = 'Low', '2' = 'Some concerns', '3' = 'High'))) %>% 
  filter(!is.na(yi)) %>% 
  select(study, robDomain1, robDomain2, robDomain3, robDomain4, robDomain5, robOverall, useMeta)

# Remove outliers (based on the results from the maDiag script)
dataNature <- dataNature %>% filter(!result %in% outliersNature)
dataSocial <- dataSocial %>% filter(!result %in% outliersSocial)

# # Within-subjects design, ES based on t-distribution ----------------------
# 
# # Identify within-subjects design reporting F and compute t
# dat[dat$Use.for.Meta == "Yes" & dat$Design == "Within" & !is.na(dat$F) & !is.na(dat$n.rep), "finalDesign"] <- "within.t"
# dat[dat$Use.for.Meta == "Yes" & dat$Design == "Within" & !is.na(dat$F) & !is.na(dat$n.rep),]$t <- sqrt(dat[dat$Use.for.Meta == "Yes" & dat$Design == "Within" & !is.na(dat$F) & !is.na(dat$n.rep),]$F)
# 
# # Compute ES using t
# dat[dat$Use.for.Meta == "Yes" & dat$Design == "Within" & !is.na(dat$t) & !is.na(dat$n.rep), "finalDesign"] <- "within.t"
# dat[dat$finalDesign == "within.t", ] %<>% mutate(
#   dCalc = abs(t)*sqrt((2 * (1 - rmCor)) / n.rep),
#   gConv = (1 - (3/(4*n.rep - 3))) * dCalc,
#   gVarConv = (1 - (3/(4*n.rep - 3)))^2 * ((1 / n.rep) + ((dCalc^2) / (2 * n.rep))) * 2 * (1 - rmCor),
#   ni = n.rep,
#   p = 2*pt(abs(t), n.rep - 1, lower.tail = FALSE)
# )
# 
# # Show the converted ESs
# dat %>% filter(finalDesign == "within.t") %>% select(gConv, gVarConv, Design, d.reported, ni, df2)
# 
# # Create a "result label" to be used as an input for p-curve analysis
# dat[dat$Use.for.Bias.Test == "Yes" & dat$Design == "Within" & !is.na(dat$t) & !is.na(dat$df2), "biasTest"] <- "within.t"
# dat$label[dat$biasTest == "within.t"] <- paste(dat[dat$biasTest == "within.t",]$paperID, "/", dat[dat$biasTest == "within.t",]$Study.Indicator, "/", dat[dat$biasTest == "within.t",]$Variable.Indicator, ": ",
#                                                     "t(", dat[dat$biasTest == "within.t",]$df2, ")=", dat[dat$biasTest == "within.t",]$t, sep = "")
# 
# # PaperID 71 used mixed-effects models, couldn't convert, so using the reported d (converted to g)
# dat$gConv[266:270] <- (1 - (3/(4*dat$n.rep[266:270] - 3))) * dat$d.reported[266:270]
# dat$gVarConv[266:270] = (1 - (3/(4*dat$n.rep[266:270] - 3))) * ((dat$n.rep[266:270])/(dat$n.rep[266:270]/2 * dat$n.rep[266:270]/2) + (dat$d.reported[266:270]^2)/(2 * (dat$n.rep[266:270])))
# 
# 
# # chi^2 -------------------------------------------------------------------
# 
# # Specify the design, compute ES, var, ni, and p
# dat[dat$Use.for.Meta == "Yes" & dat$Design == "Between" & !is.na(dat$Chisq) & !is.na(dat$n.rep), "finalDesign"] <- "between.chisq"
# 
# dat[dat$finalDesign == "between.chisq", ]$gConv <- esc_chisq(chisq = dat[dat$finalDesign == "between.chisq", ]$Chisq, totaln = dat[dat$finalDesign == "between.chisq", ]$n.rep, es.type = "g")$es
# dat[dat$finalDesign == "between.chisq", ]$gVarConv <- esc_chisq(chisq = dat[dat$finalDesign == "between.chisq", ]$Chisq, totaln = dat[dat$finalDesign == "between.chisq", ]$n.rep, es.type = "g")$var
# dat[dat$finalDesign == "between.chisq", ]$ni <- dat[dat$finalDesign == "between.chisq", ]$n.rep
# dat[dat$finalDesign == "between.chisq", ]$p <- with(dat[dat$finalDesign == "between.chisq", ], 1 - pchisq(Chisq, 1))
# 
# # Show the converted ESs
# dat %>% filter(finalDesign == "between.chisq") %>% select(gConv, gVarConv, Chisq, ni)
# 
# # Create a "result label" to be used as an input for p-curve analysis
# dat[dat$Use.for.Bias.Test == "Yes" & dat$Design == "Between" & !is.na(dat$Chisq) & !is.na(dat$n.rep), "biasTest"] <- "between.chisq"
# dat$label[dat$biasTest == "between.chisq"] <- paste(dat[dat$biasTest == "between.chisq",]$paperID, "/", dat[dat$biasTest == "between.chisq",]$Study.Indicator, "/", dat[dat$biasTest == "between.chisq",]$Variable.Indicator, ": ",
#                                                         "chi2(", dat[dat$biasTest == "between.chisq",]$df1, ")=", dat[dat$biasTest == "between.chisq",]$Chisq, sep = "")
# 
# # t for B in continuous designs -------------------------------------------
# 
# # Specify the design, compute ES, var, ni, and p
# dat[dat$Use.for.Meta == "Yes" & (dat$Design == "Continuous" |  dat$Design == "Correlation") & !is.na(dat$B) & !is.na(dat$t) & !is.na(dat$df2) & is.na(dat$beta), "finalDesign"] <- "t.from.B"
# 
# dat[dat$finalDesign == "t.from.B", ]$gConv <- esc_t(t = abs(dat[dat$finalDesign == "t.from.B", ]$t), totaln = dat[dat$finalDesign == "t.from.B", ]$df2 + 2, es.type = "g")$es
# dat[dat$finalDesign == "t.from.B", ]$gVarConv <- esc_t(t = dat[dat$finalDesign == "t.from.B", ]$t, totaln = dat[dat$finalDesign == "t.from.B", ]$df2 + 2, es.type = "g")$var
# dat[dat$finalDesign == "t.from.B", ]$ni <- dat[dat$finalDesign == "t.from.B", ]$df2 + 2
# dat[dat$finalDesign == "t.from.B", ]$p <- with(dat[dat$finalDesign == "t.from.B", ], 2*pt(abs(t), ni - 1, lower.tail = FALSE))
# 
# # Show the converted ESs
# dat %>% filter(finalDesign == "t.from.B") %>% select(gConv, gVarConv, B, t, ni)
# 
# # Create a "result label" to be used as an input for p-curve analysis
# dat[dat$Use.for.Bias.Test == "Yes" & (dat$Design == "Continuous" |  dat$Design == "Correlation") & !is.na(dat$B) & !is.na(dat$t) & !is.na(dat$df2) & is.na(dat$beta), "biasTest"] <- "t.from.B"
# dat$label[dat$biasTest == "t.from.B"] <- paste(dat[dat$biasTest == "t.from.B",]$paperID, "/", dat[dat$biasTest == "t.from.B",]$Study.Indicator, "/", dat[dat$biasTest == "t.from.B",]$Variable.Indicator, ": ",
#                                                     "t(", dat[dat$biasTest == "t.from.B",]$df2, ")=", dat[dat$biasTest == "t.from.B",]$t, sep = "")
# 
###
# # Create a results label for the rest of bias=yes and meta=no, based on p-values, converted to z-score (p-curve app works with z but not with p).
# dat[dat$Use.for.Bias.Test == "Yes" & is.na(dat$label) & !is.na(dat$paperID), "biasTest"] <- "z.from.p"
# dat$p <- ifelse(is.na(dat$p), dat$p.reported, dat$p)
# dat$label[dat$biasTest == "z.from.p"] <- paste(dat[dat$biasTest == "z.from.p",]$paperID, "/", dat[dat$biasTest == "z.from.p",]$Study.Indicator, "/", dat[dat$biasTest == "z.from.p",]$Variable.Indicator, ": ",
#                                               "Z=", qnorm(1-(dat[dat$biasTest == "z.from.p",]$p)/2), sep = "")
# 
# # Remove unused variables
# delvars <- names(dat) %in% c(
#   "t.from.r", "t.from.beta", "rVar", "beta.var", "dCalc")
# dat <- dat[!delvars]
# rm(delvars)

