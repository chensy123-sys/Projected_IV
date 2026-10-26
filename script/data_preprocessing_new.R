## Australian Longitudinal Survey Dataverse
## https://dataverse.ada.edu.au/dataverse/australian-longitudinal-survey
rm(list = ls())
library(tidyverse)
library(mice)
library(FactoMineR)
library(missMDA)


###################
## Time period 0 ##
###################

## Australian Longitudinal Survey, 1984: Wave 1, Level 1
## McRae, Ian; Parkinson, Geoff; Kronenberg, Naomi, 2019, "Australian Longitudinal Survey, 1984: Wave 1, Level 1",
## doi:10.26193/HE3BAO, ADA Dataverse, V3

data1984p <- read_csv("results/data_private/00377-p.csv")

# IV
a1 <- data1984p$VO1
a2 <- data1984p$VO2
a3 <- data1984p$VO3
a4 <- data1984p$VO4
a5 <- data1984p$VO5
a6 <- data1984p$VO6
a7 <- data1984p$VO7

x <- cbind(a1, a2, a3, a4, a5, a6, a7)
all.na.idx <- which(rowSums(is.na(x)) == ncol(x))
x <- x[-all.na.idx, ]
data1984p <- data1984p[-all.na.idx, ]

set.seed(123)
x <- imputePCA(x)$completeObs
attitude <- rowSums(x[, -3])
# attitude <- as.numeric(attitude > median(attitude))

# covariates
born_australia <- data1984p$VA12
born_australia <- as.numeric(born_australia == 1) # Born in Australia

capital <- as.numeric(data1984p$VB4 == 1) # Mostly live in capital city or not ?

married <- data1984p$VA9
married <- as.numeric(married == 2) # Married

uni_mem <- data1984p$VG10
uni_mem <- as.numeric(uni_mem == 1) # Union Member

gov_emp <- data1984p$VG9
gov_emp <- as.numeric(gov_emp == 1) # Government Employee

age <- data1984p$VA4 # Age

sex <- data1984p$VA3
sex <- as.numeric(sex == 2) # Female

## years of experience
f3 <- data1984p$VF3 # FULL OR PART TIME STUDY
f3[is.na(f3)] <- 0

f4 <- data1984p$VF4 # HAD FT JOB SINCE LEAVE SCHOOL
f4[is.na(f4)] <- 0

f7 <- data1984p$VF7 # STILL HAVE 1ST FT JOB?
f7[is.na(f7)] <- 0

f8 <- data1984p$VF8 # HOW LONG HAD 1ST FT JOB?
f8[is.na(f8)] <- 0

f9 <- data1984p$VF9 # HOW MANY WKS IF < 1 YR

f10 <- data1984p$VF10 # HOW MANY YRS, IF 1 YR+

f31 <- data1984p$VF31 # TIME HELD JOB

f32 <- data1984p$VF32 # NO. WKS HELD JOB IF < 1 YR

f33 <- data1984p$VF33 # NO. YRS HELD JOB IF 1 YR+

g21 <- data1984p$VG21 # TIME WORKED IN MAIN JOB

g22 <- data1984p$VG22 # WEEKS WORKED IF < 1 YR

g23 <- data1984p$VG23 # YEARS WORKED IF 1 YR+

year1.1 <- as.numeric(g21 == 1) * g22 / 52
year1.2 <- as.numeric(g21 == 2) * g23
year1.1[is.na(year1.1)] <- 0
year1.2[is.na(year1.2)] <- 0

year2.1 <- as.numeric(f3 == 0 & f4 == 1 & f7 == 2 & f8 == 1) * f9 / 52
year2.2 <- as.numeric(f3 == 0 & f4 == 1 & f7 == 2 & f8 == 2) * f10
year2.1[is.na(year2.1)] <- 0
year2.2[is.na(year2.2)] <- 0

year3.1 <- as.numeric(f31 == 1) * f32 / 52
year3.2 <- as.numeric(f31 == 2) * f33
year3.1[is.na(year3.1)] <- 0
year3.2[is.na(year3.2)] <- 0

year_expe <- year1.1 + year1.2 + year2.1 + year2.2 + year3.1 + year3.2

# treatment
e4 <- data1984p$VE4 # STATE LAST ATTENDED SCHOOL

e7 <- data1984p$VE7 # GRADE WHEN LEFT SCHOOL

e10 <- data1984p$VE10 # YEARS OF SCHOOLING
e10[is.na(e10)] <- 0

e14 <- data1984p$VE14 # ANY POST SCHOOL QUALS?

e16 <- data1984p$VE16 # TYPE OF FIRST QUALIFICATION

e23 <- data1984p$VE23 # OBTAINED OTHER QUALIFICATION SINCE

e25 <- data1984p$VE25 # TYPE OF HIGHEST OTHER QUAL.

year1 <- as.numeric(e7 == 1) * 12 + 
  as.numeric(e7 == 2) * 11 + 
  as.numeric(e7 == 3) * 10 + 
  as.numeric(e7 == 4) * 9 + 
  as.numeric(e7 == 5) * 8 + 
  as.numeric(e7 == 6) * 7
year1[is.na(year1)] <- 0

year2 <- as.numeric(e14 == 1) * (as.numeric(e16 == 1) * 4 +
                                   as.numeric(e16 == 2) * 2 +
                                   as.numeric(e16 == 3) * 4 +
                                   as.numeric(e16 == 4) * 3 +
                                   as.numeric(e16 == 5) * 2 +
                                   as.numeric(e16 == 6) * 2 +
                                   as.numeric(e16 == 7) * 1 +
                                   as.numeric(e16 == 8) * 1)
year2[is.na(year2)] <- 0

year3 <- as.numeric(e23 == 1) * (as.numeric(e25 == 1) * 4 +
                                   as.numeric(e25 == 2) * 2 +
                                   as.numeric(e25 == 3) * 4 +
                                   as.numeric(e25 == 4) * 3 +
                                   as.numeric(e25 == 5) * 2 +
                                   as.numeric(e25 == 6) * 2 +
                                   as.numeric(e25 == 7) * 1 +
                                   as.numeric(e25 == 8) * 1)
year3[is.na(year3)] <- 0

year_edu <- pmax(year1 + year2 + year3, e10)
year_edu[year_edu == 0] <- NA


# outcome
g3 <- data1984p$VG3 # MORE THAN 1 JOB CURRENTLY

g4 <- data1984p$VG4 # WEEKLY HRS WORKED MAIN JOB

g5 <- data1984p$VG5 # WEEKLY HRS WORKED TOTAL

g7 <- data1984p$VG7 # TAKE-HOME PAY
g7[g7 == 99] <- NA

g8 <- data1984p$VG8 # GROSS PAY IN MAIN JOB
g8[g8 == 99] <- NA

wage1 <- as.numeric(g7 == 1) * 0 +
  as.numeric(g7 == 2) * 20 +
  as.numeric(g7 == 3) * 60 +
  as.numeric(g7 == 4) * 100 +
  as.numeric(g7 == 5) * 130 +
  as.numeric(g7 == 6) * 150 +
  as.numeric(g7 == 7) * 170 +
  as.numeric(g7 == 8) * 190 +
  as.numeric(g7 == 9) * 210 +
  as.numeric(g7 == 10) * 230 +
  as.numeric(g7 == 11) * 250 +
  as.numeric(g7 == 12) * 270 +
  as.numeric(g7 == 13) * 300 +
  as.numeric(g7 == 14) * 330 +
  as.numeric(g7 == 15) * 350
wage1[is.na(wage1)] <- 0

wage2 <- as.numeric(g8 == 1) * 0 +
  as.numeric(g8 == 2) * 20 +
  as.numeric(g8 == 3) * 60 +
  as.numeric(g8 == 4) * 100 +
  as.numeric(g8 == 5) * 130 +
  as.numeric(g8 == 6) * 150 +
  as.numeric(g8 == 7) * 170 +
  as.numeric(g8 == 8) * 190 +
  as.numeric(g8 == 9) * 210 +
  as.numeric(g8 == 10) * 230
wage2[is.na(wage2)] <- 0

wage <- wage1 + wage2
wage[wage == 0] <- NA

hour1 <- as.numeric(g3 == 1) * (as.numeric(g4 == 1) * 35 +
                                  as.numeric(g4 == 2) * 32 +
                                  as.numeric(g4 == 3) * 27 +
                                  as.numeric(g4 == 4) * 22 +
                                  as.numeric(g4 == 5) * 17 +
                                  as.numeric(g4 == 6) * 12 +
                                  as.numeric(g4 == 7) * 7 +
                                  as.numeric(g4 == 8) * 2.5)
hour1[is.na(hour1)] <- 0

hour2 <- as.numeric(g3 == 2) * (as.numeric(g5 == 1) * 35 +
                                  as.numeric(g5 == 2) * 32 +
                                  as.numeric(g5 == 3) * 27 +
                                  as.numeric(g5 == 4) * 22 +
                                  as.numeric(g5 == 5) * 17 +
                                  as.numeric(g5 == 6) * 12 +
                                  as.numeric(g5 == 7) * 7 +
                                  as.numeric(g5 == 8) * 2.5)
hour2[is.na(hour2)] <- 0
hour <- hour1 + hour2
hour[hour == 0] <- NA

wage_hour <- wage / hour


DATA0 <- data.frame(t = 0, # time
                    wage_hour = wage_hour, # outcome
                    year_edu = year_edu, # treatment
                    attitude = attitude, # IV
                    sex = sex,
                    capital = capital,
                    born_australia = born_australia,
                    married = married,
                    uni_mem = uni_mem,
                    gov_emp = gov_emp,
                    age = age,
                    year_expe = year_expe)

DATA0$born_australia <- as.factor(DATA0$born_australia)
DATA0$married <- as.factor(DATA0$married)
DATA0$uni_mem <- as.factor(DATA0$uni_mem)
DATA0$gov_emp <- as.factor(DATA0$gov_emp)
DATA0 <- DATA0 %>% tidyr::drop_na(wage_hour)



###################
## Time period 1 ##
###################

## Australian Longitudinal Survey, 1985: Wave 1, Level 2
## McRae, Ian; Parkinson, Geoff; Woyzbun, Lyn, 2019, "Australian Longitudinal Survey, 1985: Wave 1, Level 2", 
## doi:10.26193/HMBXJV, ADA Dataverse, V3

data1985p1 <- read_csv("results/data_private/00413-p.csv")
# IV
a1 <- data1985p1$VO1
a1[a1 == 8] <- NA
a2 <- data1985p1$VO2
a2[a2 == 8] <- NA
a3 <- data1985p1$VO3
a3[a3 == 8] <- NA
a4 <- data1985p1$VO4
a4[a4 == 8] <- NA
a5 <- data1985p1$VO5
a5[a5 == 8] <- NA
a6 <- data1985p1$VO6
a6[a6 == 8] <- NA
a7 <- data1985p1$VO7
a7[a7 == 8] <- NA

x <- cbind(a1, a2, a3, a4, a5, a6, a7)
all.na.idx <- which(rowSums(is.na(x)) == ncol(x))
x <- x[-all.na.idx, ]
data1985p1 <- data1985p1[-all.na.idx, ]

set.seed(2025)
x <- imputePCA(x)$completeObs
attitude <- rowSums(x[, -3])
# attitude <- as.numeric(attitude > median(attitude))

# covariates
born_australia <- data1985p1$VB3
born_australia <- as.numeric(born_australia == 1) # Born in Australia

capital <- as.numeric(data1985p1$VB7 == 1) # Mostly live in capital city or not ?

married <- data1985p1$VA7
married <- as.numeric(married == 2) # Married

uni_mem <- data1985p1$VG11
uni_mem <- as.numeric(uni_mem == 1) # Union Member

gov_emp <- data1985p1$VG10
gov_emp <- as.numeric(gov_emp == 1) # Government Employee

age <- data1985p1$VA4 # Age

sex <- data1985p1$VA3
sex <- as.numeric(sex == 2) # Female

## years of experience
f3 <- data1985p1$VF3 # FULL OR PART TIME STUDY
f3[is.na(f3)] <- 0
f3[f3 == 8] <- 0

f4 <- data1985p1$VF4 # HAD FT JOB SINCE LEAVE SCHOOL
f4[is.na(f4)] <- 0

f7 <- data1985p1$VF7 # STILL HAVE 1ST FT JOB?
f7[is.na(f7)] <- 0

f8 <- data1985p1$VF8 # HOW LONG HAD 1ST FT JOB?
f8[is.na(f8)] <- 0
f8[f8 == 8] <- 0

f9 <- data1985p1$VF9 # HOW MANY WKS IF < 1 YR
f9[is.na(f9)] <- 0

f10 <- data1985p1$VF10 # HOW MANY YRS, IF 1 YR+
f10[is.na(f10)] <- 0

f31 <- data1985p1$VF31 # TIME HELD JOB
f31[is.na(f31)] <- 0

f32 <- data1985p1$VF32 # NO. WKS HELD JOB IF < 1 YR
f32[is.na(f32)] <- 0

f33 <- data1985p1$VF33 # NO. YRS HELD JOB IF 1 YR+
f33[is.na(f33)] <- 0

g23 <- data1985p1$VG23 # TIME WORKED IN MAIN JOB
g23[is.na(g23)] <- 0
g23[g23 == 8] <- 0

g24 <- data1985p1$VG24 # WEEKS WORKED IF < 1 YR
g24[is.na(g24)] <- 0

g25 <- data1985p1$VG25 # YEARS WORKED IF 1 YR+
g25[is.na(g25)] <- 0

year_expe <- as.numeric(g23 == 1) * g24 / 52 + as.numeric(g23 == 2) * g25 +
  as.numeric(f3 == 0 & f4 == 1 & f7 == 2) * (as.numeric(f8 == 1) * f9 / 52 + as.numeric(f8 == 2) * f10) +
  as.numeric(f31 == 1) * f32 / 52 + as.numeric(f31 == 2) * f33


# treatment
e3 <- data1985p1$VE3 # STATE LAST ATTENDED SCHOOL

e5 <- data1985p1$VE5 # GRADE WHEN LEFT SCHOOL

e8 <- data1985p1$VE8 # YEARS OF SCHOOLING
e8[is.na(e8)] <- 0

e12 <- data1985p1$VE12 # ANY POST SCHOOL QUALS?

e14 <- data1985p1$VE14 # TYPE OF FIRST QUALIFICATION

e21 <- data1985p1$VE21 # OBTAINED OTHER QUALIFICATION SINCE

e23 <- data1985p1$VE23 # TYPE OF HIGHEST OTHER QUAL.

year1 <- as.numeric(e5 == 1) * 12 + 
  as.numeric(e5 == 2) * 11 + 
  as.numeric(e5 == 3) * 10 + 
  as.numeric(e5 == 4) * 9 + 
  as.numeric(e5 == 5) * 8
year1[is.na(year1)] <- 0

year2 <- as.numeric(e12 == 1) * (as.numeric(e14 == 10) * 4 +
                                   as.numeric(e14 == 20) * 2 +
                                   as.numeric(e14 == 30) * 4 +
                                   as.numeric(e14 == 40) * 3 +
                                   as.numeric(e14 == 50) * 2 +
                                   as.numeric(e14 == 60) * 2 +
                                   as.numeric(e14 == 70) * 1 +
                                   as.numeric(e14 == 80) * 1)
year2[is.na(year2)] <- 0

year3 <- as.numeric(e21 == 1) * (as.numeric(e23 == 10) * 4 +
                                   as.numeric(e23 == 20) * 2 +
                                   as.numeric(e23 == 30) * 4 +
                                   as.numeric(e23 == 40) * 3 +
                                   as.numeric(e23 == 50) * 2 +
                                   as.numeric(e23 == 60) * 2 +
                                   as.numeric(e23 == 70) * 1 +
                                   as.numeric(e23 == 80) * 1)
year3[is.na(year3)] <- 0

year_edu <- pmax(year1 + year2 + year3, e8)
year_edu[year_edu == 0] <- NA


# outcome
g3 <- data1985p1$VG3 # MORE THAN 1 JOB CURRENTLY
g3[is.na(g3)] <- 0

g4 <- data1985p1$VG4 # WEEKLY HRS WORKED MAIN JOB

g5 <- data1985p1$VG5 # WEEKLY HRS WORKED TOTAL
g5[g5 == 88 | g5 == 99] <- NA

g7 <- data1985p1$VG7 # HOW OFTEN PAID

g8 <- data1985p1$VG8 # GROSS PAY IN MAIN JOB
g8[g8 == 9999] <- NA

wage <- g8 * (as.numeric(g7 == 1) * 1 + as.numeric(g7 == 2) / 2 + as.numeric(g7 == 1) / 4)
wage[wage == 0] <- NA

hour1 <- as.numeric(g3 == 1) * g4
hour1[is.na(hour1)] <- 0
hour2 <- as.numeric(g3 == 2) * g5
hour2[is.na(hour2)] <- 0
hour <- hour1 + hour2
hour[hour == 0] <- NA

wage_hour <- wage / hour
wage_hour[wage_hour > 100] <- NA


DATA1 <- data.frame(t = 1, # time
                    wage_hour = wage_hour, # outcome
                    year_edu = year_edu, # treatment
                    attitude = attitude, # IV
                    sex = sex,
                    born_australia = born_australia,
                    capital = capital,
                    married = married,
                    uni_mem = uni_mem,
                    gov_emp = gov_emp,
                    age = age,
                    year_expe = year_expe)

DATA1$born_australia <- as.factor(DATA1$born_australia)
DATA1$married <- as.factor(DATA1$married)
DATA1$uni_mem <- as.factor(DATA1$uni_mem)
DATA1$gov_emp <- as.factor(DATA1$gov_emp)
DATA1 <- DATA1 %>% tidyr::drop_na(wage_hour)


################
## Final dataset
DATA <- rbind(DATA0, DATA1) %>% mice() %>% complete() %>% sample() %>% mutate(
  born_australia = as.numeric(born_australia == 1),
  married = as.numeric(married == 1),
  uni_mem = as.numeric(uni_mem == 1),
  gov_emp = as.numeric(gov_emp == 1)
)
save(DATA, file = "ALS_DATA.RData")



# Remark：协变量向量 X 包含以下六个变量：
# born_australia 表示是否出生于澳大利亚（1=是，0=否）；
# married 表示是否已婚（1=是，0=否）；
# uni_mem 表示是否为工会会员（1=是，0=否）；
# gov_emp 表示是否为政府雇员（1=是，0=否）；
# age 表示受访者年龄；year_expe 表示以年为单位计算的工作经验。
# captital: 居住城市是否在省会城市

# Y: wage_hour: 小时工资
# A: year_edu: 教育水平
# Z: attitude, t: 态度，时间
