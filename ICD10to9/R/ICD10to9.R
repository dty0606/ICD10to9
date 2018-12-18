#' ICD Transformation and mortality calculation
#'
#' This function is able to identify ICD9 code(s) matched ICD10 code(s) of interest, and then
#'  calculate the disease-specific age-adjusted mortality rate and
#' years of potential life lost (YPLL).
#'
#' Three datasets should be prepared for this function. \cr
#'      1. a death dataset which contains gender, age, and ICD code of death reason;\cr
#'      2. a population dataset which contains numbers of people for each gender and age; \cr
#'      3. a dictionary which contains ICD10 codes and their corresponding ICD9 codes.\cr \cr
#' After inputting any ICD10 code, 'ICD10to9' will search for death records of both ICD10 and corresponding ICD9 code.
#' And then, the calculation will be based on these death records. One limitation of the current version is
#' that the datasets should be preprocessed according to the description in Arguments section. Our team is trying to make it more flexible. The package will be updated later.\cr \cr
#' Within this package, a sample dataset which is death records in Chile, Chile population dataset from UN, and
#' a dictionary for cardiovascular disease are provided.

#' @param sample   A dataset of death records.\cr
#' It should contain 4 columns: 1. SEXO, 2. EDAD_TIPO, 3. EDAD_CANT, 4. dfr_col.\cr
#' SEXO: 1=males, 2=females; EDAD_TIPO: Type of age measure, 1=years, 2= months, 3 =days, 4=hours and minutes; EDAD_CANT: Quantity of age (in EDAD_TIPO measure); dfr_col: ICD code for death reason of this individual.
#' @param dfr_col         The name of ICD9/10 code column in SAMPLE dataset.
#' @param Diag.ICD9       Diagnosis ICD10 code(s) for disease(s) of interest.\cr
#' It should be the same format as 'ICD10' in DICTIONARY dataset.
#' @param dictionary      A dataset with two columns which users prepared for ICD transformation.\cr
#' It should contain 2 columns: 1. ICD10, 2. ICD9.\cr
#' ICD10: ICD10 code, ICD9: corresponding codes in ICD9.\cr
#' IMPORTANCE: this package can only identify ICD9 codes included in this DICTIONARY dataset. Make sure that ICD9 code in the DICTIONARY has all possible formats of ICD9 in the SAMPLE dataset.\cr
#' @param population      A Dataset of population.\cr
#' It should contain 3 columns: 1. Sex, 2. Age, 3. Value.\cr
#' Sex: 'Both Sexes', 'Male', 'Female'; Age: numeric; Value: numeric, the number of individuals of a certain age and gender
#' @param age.cut         User-defined age group cuts.
#' @param YPLL.yr         Expected years of life.
#' @return I. unadjusted.death.rate$total.death.rate: unadjusted death rate for total population by user-defined age group;\cr
#' II. unadjusted.death.rate$male.death.rate: unadjusted death rate for males by user-defined age group;\cr
#' III. unadjusted.death.rate$female.death.rate: unadjusted death rate for females by user-defined age group;\cr
#' IV. direct: age-adjusted death rates using direct method by gender\cr
#' V. indirect: standardized mortality ratio (SMR) by gender\cr
#' VI. YPLL: Years of potential life lost
#' @export
#' @examples 
#' data(Chile.death.1994.ICD9)
#' data(Chile.population.1998)
#' data(ICD10to9CVD)
#' ICD10to9(sample = Chile.death.1994.ICD9, dfr_col = 'DIAG1', Diag.ICD10 = c('I21.0'), 
#' dictionary = ICD10to9CVD, pop.data = Chile.population.1998, age.cut = c(1, 5, 15, 25, 45, 65, 80), YPLL.yr = 80)

# ICD10 Details Aware that the retranslation might not produce the same original data because different codes from ICD9
# might be matched with the same ICD10 code.

# sample <- read.table('C:/Users/chc217/Box Sync/R package project/DeathsChile2004ICD10.csv',header = TRUE, sep = ',')
# pop.data <- read.table('C:/Users/chc217/Box Sync/R package project/Chile2014Population.csv',header = TRUE, sep = ',')
# dictionary <- read.table('C:/Users/chc217/Box Sync/R package project/demo.csv',header = TRUE, sep = ',') icd10to9(sample,
# dfr_col = 'DIAG1', Diag.ICD10 = c('I210'), dictionary, pop.data,age.cut = c(1, 5, 15, 25, 45, 65, 80), YPLL.yr = 80)


ICD10to9 <- function(sample, dfr_col, Diag.ICD10, dictionary, pop.data, age.cut, YPLL.yr) {
    
    ## ask user to check input documents, if not, stop
    fun <- function() {
        answer = readline(cat("Are you using correct documents(y/n)?\n"))
        while (!answer %in% c("y", "n")) answer = readline(cat("Please enter y/n:\n"))
        return(answer)
    }
    check <- fun()
    if (check != "y") 
        stop("Please use correct documents, quit now.\n")
    
    ### warning msg warning msg for YPLL.yr
    if (!is.numeric(YPLL.yr)) 
        warning("YPLL.yr is not numeric\n") else if (YPLL.yr < 0) 
        warning("YPLL.yr is less than 0\n")
    ## warning msg for dfr_col
    if (!dfr_col %in% colnames(sample)) 
        warning("dfr_col is not in sample column list\n")
    ## warning msg for age.cut
    if (age.cut[length(age.cut)] > 80) {
        warning("Age > 80, will all be coerced to age group: 80+\n")
        age.cut = c(age.cut[age.cut <= 80])
    }
    
    ## translate input ICD10 to ICD9
    if (missing(dictionary)) 
        stop("Please load dictionary")
    if (missing(Diag.ICD10)) {
        ICD9.dictionary <- dictionary
    } else {
        ICD9.dictionary <- dictionary[dictionary$ICD10 %in% c(Diag.ICD10), 2]
    }
    
    if (sum(ICD9.dictionary %in% as.character(unique(sample[, dfr_col]))) == 0) 
        stop("Sample doesn't include any death record with disease matching to input ICD10 code. \n")
    
    ## extract all death records matching to input ICD10 code for following mortality calculation
    sample = sample[sample[, dfr_col] %in% ICD9.dictionary, ]
    
    ## age group creation
    age = ifelse(sample$EDAD_TIPO == 1, sample$EDAD_CANT, 0)
    age.matrix = as.data.frame(cbind(age, 0))
    names(age.matrix) = c("Age", "Age.cat")
    ## assign age group to each age
    for (i in 1:(length(age.cut) + 1)) {
        if (i == 1) 
            age.matrix[age < age.cut[i], 2] = i else if (i <= length(age.cut)) 
            age.matrix[age.cut[i - 1] <= age & age < age.cut[i], 2] = i else age.matrix[age.cut[i - 1] <= age, 2] = i
    }
    ## calculate number of deaths in each age group
    sample$age.cat = age.matrix$Age.cat
    sample.both.group = as.numeric(tabulate(sample$age.cat))
    sample.male.group = as.numeric(tabulate(sample[sample$SEXO == 1, ]$age.cat))
    sample.female.group = as.numeric(tabulate(sample[sample$SEXO == 2, ]$age.cat))
    ## preprocess UN data: coerce character to numeric, change 80+ to 80, only keep rows for each age level
    pop.data = pop.data[pop.data$Area == "Total", ]
    pop.both <- pop.data[pop.data$Sex == "Both Sexes", ]
    pop.male <- pop.data[pop.data$Sex == "Male", ]
    pop.female <- pop.data[pop.data$Sex == "Female", ]
    pop.both$Age = suppressWarnings(as.numeric(as.character(pop.both$Age)))
    pop.both$Age[nrow(pop.both)] = 80
    pop.both = na.omit(pop.both[, -10])
    pop.male$Age = suppressWarnings(as.numeric(as.character(pop.male$Age)))
    pop.male$Age[nrow(pop.male)] = 80
    pop.male = na.omit(pop.male[, -10])
    pop.female$Age = suppressWarnings(as.numeric(as.character(pop.female$Age)))
    pop.female$Age[nrow(pop.female)] = 80
    pop.female = na.omit(pop.female[, -10])
    
    ## calculate age group population
    pop.sum <- function(data) {
        g.pop = vector()
        for (i in 1:(length(age.cut) + 1)) {
            if (i == 1) 
                g.pop[i] = sum(data$Value[data$Age < i]) else if (i <= length(age.cut)) 
                g.pop[i] = sum(data$Value[age.cut[i - 1] <= data$Age & data$Age < age.cut[i]]) else g.pop[i] = sum(data$Value[age.cut[i - 1] <= data$Age])
        }
        return(g.pop)
    }
    pop.both.group = pop.sum(pop.both)
    pop.male.group = pop.sum(pop.male)
    pop.female.group = pop.sum(pop.female)
    
    ## rate calculation and age group labels
    age.name <- function(age.cut) {
        name = vector()
        for (i in 1:(length(age.cut) + 1)) {
            if (i == 1) 
                name[i] = paste("Age:<", age.cut[i], sep = "") else if (i <= length(age.cut)) 
                name[i] = paste("Age:", age.cut[i - 1], "-", age.cut[i] - 1, sep = "") else name[i] = paste("Age:", age.cut[i - 1], "+", sep = "")
        }
        return(name)
    }
    age.group.level = c(age.name(age.cut), "Total")
    both.rate = c(sample.both.group/pop.both.group, sum(sample.both.group)/sum(pop.both.group)) * 1000
    names(both.rate) = age.group.level
    male.rate = c(sample.male.group/pop.male.group, sum(sample.male.group)/sum(pop.male.group)) * 1000
    names(male.rate) = age.group.level
    female.rate = c(sample.female.group/pop.female.group, sum(sample.female.group)/sum(pop.female.group)) * 1000
    names(female.rate) = age.group.level
    
    ## Q1 direct method
    exp.male.count = pop.both.group * male.rate[-length(male.rate)]/1000
    exp.female.count = pop.both.group * female.rate[-length(female.rate)]/1000
    age.adj.male.rate = sum(exp.male.count)/sum(pop.both.group) * 1000
    age.adj.female.rate = sum(exp.female.count)/sum(pop.both.group) * 1000
    age.adj.mf.rate = age.adj.male.rate/age.adj.female.rate
    direct = c(age.adj.male.rate, age.adj.female.rate, age.adj.mf.rate)
    names(direct) = c("Age-adjusted death rate for male (per 1000)", "Age-adjusted death rate for female(per 1000)", "Age-adjusted rate ratio M/F=")
    
    ## Q2 indirect method
    exp.male.count.ind = pop.male.group * both.rate[-length(both.rate)]/1000
    smr.male = sum(sample.male.group)/sum(exp.male.count.ind)
    exp.female.count.ind = pop.female.group * both.rate[-length(both.rate)]/1000
    smr.female = sum(sample.female.group)/sum(exp.female.count.ind)
    indirect = c(smr.male, smr.female)
    names(indirect) = c("Standardized mortality ratio (SMR) for Male", "Standardized mortality ratio (SMR) for Female")
    
    ## Q3 Years of potential life lost (YPLL)
    mid.point = vector()
    new.age.cut = age.cut[age.cut <= YPLL.yr]
    for (i in 1:length(new.age.cut)) {
        if (i == 1) 
            mid.point[i] = new.age.cut[i]/2 else mid.point[i] = (new.age.cut[i - 1] + new.age.cut[i])/2
    }
    yr.till = YPLL.yr - mid.point
    obs.dcount = sample.both.group[1:length(yr.till)]  ##no midpoint for last age level,warning here
    YPLL = sum(yr.till * obs.dcount)
    names(YPLL) = paste("Years of potential life lost (YPLL)")
    ## result
    result = list(unadjusted.death.rate = list(total.death.rate = both.rate, male.death.rate = male.rate, female.death.rate = female.rate), 
        direct = direct, indirect = indirect, YPLL = YPLL)
    cat("End\n")
    return(result)
}
