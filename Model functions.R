
a=as.numeric(Sys.time())
set.seed(a)
select <- dplyr::select

Pivot <- function(alldata1){ 
  #This takes the raw TOJ data and summarizes it for each participant
  pivot_full<-alldata1 %>% 
    filter(Actual>=0) %>% #remove timeouts
    group_by(Participant,Condition,Hands,SOA) %>%
    summarise(numTrials = length(Actual), numRF = sum(Actual), Data = mean(Actual)) %>% 
    ungroup()
  return(pivot_full)
}


############# Hierarchical Model Functions 
log_gaussian_prior <- function(x,mean,sigma){
  #This returns the log value of the probability density function for a normal distribution. 
  #This is equivalent to taking the log value of the dnorm function in R at larger values, but will avoid underflows when values get smaller
  answer <- -log(sigma) - 0.5*log(2*pi) - (0.5*((x-mean)/(sigma))^2)
  return(answer)
}

genHypothesis <- function(wInt,wExt,SOA,taskInt,taskExt, numTrialsCondition){
  #this takes the chosen internal and external weight and creates the psychometric curve for the crossed and uncrossed posture at each SOA.
  #uncrossed: sum of the weights, crossed: difference between weights
  
  #The data were sampled across SOA while the hands and task paraemeter were varied.
  #Condition: Allocentric and Somatotopic
  #wInt, wExt --> the weight applied to the two reference frames in the model.
  #SOA --> a list of values sampled in the experiment (coded in seconds)
  #numTrials --> when the probability of making a particular response is zero or one, 
  #              the smallest or largest value must be coverted to 1/2 the smallest unit
  #              of measure, which in this case is the trial--> so p(zero) = .5/numTrials.
  
  #empty list to be filled with the exponent for the psychometric function
  t <- c() 
  #A crossed first becase R is alphabetical
  t[1:length(wInt)] = wInt - wExt 
  #then A uncrossed
  t[(length(wInt)+1):(length(wInt)*2)] = wInt + wExt
  #then B crossed
  t[(length(wInt)*2+1):(length(wInt)*3)] = wInt*taskInt - wExt*taskExt
  #finally B uncrossed
  t[(length(wInt)*3+1):(length(wInt)*4)] = wInt*taskInt + wExt*taskExt 

  #takes the exponent and SOA and returns the hypothesized probability
  hyp<-function(SOA, t){
    x = list((1/(1+exp(SOA*t*-1))))
    return(x)
  }
  
  #creates a list of the hypothesized probability of the hypothesis
  p <-c()
  p <- sapply(SOA, hyp, t=t)
  p <- do.call(c, p)

  #removes 1 and 0 from probabilities as we are using log-likelihoods
  p = ifelse(p == 1, 1-(0.5/numTrialsCondition), ifelse(p == 0, 0.5/numTrialsCondition, p))
  return(p)
}

compLikelihood <- function(Hypothesis,numTrials,numRF){
  #this uses the binomial equation to calculate the likelihood of the data 
  #given the hypothesis
  #Kaian Unwalla
  #MSP Lab, 2020
  #dis May 24, 2020. 
  #Hypothesis <- a vector of proportions against which to compare the actual
  #              number of responses and trials. These are calculated in the genHypothesis function.
  #numTrials <- actual number of trials in the analysis. Some trials are removed
  #              because the participant timed out etc.
  #numRF <- the number of "right first" responses made by the participant
  #p = probability of the hypothesis, given the data. 
  
  p = (numRF*log(Hypothesis))+((numTrials-numRF)*log(1-Hypothesis))
  return(p)
}


PPMC <- function(PCD_easy, pRfirst){
  #480 observation (1 for each row of hypothesis), with 30 trials in it (number of trials in experiment)
  pRfirst <- c(rbinom(length(hypothesis_current),numTrialsCondition,hypothesis_current))/numTrialsCondition
  calcPPMC <- cbind(PCD_easy, pRfirst) %>%
    group_by(Participant,Condition) %>%
    spread(key = Hands, pRfirst) %>%
    mutate(PCD = ifelse(SOA < 0, (Crossed - Uncrossed),(Uncrossed - Crossed))) %>%
    summarize(PCD = sum(PCD)) %>%
    spread(key = Condition, PCD)
  
  A <- calcPPMC$A
  B <- calcPPMC$B
  
  aPPMC <- mean(A)
  bPPMC <- mean(B)
  allPPMC <- c(A = aPPMC, B = bPPMC)
  return(allPPMC)
}


############# Participant Specific Model Functions 

#This code calculates the maximum likelihood estimated weights
calc.MLE <- function(alldata1, weightMin, weightMax, spacing){
  #Summarizes the data, and nest so that all data fitted into one row
  participant_data<-Pivot(alldata1) %>% nest()
  
  #Selects range of internal and external weights for MLE
  wInt <- seq(from = weightMin, to = weightMax, by = spacing)
  wExt <- seq(from = weightMin, to = weightMax, by = spacing)
  
  #Calculates likelihood of the data
  Likelihoods <- expand.grid(wInt = wInt,wExt = wExt) %>% #creates tibble with all combinations of internal and external weights
    add_column(participant_data) %>% #adds nested participant data
    unnest_dt(data) %>% #this makes sure that all the data is available for each weight combinatin
    mutate(Hypothesis = genHypothesisMLE(wInt,wExt,Hands,SOA), 
           LogLikelihood = compLikelihood(Hypothesis,numTrials,numRF)) %>% 
    group_by(wInt,wExt,Participant,Condition) %>% 
    summarize(LogLikelihood = sum(LogLikelihood)) 
  
  
  #Converts the log likelihoods into likelihoods
  Normalized_likelihoods <- Likelihoods %>%
    group_by(Participant, Condition) %>%
    mutate(new_LogLikelihood = LogLikelihood - max(LogLikelihood), 
           Likelihood = exp(new_LogLikelihood)) %>%  #finds maximum log likelihood and subtracts it from all, then normalizes
    ungroup()
  
  return(Normalized_likelihoods)
}

genHypothesisMLE <- function(wInt,wExt,Hands,SOA,numTrialsCondition){
  #this takes the chosen internal and external weight and creates the psychometric curve for the crossed and uncrossed posture at each SOA.
  #uncrossed: sum of the weights, crossed: difference between weights
  
  #The data were sampled across SOA while the hands and task parameter were varied.
  #Condition: Allocentric and Somatotopic
  #wInt, wExt --> the weight applied to the two reference frames in the model.
  #SOA --> a list of values sampled in the experiment (coded in seconds)
  #numTrials --> when the probability of making a particular response is zero or one, 
  #              the smallest or largest value must be coverted to 1/2 the smallest unit
  #              of measure, which in this case is the trial--> so p(zero) = .5/numTrials.
  
  
  t = ifelse(Hands == 'Crossed', wInt - wExt, wInt + wExt)
  p = (1/(1+exp(SOA*t*-1)))
  p = ifelse(p == 1, 1-(0.5/numTrialsCondition), ifelse(p == 0, 0.5/numTrialsCondition, p))
  return(p)
}

max.MLE <- function(Normalized_likelihoods){
  #Finds most likely internal and external weight pair
  Max_likelihoods <- Normalized_likelihoods %>%
    group_by(Participant,Condition,wInt,wExt) %>%
    summarize(Likelihood = mean(Likelihood)) %>%
    group_by(Participant,Condition) %>%
    mutate(rank = rank(Likelihood, ties.method = "first")) %>% #sorts by most likely weight pair for each participant and condition. When there is a tie, it takes lowest combined values for internal and external weight
    filter(rank == max(rank)) %>% 
    select(Participant,Condition,wInt,wExt,Likelihood) %>% 
    ungroup()
  
  return(Max_likelihoods)
}


############## Data Visualization Functions

kaian_theme <- function(base_size = base_size){
  theme_minimal(base_size = base_size) %+replace% 
    theme(legend.title = element_blank(), 
          panel.border = element_rect(color = 'black', linetype = 'solid', fill=NA), 
          panel.background = element_rect(color='white'),
          panel.grid.major = element_line(size = 0.5),
          panel.grid.minor = element_line(size = 0.25),
          legend.position = 'bottom', 
          #legend.position = c(0, 1),
          #legend.justification = c('left','top'),
          legend.spacing = unit(0.5,'cm'),
          legend.key.width = unit(1,'cm'),
          plot.caption = element_text(hjust=0)
    )
}

theme_set(kaian_theme(base_size = 10))


Pivot_overall <- function(alldata1){
  #this takes the raw TOJ data and summarizes it based on the experiment condition and hand posture and SOA. It also calculates morey within-subject error bars
  pivot_full<-alldata1 %>% 
    group_by(Participant,Condition,Hands,SOA) %>%
    summarise(Data = mean(Actual)) %>%
    ungroup()  %>% 
    mutate(grand_mean = mean(Data), n = length(unique(Participant))) %>% 
    group_by(Participant) %>% 
    mutate(part_mean = mean(Data), norm_mean = Data - part_mean + grand_mean, Length = length(Participant)) %>% 
    group_by(Condition,Hands,SOA) %>% 
    summarise(Data = mean(Data), n = mean(n), morey = (var(norm_mean)*mean(Length))/(mean(Length)-1)) %>% 
    mutate(SE = sqrt(morey)/sqrt(n)) %>% 
    ungroup()
  return(pivot_full)
}

PCD <- function(alldata1){
  #this calculates the PCD score for each participant
  pivot_full <- Pivot(alldata1) %>% 
    group_by(Participant,Condition) %>% 
    select(Participant,Condition,Hands,SOA,Data) %>% 
    spread(key=Hands, Data) %>% 
    mutate(PCD = ifelse(SOA < 0, (Crossed - Uncrossed),(Uncrossed - Crossed))) %>% 
    summarize(PCD = sum(PCD)) %>% 
    ungroup()
  return(pivot_full)
}

PCD_error <- function(alldata1){
  #this calculates overall PCD score for each condition, and within-subject morey error bars
  pivot_full<-alldata1 %>% 
    PCD() %>% 
    mutate(grand_mean = mean(PCD), n = length(unique(Participant))) %>%
    group_by(Participant) %>%
    mutate(part_mean = mean(PCD), norm_mean = PCD - part_mean + grand_mean, Length = length(Participant)) %>%
    group_by(Condition) %>% 
    summarise(PCD = mean(PCD), n = mean(n), morey = (var(norm_mean)*mean(Length))/(mean(Length)-1)) %>% 
    mutate(SE = sqrt(morey)/sqrt(n)) %>% 
    ungroup()
  return(pivot_full)
}
