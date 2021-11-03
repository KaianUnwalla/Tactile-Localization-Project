

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
theme_set(kaian_theme(base_size = 20))

Pivot <- function(alldata1){ 
  #This takes the raw TOJ data and summarizes it for each participant
  pivot_full<-alldata1 %>% 
    filter(Actual>=0) %>% #remove timeouts
    group_by(Participant,Condition,Hands,SOA) %>%
    summarise(numTrials = length(Actual), numRF = sum(Actual), Data = mean(Actual)) %>% 
    ungroup()
  return(pivot_full)
}

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

compLikelihood <- function(Hypothesis,numTrials,numRF){
  #this uses the binomial equation to calculate the likelihood of the data given the hypothesis
  p = (numRF*log(Hypothesis))+((numTrials-numRF)*log(1-Hypothesis))
  return(p)
}

genHypothesisMLE <- function(wInt,wExt,Hands,SOA,numTrialsCondition, Condition=NULL){
  #this takes the chosen internal and external weight and creates the psychometric curve for the crossed and uncrossed posture at each SOA.
  #uncrossed: sum of the weights, crossed: difference between weights
  
  #The data were sampled across SOA while the hands and task paraemeter were varied.
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
