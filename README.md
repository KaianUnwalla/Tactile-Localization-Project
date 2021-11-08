## Project Description

This is code I made for a project during my PhD. I have included some hypothetical data as a sample to showcase the code. 

The data were collected to better understand how we locate touch to our hands. We locate a touch by integrating internal (touch) information and external (largely spatial) information. We measure this integration process using a tactile temporal order judgement task. Here, participants are asked to indicate which hand received a vibration first. This is completed with their hands uncrossed and with the hands crossed. Consistently, participants are less accurate and slower when their hands are in a crossed posture. This crossed-hands deficit is caused by a conflict between the internal and external information when the hands are crossed (touch tells you left hand, whereas spatially the hand is on the right side of space).

This effect is moderated by many other factors. Some such factors are visual information (performance improves with the eyes closed), and body posture (performance improves when lying down).

I have included code showing various analyses and visualizations that I  used. 

- Hypothetical data.Rmd allows you to create your own hypothetical data. You can then see how changes to certain parameters changes performance on the task.

- Basic analyses.Rmd provides the code for the typical visualizations used in most crossed-hands papers. It also conducts the typical frequentist statistics employed in this task.

- Participant specific model.Rmd employs a Bayesian maximum likelihood estimation in order to estimate the slope of the psychometric curves. The value of the slope is broken down into an internal weight and an external weight where the sum represents the slope of the uncrossed curve, and the difference is the slope of the crossed curve (see Unwalla, Goldreich & Shore, 2021 or Unwalla, Cadieux and Shore, 2021 for a full explanation of these weights).
    - Unwalla, K., Goldreich, D., & Shore, D. I. (2021). Exploring reference frame integration using response demands in a tactile temporal-order judgement task. Multisensory Research, 34(8), 807-838.
    - Unwalla, K., Cadieux, M. L., & Shore, D. I. (2021). Haptic awareness changes when lying down. Scientific Reports, 11(1), 1-7.

- Hierarchical model.Rmd employs a Markov Chain Monte Carlo simulation using the Bayesian maximum likelihood estimation to determine individual participant weights, but also the internal and external weights of the population. 

- Model functions.R is a file containing all the functions used for this project. 