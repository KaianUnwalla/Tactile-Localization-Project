# Project Description

## This is code I made for a project during my PhD. I have included some hypothetical data as a sample to showcase the code. 

## The data were collected to better understand how we locate touch to our hands. This was done using a tactile temporal order judgement task. Here, participants are asked to indicate which hand received a vibration first. This is completed with their hands uncrossed and with the hands crossed. Consistently, participants are less accurate and slower when their hands are in a crossed posture. This is called the crossed-hands deficit. 

## This effect is moderated by many other factors. Some such factors are visual information (performance improves with the eyes closed), and body posture (performance improves when lying down).

## I have included code showing various analyses and visualizations that I  used. 
### - Hypothetical data.Rmd allows you to create your own hypothetical data. You can then see how changes to certain parameters changes performance on the task.
### - Basic analyses.Rmd provides the code for the typical visualizations used in most crossed-hands papers. It also conducts the typical frequentist statistics employed in this task.
### - Participant specific model.Rmd employs a Bayesian maximum likelihood estimation in order to estimate the slope of the psychometric curves. The value of the slope is broken down into an internal weight and an external weight. s