---
output:
  html_document: default
  word_document: default
  pdf_document: default
---
# How To Test & Roll 

Elea McDonnell Feit, Drexel University  
Ron Berman, The Wharton School

*Prepared for the July 2020 R Ladies Philly Meetup*

## Workshop Description
Marketers often use A/B testing as a tool to compare marketing treatments in a test stage and then deploy the better-performing treatment to the remainder of the consumer population. These tests have traditionally been analyzed using hypothesis testing. In this workshop, Elea and Ron will explain a new approach to A/B testing that they developed called "Test & Roll" that focuses on the profit earned during and after the test. It is based on a Bayesian decision theory and the (very technical) Test & Roll paper explains why marketers should use this approach. This workshop will focus on how to Test & Roll using hands-on examples in R. Elea and Ron will cover:

- When you should use "Test & Roll" versus traditional hypothesis test
- How to compute the Test & Roll sample size that maximizes profits
- How to estimate priors from your own data to help plan your Test & Roll
- How to analyze a Test & Roll experiment (It's easy!)

## Materials
If you want follow along, you can view the [slides]() in your browser. 

If you want to run the code as we go through the workshop, You will need to have R, R Studio and the `rstan` package installed. All the code is in the [R Markdown Code for the slides]() and you  will also need the file [`nn_functions.R`](). 

If you use git, you can also clone the entire workshop repository at (github.com/eleafeit/testandroll). The `howto` folder contains the R Markdown files for hte workshop. (There is other [good stuff in the repo](), including a copy of the paper and R scripts for replicating the analysis in the paper.) 

## Expectations
- We will assume you have a general idea of what a probability distribution (like the normal distribution) is. 
- We expect you are comfortable reading R. We will use R for mathematical calculations, for loops, and plotting. 
- You do not need to know how to plan and analyze an A/B test using hypothesis testing, but see Elea's tutorial on "Advanced A/B Testing() which was [recorded](https://www.youtube.com/watch?v=QXpYtM-Zlxg&t=4s) in a previous R Ladies Philly workshop, if you are interested.
