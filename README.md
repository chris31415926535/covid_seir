
<!-- README.md is generated from README.Rmd. Please edit that file -->

# covid\_seir

<!-- badges: start -->

<!-- badges: end -->

This simple Shiny app is an interactive simulation of how physical
distancing affects the spread of COVID-19 at the population level. It
puts you in the role of choosing how much physical distancing the whole
society will do: choose the level of interaction compared to “normal”
(100 is business as usual, 0 is no one interacts with anyone), then
press the button to see what happens over the next 30 days.

It’s based on the epidemiological model developed by Anderson SC et
al. (2020) and presented in [this
paper](https://www.medrxiv.org/content/10.1101/2020.04.17.20070086v1.full.pdf),
which is used by the Public Health Agency of Canada in its [official
reports](https://www.canada.ca/content/dam/phac-aspc/documents/services/diseases-maladies/coronavirus-disease-covid-19/epidemiological-economic-research-data/update-covid-19-canada-epidemiology-modelling-20210115-en.pdf)
(p. 9). The authors have created **covidseir**, [a wonderful R package
that you can download and install from
GitHub\!](https://github.com/seananderson/covidseir)

This simulation is very simple: you choose the degree of physical
distancing as a percentage of interactions (100% is business as usual,
0% is no interactions whatsoever) and click “Go,” and the simulation
shows you what happens over the next 30 days.

# How does it work?

The model is a variant of [the epidemiological SEIR
model](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SEIR_model)
of how diseases spread. Please see [the
paper](https://www.medrxiv.org/content/10.1101/2020.04.17.20070086v1.full.pdf)
for details, but at a high level it tracks how people move through six
states over time:

  - \(S\) : Susceptible to the virus;
  - \(E_1\) : Exposed to the virus;
  - \(E_2\) : Exposed, pre-symptomatic, and able to infect others;
  - \(I\) : Symptomatic and able to infect others;
  - \(Q\) : quarantined; and
  - \(R\) Recovered or deceased.

There are also six variables for people practicing physical distancing:
\(S_d\), \(E_1_d\), \(E_2_d\), \(I_d\), \(Q_d\), and \(R_d\). People are
assumed to be immune once they enter state \(R\).

The model describes how these groups change over time with differential
equations. I won’t repeat them here (you’re welcome), but I *will* [link
to the paper one more time so you can see them yourselves on pages 6
and 7\!](https://www.medrxiv.org/content/10.1101/2020.04.17.20070086v1.full.pdf)

In brief, people start out \(S\)usceptible, then get \(E_1\)xposed, then
after a while they are both \(E_2\)xposed and able to infect others,
then they become visibly \(I\)nfected and symptomatic, then either go
straight to \(R\)ecovery, or into \(Q\)uarantine and then into
\(R\)ecovery.

At the same time, people can either start or stop physical distancing
when they’re in any of these states; the model assumes that there is
some flow back and forth from distancing and not-distancing, but that in
general more people start distancing than stop.

# Please note\!

  - Policy changes in physical distancing happen right away: we assume
    society as a whole turns on a dime, responsive to your whim.
  - The simulation uses Anderson et al.’s (2020) initial conditions,
    which represent British Columbia.
  - The simulation solves the differential equations numerically (using
    the **deSolve** package) using a resolution of 0.1 days, and then
    keeps results for whole-numbered time-steps for plotting.
  - **I’m not an epidemiologist** and am not affiliated with Anderson et
    al. in any way. I just like differential equations.
  - **This is for edutainment purposes only** to learn about how
    physical distancing can affect COVID-19 prevalence on the population
    level.
      - It comes with no warranties, guarantees, or manatees\!
