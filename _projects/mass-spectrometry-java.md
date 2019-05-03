---
title: Some Java utilities for chemical ionization mass spectrometry
layout: project
name: Some Java utilities for chemical ionization mass spectrometry
---

GUI interface can be come handy when perform simple calculator operations. Here I wrote two Java uitilities in Swing to perform ion cluster composition calculation and mass concentration calculation

1. [Combinatoric ion cluster composition calculator](#utility-1)
2. [Concentration calculator for ion-drift chemical ionizaiton mass spectrometry](#utility-2)
<!--more-->

---
## Combinatoric ion cluster composition calculator
{:#utility-1}
Reasoning mass spectrum for chemical ionzation mass spectrometry for unit mass resolution could be challenge. Even for high resolution mass spectrum, the number of combinatins of elemental formulas could explode when the ion clusters are large. I have developed a easy to use Java GUI to provide combinatoric calculation for possible combination of molecules and ions in a givien chemical system defined by the user. Basically, the user could provide possible molecules in their system. For a given unit mass of an ion, the combinations will be populated into the table, once the calculate button is clicked.
<!--more-->

![alt java-utility](/assets/img/java-utility.png)

---
## Concentration and E/N ratio calculator
{:#utility-2}
Ion drift chemical ionization mass spectrometry data is simple to analyze, but to make sure the data is of high quality, the E/N ratio (which discribed relative ratio between electrostatic kinetic energy and thermal kinetic energy). The larger E/N ratio is, the higher probability the ion cluster will get fragmented during movement with colliding with buffer gas.

![alt concentration-calculator](/assets/img/concentration-calculator.png)