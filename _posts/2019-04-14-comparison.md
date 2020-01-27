---
layout: post
title: Synthetic-Endogenous RSS Comparison 
description: Synthetic-Endogenous RSS Comparison 
permalink: comparison
---

---

We investigated how each point mutation in the reference RSS (V4-57-1) compared
to the various endogenous sequences. In the figure below, we show the behavior
of the chosen endogenous sequence along with the behavior of each point mutation
in the reference sequence needed to recover the endogenous sequence. To use this
tool, choose your desired endogenous sequence from the drop down menu. This will
bring up the the endogenous sequence where each point mutation is colored
differently. The plots below will show the behavior of each individual point
mutation. To look at a single mutation in isolation, hover your mouse over the
mutation in the endogenous sequence. The selected mutation will retain its color
while all other point mutants will be displayed in a faint grey color. The
looping frequency of the endogenous sequence (top-left panel) is shown in grey
at the far right of the plot. 

This figure was created using the [Bokeh plotting
library](http://bokeh.pydata.org). The code used to generate this figure can
be found on the [code page]({{site.baseurl}}/code/) of this website. An offline version of the figure [can be downloaded here]({{site.baseurl}}/figures/point_endogenous_comparison.html).

<center>

{% include point_endogenous_comparison.html %}

</center>
