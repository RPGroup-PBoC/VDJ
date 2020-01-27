---
layout: post
title: Cutting Probability Model Explorer
date: "2017-09-12 13:32:20 +0300"
description: pcut_explorer
permalink: pcut_model_explorer
---

---

In this work, we computed the probability of a given sequence resulting in
successful cleavage using a Bayesian definition of probability. Briefly, we modeled 
that the number of cleavage events given a number of observed looping events should 
be Binomially distributed with a cutting probability p. We assume that, *a priori*, 
the cleavage probability can be **any** value between 0 and 1. With this uniform prior 
probability in hand, we can enumerate the posterior probability of the cleavage probability 
as 

$$
P(p_{cut}\,\vert\, n_{loops}, n_{cuts}) = \frac{(n_{loops} + 1)!}{n_{cuts}!(n_{loops} - n_{cuts})!}p_{cut}^{n_{cuts}}(1 - p_{cut})^{n_{loops} - n_{cuts}}.
$$

The interactive figure given below allows you to explore how different number of 
observed cutting events and total observed looping events influences our estimation
of the probability of a given value of the cutting probability. This figure was created using
the [Bokeh plotting library](http://bokeh.pydata.org). The code used to
generate it can be downloaded from the [Code]({{site.baseurl}}/code) page of
this website. An offline version of the figure [can be downloaded
here]({{site.baseurl}}/figures/pcut_model_explorer.html).

<center>

{% include pcut_model_explorer.html %}

</center>
