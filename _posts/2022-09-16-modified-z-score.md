---
title: Modified Z-scores
toc: true
toc_sticky: false
tags: ["outlier-detection"]
comments: true
date: 2022-09-17
header:
  teaser: "../images/bocpd_graphical_model.png"

---

Z-score is the most widely used method to measure outlier-strength. However,
there exists modified versions that are more robust to the presence of other outliers.
IBM suggests the following estimator in their <a href="https://www.ibm.com/docs/en/cognos-analytics/11.1.0?topic=terms-modified-z-score">Cognos Analytics page</a>:
$$ z(x) = \begin{cases} 
\frac{x - \text{median}(x_{1},..., x_{n})}{1.486 \times \text{MedianAD}(x_{1},..., x_{n})}  \text{     if     } \text{  MedianAD}(x_{1},..., x_{n}) \neq 0 \\
\frac{x - \text{median}(x_{1},..., x_{n})}{1.253314 \times \text{MeanAD}(x_{1},..., x_{n})} \text{     else}  \end{cases} $$

While the first part has been well covered, e.g. in {% cite leys2013detecting%}. I liked the fact that IBM documentation in addition accounts for 
the case where the MedianAD (median absolute deviation) is equal to 0.

In this post, I am going to explain the rationale behind this modified version of the z-score. <br />
<br />First, I will outline 3 limitations of the non modified z-score.
<br />We will then see how we can solve for two of them. 
<br />By then we will be able to understand where those number 1.486 and 1.253314 come from (we will actually see that IBM made a typo,
it should be 1.4826 instead of 1.486 🙂).



# Limitations of Z-score
## Limitation 1: Normal distribution assumption
The idea behind Z-score is that if X follows a normal distribution $$\mathcal{N}(\mu, \sigma^{2})$$.
Then $$\frac{X -\mu}{\sigma}\sim \mathcal{N}(0, 1)$$.<br/>
An outlier can be loosely defined as a point that deviates from expected behaviour.<br/>
The probability of drawing an instance with a z-score in absolute value as extreme or higher than 3 would be of 0.03%.
The probability of drawing an instance with a z-score in absolute value as extreme or higher than 4 would be of .


In practice $$\mu$$ and $$\sigma$$ are unknown. <br/>
In finite sample, natural estimators of $$\mu$$ and $$\sigma$$ are respectively.
Let $$\{x_{1},..., x_{n}\}$$ denote our sample, Z-score is then defined as follow:
$$Z(x) = \frac{x - \bar{x}}{s}$$[^1] <br/>
Where $$\bar{x} = \frac{\sum_{i=1}^{n} x_{i}}{n}$$ and $$s = \sqrt{\frac{\sum_{i=1}^{n} (x_{i} -\bar{x})^{2}}{n -1}}$$ 

## Limitation 2: Z-score can perform very badly in small samples
<u> Property:</u>
Let $$\{x_{1},..., x_{n}\}$$ denote our sample. 
$$\forall x \in \{x_{1},..., x_{n}\} :  |Z(x)| \leq \frac{n-1}{\sqrt{n}}$$.

The proof of this property can be found in {% cite shiffler1988maximum%}. It is shown that $$\frac{n-1}{\sqrt{n}}$$ is 
the maximum and not only an upper-bound.
To convince yourself, if we take the following example from {% cite iglewicz1993detect%}: <br/>
$$ \forall i \in \left[\!\!\left[1 , n-1 \right]\!\!\right]x_{i} = x $$, and $$x_{n} = x + ny$$, we get
$$Z(x_{n}) = \frac{n-1}{\sqrt{n}}$$

This property implies that 

## Limitation 3: The mean and the standard deviation are strongly impacted by outliers

The example that we have used above illustrates a third point. $$\bar{x}$$ and $$s$$ are not robust statistics.
As <a href="https://en.wikipedia.org/wiki/Robust_statistics">wikipedia page</a>

# Here comes the modified Z-score to the rescue
## Solving for limitation 2 and 3

## The issue with issue the Median Absolute Deviation



[^1]: Please note that the z-score actually follows a t-distribution. However, in large samples the t-distribution very close to a normal distribution.