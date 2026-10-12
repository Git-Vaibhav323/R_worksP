# Final project summary — statistical concepts demonstrated

**Title.** Statistical Analysis of Rainfall Variability and Regional Monsoon Patterns in India  
**Data.** IMD all-India area-weighted rainfall (1901–2015) and 36 IMD subdivisions (1901–2017)  
**Implementation.** Base R 4.6.1; scripts `R_code/01`–`11`

## Concepts covered (with where they appear)

| Concept | Where |
| --- | --- |
| Data quality, missing vs zero, duplicates, consistency of seasonal totals | `01`–`02`; report §8 |
| Mean, median, mode, range, variance, SD, quartiles, percentiles, IQR, CV, skewness, kurtosis, SE | `03`; report §10 |
| Empirical and conditional probability; IMD-style deficient/normal/excess | `03`; report §11 |
| Histogram, density, ECDF, Q–Q plots | `04`; Figures 4–8 |
| Normal, gamma, exponential; MLE; AIC; Shapiro–Wilk | `04`; Table 5 |
| Parametric-bootstrap KS (not naive ks.test p-values) | `04` |
| Zero rainfall hurdle (Bernoulli + gamma on positives) | `04`; West Rajasthan December |
| Sampling distribution of the mean; CLT simulation (fixed seed) | `05`; Figure 9 |
| *t* confidence intervals 90/95/99%; Wald interval for a proportion | `05`; Table 7 |
| One-sample *t* and Wilcoxon signed-rank | `06` |
| Welch two-sample *t* and Wilcoxon rank-sum | `06` |
| Paired *t* and paired Wilcoxon | `06` |
| Chi-square GOF on categories; chi-square independence | `06` |
| Pearson and Spearman correlation; covariance; heatmap | `07`; Figure 11 |
| Simple and multiple linear regression; *R*²; coefficient tests | `07` |
| Residual plots, Cook’s distance, Durbin–Watson, heteroscedasticity screen | `07`; Figure 14 |
| One-way ANOVA, Welch ANOVA, Kruskal–Wallis, Fligner–Killeen | `08` |
| Tukey HSD and Bonferroni pairwise Welch/Wilcoxon | `08` |
| Regional means, SD, CV | `08`; Figures 15–16 |
| Time-series plot, moving average, ACF/PACF | `09`; Figures 18–20 |
| OLS trend vs Mann–Kendall vs Sen’s slope | `09`; Table 11 |

## Analyses deliberately not done

Chi-square on raw millimetres; independent *t* on two seasons of the same year; correlating annual totals with their monsoon component; deleting extremes; imputing missing months; climate-change attribution.

## Headline numerical results (from saved tables)

- All-India annual mean 1182.03 mm (n = 115); 95% *t* CI 1161.59–1202.48 mm  
- June–September mean 890.26 mm (about 75.3% of the annual mean)  
- Deficient / normal / excess monsoon years: 25 / 74 / 16  
- Shapiro–Wilk on annual totals: *p* = 0.8014; exponential bootstrap KS *p* = 0  
- One-sample *t* vs 1180 mm: *t* = 0.197, *p* = 0.8441; do not reject H0  
- Kerala vs West Rajasthan Welch *t* = 64.77, *p* < 0.0001  
- National annual Mann–Kendall *p* = 0.3895; Sen slope −0.27 mm/year  
- Annual ~ year: *R*² = 0.008, slope *p* = 0.3423  
