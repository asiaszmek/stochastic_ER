import sys

import pandas as pd


import statsmodels.api as sm
import statsmodels.formula.api as smf

alldata = pd.read_csv(sys.argv[1])
for diam in [1.2, 2.4, 6.0]:
    new_data = alldata[alldata["model"]==diam]
    results2 = smf.ols('y ~ x*condition', data=new_data).fit()
    print(sm.stats.anova_lm(results2, typ=2))
    print(results2.summary())

