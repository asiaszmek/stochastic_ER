import sys

import pandas as pd
import pingouin as pg

alldata = pd.read_csv(sys.argv[1])
for diam in [1.2, 2.4, 6.0]:
    new_data = alldata[alldata["model"] == diam]
    print(pg.anova(dv="y", between=["condition", "x"], data=new_data, detailed=True))
    print(pg.pairwise_tests(data=new_data, dv="y",  between=["condition", "x"], padjust="holm"))
