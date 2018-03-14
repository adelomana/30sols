import seaborn as sns

import matplotlib,matplotlib.pyplot

tips = sns.load_dataset("tips")
axo = sns.violinplot(x="day", y="total_bill", data=tips,inner=None, color=".8")
#axo = sns.stripplot(x="day", y="total_bill", data=tips, jitter=True)
axo = sns.swarmplot(x="day", y="total_bill", data=tips)
fig = axo.get_figure()
fig.savefig("output.pdf")

