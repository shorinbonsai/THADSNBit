import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import matplotlib

#SMALL_SIZE = 10
#EXTRA_SMALL = 5
#matplotlib.rc('font', size=SMALL_SIZE)
#matplotlib.rc('axes', titlesize=SMALL_SIZE)

dataframe3 = pd.read_csv("./ED30alphaT.csv", error_bad_lines=False, skiprows=[0])
fig, ax = plt.subplots(figsize=(14, 8))
ax.set_ylim(10, 75)
plt.xlabel("Parameter settings 1 through 8, each group of three plots Left -> Right  1, 3, 5 Starting edges")
plt.ylabel("Fitness")
plt.title("Epidemic Duration - 30 Alpha")

ax = sns.boxplot(data=dataframe3, notch=True)
[ax.axvline(x+.5,color='k') for x in ax.get_xticks()]

plt.savefig('ED30alpha.png')
plt.show()
plt.close()

dataframe3 = pd.read_csv("./ED40alphaT.csv", error_bad_lines=False, skiprows=[0])
fig, ax = plt.subplots(figsize=(14, 8))
ax.set_ylim(10, 75)
plt.xlabel("Parameter settings 1 through 8, each group of three plots Left -> Right  1, 3, 5 Starting edges")
plt.ylabel("Fitness")
plt.title("Epidemic Duration - 40 Alpha")

ax = sns.boxplot(data=dataframe3, notch=True)
[ax.axvline(x+.5,color='k') for x in ax.get_xticks()]

plt.savefig('ED40alpha.png')
plt.show()
plt.close()
