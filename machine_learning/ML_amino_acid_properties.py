import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text
import numpy as np

helix_propensity = {
    'A': 1.45, 'C': 0.77, 'D': 0.98, 'E': 1.53,
    'F': 1.12, 'G': 0.53, 'H': 1.24, 'I': 1.00,
    'K': 1.07, 'L': 1.34, 'M': 1.20, 'N': 0.73,
    'P': 0.59, 'Q': 1.17, 'R': 0.79, 'S': 0.79,
    'T': 0.82, 'V': 1.14, 'W': 1.14, 'Y': 0.61
}
acetylation = {
    "A": 1, "C": 1, "T": 1, "S": 1, "V": 1, 
    "G": 1, "D": 1, "E": 1, "N": 1, "Q": 1, 
    "L": 1, "I": 1, "F": 1, "Y": 1, "K": 1,
    "P": 0, "H": 0, "M": 0, "R": 0, "W": 0}

import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt

# 1. Logistische Regressionsgewichte für Position 2 
'''logreg_weights = {
    'A': 0.257219780960474, 'C': -0.0413896992435331, 'D': 0.00364784309706693, 'E': -0.028789513900349,
    'F': 0.0764302840296069, 'G': 0.0423468601965782, 'H': -0.0363504607015354, 'I': 0.0680002439953118,
    'K': -0.092439157225919, 'L': 0.190410346507222, 'M': 0.0279109389964102, 'N': -0.0818165508890066,
    'P': 0.0298756635173394, 'Q': 0.0226925298464106, 'R': -0.00856512370737523, 'S': 0.0479925854507557,
    'T': 0.0312569063095202, 'W': 0.101211871882227, 'Y': 0.00137822826722116
}'''
logreg_weights = {
    'G': 0.13877080199116,
    'E': 0.0749801254486483,
    'D': 0.0710235618197539,
    'R': 0.0642563387617163,
    'T': 0.0491423221436578,
    'Y': 0.0365315350644674,
    'Q': 0.0365313728589567,
    'C': 0.0342373890977853,
    'K': 0.0265010807097719,
    'N': 0.0259348467169903,
    'A': 0.0258448816760883,
    'F': 0.0252048195745963,
    'L': 0.0248885846054995,
    'S': 0.0160337033857681,
    'W': 0.0096765601312246,
    'P': 0.0061863312026707,
    'I': 0.0052557176092602,
    'M': 0.0020424450498045,
    'H': 0.0000562618023934408
}


# V nicht mit einbezogen
logreg_df = pd.DataFrame(list(logreg_weights.items()), columns=["AA", "logreg_weight"])

# 2. Aminosäure-Eigenschaften hinzufügen
aa_properties = {
    'A': {'hydrophobicity': 1.8, 'pI': 6.0, 'volume': 88.6},
    'C': {'hydrophobicity': 2.5, 'pI': 5.1, 'volume': 108.5}, 
    'D': {'hydrophobicity': -3.5, 'pI': 2.8, 'volume': 111.1},
    'E': {'hydrophobicity': -3.5, 'pI': 3.2, 'volume': 138.4},
    'F': {'hydrophobicity': 2.8, 'pI': 5.5, 'volume': 189.9},
    'G': {'hydrophobicity': -0.4, 'pI': 6.0, 'volume': 60.1},
    'H': {'hydrophobicity': -3.2, 'pI': 7.6, 'volume': 153.2},
    'I': {'hydrophobicity': 4.5, 'pI': 6.0, 'volume': 166.7},
    'K': {'hydrophobicity': -3.9, 'pI': 9.7, 'volume': 168.6},
    'L': {'hydrophobicity': 3.8, 'pI': 6.0, 'volume': 166.7},
    'M': {'hydrophobicity': 1.9, 'pI': 5.7, 'volume': 162.9},
    'N': {'hydrophobicity': -3.5, 'pI': 5.4, 'volume': 114.1},
    'P': {'hydrophobicity': -1.6, 'pI': 6.3, 'volume': 112.7},
    'Q': {'hydrophobicity': -3.5, 'pI': 5.7, 'volume': 143.8},
    'R': {'hydrophobicity': -4.5, 'pI': 10.8, 'volume': 173.4},
    'S': {'hydrophobicity': -0.8, 'pI': 5.7, 'volume': 89.0},
    'T': {'hydrophobicity': -0.7, 'pI': 5.6, 'volume': 116.1},
    ''''V': {'hydrophobicity': 4.2, 'pI': 6.0, 'volume': 140.0},'''
    'W': {'hydrophobicity': -0.9, 'pI': 5.9, 'volume': 227.8},
    'Y': {'hydrophobicity': -1.3, 'pI': 5.7, 'volume': 193.6}
}
natA = {'A', 'S', 'T', 'V', 'C', 'G'}
natC = {'L', 'I', 'F', 'W'}
natB = {'D', 'E', 'N', 'Q'}

# 3. Erstelle DataFrame mit Eigenschaften
prop_df = pd.DataFrame(aa_properties).T
prop_df["AA"] = prop_df.index
prop_df["helix_propensity"] = prop_df["AA"].map(helix_propensity)
prop_df["acetylation"] = prop_df["AA"].map(acetylation)
# integriere NatA, NatB, NatC
prop_df["NatA"] = prop_df["AA"].apply(lambda x: 1 if x in natA else 0)
prop_df["NatB"] = prop_df["AA"].apply(lambda x: 1 if x in natB else 0)
prop_df["NatC"] = prop_df["AA"].apply(lambda x: 1 if x in natC else 0)

# 4. Merge: Logreg-Daten + Eigenschaften
df = pd.merge(logreg_df, prop_df, on="AA")

# 5. Lineare Regression
X = df[["hydrophobicity", "pI", "helix_propensity", "acetylation"]]
X = sm.add_constant(X)
y = df["logreg_weight"]
model = sm.OLS(y, X).fit()
# 6. Ergebnisse anzeigen
print(model.summary())
# plot -log10(p-Werte) gegen coeffizienten
plt.figure(figsize=(10, 6))

# Farben basierend auf Signifikanz
colors = ['red' if p < 0.05 else 'gray' for p in model.pvalues[1:]]

plt.scatter(model.params[1:], -np.log10(model.pvalues[1:]), color=colors, s=100, alpha=0.7)
# Punkte beschriften mit adjust_text
texts = []
for i, txt in enumerate(model.params.index[1:]):
    texts.append(plt.text(model.params[i + 1], -np.log10(model.pvalues[i + 1]), txt, fontsize=10))

adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
plt.xlabel('coefficients')
plt.ylabel('-log10(p-value)')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()