import pandas as pd

summary = pd.read_csv("data/interactions/summaries.csv") 

summary["diversity bacteria"] = 2**summary["entropy_bacteria"]
summary["diversity phages"] = 2**summary["entropy_phages"]
summary["V / B"] = summary["number_phages"] / summary["number_bacteria"]

pivot = pd.pivot_table(summary, values=["number_bacteria", "V / B", "diversity bacteria",
                                        "diversity phages"],
                                index=["infecttype", "nspecies"])
pivot = pivot[["number_bacteria", "V / B", "diversity bacteria", "diversity phages"]]
print(pivot.to_latex())

# table for interactions

# only use 10 species
nsp = 5
#summary = summary.loc[summary.nspecies==nsp,:] 

summary.loc[:,'ones'] = 1 

tab_bact = summary.pivot_table(index=["infecttype", "nspecies"], values="ones", columns=["final_bact_sp"], aggfunc="count")
tab_phage = summary.pivot_table(index=["infecttype", "nspecies"], values="ones", columns=["final_phage_sp"], aggfunc="count")

tab_both = pd.concat([tab_bact, tab_phage], axis=1)
print(tab_both.to_latex())