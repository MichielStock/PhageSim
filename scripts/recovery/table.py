import pandas as pd

summary = pd.read_csv("data/recovery/summaries.csv") 

summary["ones"] = 1

tab_bact = summary.pivot_table(index=["infecttype", "nspecies"], values="ones", columns=["final_bact_sp"], aggfunc="count")

# pick only those with bacteria
summary = summary.loc[summary.number_bacteria>0,:]

summary["diversity bacteria"] = 2**summary["entropy_bacteria"]
summary["diversity phages"] = 2**summary["entropy_phages"]
summary["V / B"] = summary["number_phages"] / summary["number_bacteria"]

pivot = pd.pivot_table(summary, values=["number_bacteria", "V / B", "diversity bacteria",
                                        "diversity phages"],
                                index=["infecttype", "nspecies"])
pivot = pivot[["number_bacteria", "V / B", "diversity bacteria", "diversity phages"]]


table = pd.concat([pivot, tab_bact], axis=1)
print(table.to_latex())