import pandas as pd

df1 = pd.read_csv("lumi_data.csv")
df2 = pd.read_csv("lumi_data_LD2.csv")
#df3 = pd.read_csv("16718lumi_data.csv")
#df4 = pd.read_csv("16719lumi_data.csv")
#df5 = pd.read_csv("16721lumi_data.csv")
#df6 = pd.read_csv("16722lumi_data.csv")
#df7 = pd.read_csv("16723lumi_data.csv")
#df8 = pd.read_csv("16725lumi_data.csv")
#df9 = pd.read_csv("16726lumi_data.csv")

combined = pd.concat([df1, df2], ignore_index=True)

combined.to_csv("lumi_data_LD2.csv", index=False)
