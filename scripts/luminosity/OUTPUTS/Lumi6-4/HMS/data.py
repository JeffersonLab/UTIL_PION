import pandas as pd 

df = pd.read_csv("yield_data_LD2.csv")
new_df = df [["current","yieldRel_HMS_scaler","uncern_yieldRel_HMS_scaler"]]
new_df.to_csv("yield_data_LD2_Sc.csv")
print(new_df)
