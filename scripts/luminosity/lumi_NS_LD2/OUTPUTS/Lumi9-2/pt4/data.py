import pandas as pd 

df = pd.read_csv("yield_data_LD2.csv")
new_df = df [["current","yieldRel_SHMS_scaler","uncern_yieldRel_SHMS_scaler"]]
new_df.to_csv("yield_data_LD2_SHMS_Sc.csv")
print(new_df)
