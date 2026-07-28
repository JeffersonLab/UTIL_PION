import pandas as pd 

df = pd.read_csv("yield_data_LD2.csv")
new_df = df [["rate_SHMS","yieldRel_SHMS_scaler","uncern_yieldRel_SHMS_scaler"]]
new_df.to_csv("yield_data_LD2_r.csv")
print(new_df)
