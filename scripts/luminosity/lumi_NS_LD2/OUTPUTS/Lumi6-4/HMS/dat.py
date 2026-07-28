import pandas as pd 

df = pd.read_csv("yield_data_LD2.csv")
new_df = df [["current","yieldRel_HMS_track","uncern_yieldRel_HMS_track"]]
new_df.to_csv("yield_data_LD2_track.csv")
print(new_df)
