import pandas as pd
import csv
import time

source_path = r'C:\Users\galch\OneDrive\Документы\GitHub\recalculation_overall_statistics_with_high_res_from_refinement\project\full_automation\test_files\test_results.csv'
live_path = 'live.csv'

df = pd.read_csv(source_path, delimiter=';')
df = df.T.reset_index() 

with open(live_path, 'w', encoding='utf-8') as f:
    f.write(df.head(0).to_csv(sep=';', index=False)) 
for i, (idx, row) in enumerate(df.iterrows(), start=1):
    with open(live_path, 'a', encoding='utf-8') as f:
        row.to_frame().T.to_csv(f, header=False, sep=';', index=False)
    print(f"Written row {i}/{len(df)}")
    time.sleep(1)
