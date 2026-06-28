# scripts/02b_inspect_candidates.py
import pandas as pd

df = pd.read_csv("data/processed/candidate_grnas_raw.csv")
print(df.shape)
print(df["pre_filter_pass"].value_counts())
print(df[df["pre_filter_pass"]]["gc_content"].describe())

# Export filtered-only set
df_pass = df[df["pre_filter_pass"]].copy()
df_pass.to_csv("data/processed/candidate_grnas_filtered.csv", index=False)
print(f"\n[OK] {len(df_pass)} candidates passed pre-filter and saved.")
