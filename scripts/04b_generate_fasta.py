import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()

df = pd.read_csv(args.input)
if 'pre_filter_pass' in df.columns:
    df = df[df['pre_filter_pass']]

with open(args.output, 'w') as f:
    for _, row in df.iterrows():
        f.write(f">{row['grna_id']}\n{row['spacer']}\n")
