import pandas as pd
df = pd.read_csv('participants.tsv', sep='\t')

# 1. Diagnosis counts
print(df['diagnosis'].value_counts(), '\n')

# 2. Prefix rule check
print(pd.crosstab(df['participant_id'].str[4:6], df['diagnosis']), '\n')

# 3. What the modality flags actually contain
print("dwi column raw values:", df['dwi'].unique())
print("T1w column raw values:", df['T1w'].unique(), '\n')

# 4. The number that matters: who has both DWI and T1
has_dwi = df['dwi'].notna()
has_t1  = df['T1w'].notna()
print(pd.crosstab(df['diagnosis'], has_dwi & has_t1), '\n')

# 5. Artifact flag distribution
print(pd.crosstab(df['diagnosis'], df['ghost_NoGhost']))