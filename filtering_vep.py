import pandas as pd

# Load annotated variant file
df = pd.read_csv("/content/TNBC01.VEP.output.txt", sep="\t")
print(df.shape)
print(df.columns.tolist())
