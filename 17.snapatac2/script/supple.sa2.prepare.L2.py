import pandas as pd

barcode2id_file = "../resource/barcode2id.csv"
barcode2id = pd.read_csv(barcode2id_file, header = 0)
cluster2size = barcode2id.L1.value_counts()
with open("../resource/sa2_L1_cluster2size.csv", 'w') as f:
    f.writelines([f"{c},{s}\n" for c, s in zip(
        cluster2size.index.to_list(), cluster2size)])


