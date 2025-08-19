import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import copy
import time
import os
import scipy.stats as stats
import glob


IN="/users/win-fmrib-analysis/gux928/biobank/tmp/PANDORA_testing/BHC/GWAS_2025_08_11"
OUT="/users/win-fmrib-analysis/gux928/biobank/tmp/PANDORA_testing/BHC/GWAS_2025_08_11"


Nphenos=len(glob.glob(IN+"/[0:9]???"))
for path_cnt in range(1, Nphenos+1):
    try:
        padded_number = f"{path_cnt:04}"
        start_time = time.time()
        inpath=IN+"/"+padded_number+"/"
        out_path=OUT+"/"
        in_file=inpath+"p1/0001.txt"
        hit_file = inpath+"hits/0001.txt"
        table_file = inpath+"Table1.txt"
        out_file = out_path+padded_number
        pheno = path_cnt
        bf_level = -np.log10( (10**-7.5) / Nphenos )
        gwas_level = 7.5  # Default GWAS significance level

        # Read significant SNPs
        hits = pd.read_csv(hit_file, sep=" ", header=0)
        hits.rename(columns={"chr": "chromosome", "rsid": "SNP", "pos": "position", "P.value": "p_value"}, inplace=True)
        hits["chromosome"] = hits["chromosome"].astype(str).replace({"X": "23", "0X": "23"}).astype(int)
        hits["chromosome"] = hits["chromosome"].apply(lambda x: f"{x:02d}")
        
        data_table = pd.read_csv(table_file, sep=" ", header=0)
        data_table.rename(columns={"chr": "chromosome", "rsid": "SNP", "p-value": "p_value"}, inplace=True)
        data_table["chromosome"] = data_table["chromosome"].astype(str).replace({"X": "23", "0X": "23"}).astype(int)
        data_table["chromosome"] = data_table["chromosome"].apply(lambda x: f"{x:02d}")

        old_hits=hits.copy()
        hits = hits.iloc[np.where(data_table['info'] > 0.3)[0], :].copy()
        # bad_snp_chr=data_table[data_table["info"] <= 0.3][["chromosome", "SNP"]]
        bad_snp_chr = data_table.loc[data_table['info'] <= 0.3, ['chromosome', 'SNP']]

        if hits.shape[0]>=1:
            if not os.path.exists(out_path+"p"):
                try:
                    os.makedirs(out_path)
                except:
                    print("boo! directory exists")
                print("there are at least 1 significant hit, let's do Manhattan plots")
                print(hits.shape)
                # Read GWAS data
                data = pd.read_csv(in_file, sep=" ", header=0)
                data.rename(columns={"chr": "chromosome", "MarkerName": "SNP", "pos": "position", "P.value": "p_value"}, inplace=True)

                # Convert chromosome X to numeric representation
                data["chromosome"] = data["chromosome"].astype(str).replace({"X": "23", "0X": "23"}).astype(int)
                data["chromosome"] = data["chromosome"].apply(lambda x: f"{x:02d}")  # Ensure 2-digit format
                data = data.sort_values(["SNP", "p_value"], ascending=[True, False]).drop_duplicates("SNP")
                
                # Filter data for matching chr/rsid pairs
                old_data=data.copy()
                data = data.merge(bad_snp_chr, on=['chromosome', 'SNP'], how='left', indicator=True).query('_merge == "left_only"').drop(columns=['_merge'])

                # Mark significant hits
                # data["is_hit"] = data.apply(lambda row: f"{row.chromosome} {row.SNP}" in hits.apply(lambda x: f"{x.chromosome} {x.SNP}", axis=1).values, axis=1)
                # Create a set of "chromosome SNP" strings from hits
                hit_set = set(hits["chromosome"].astype(str) + " " + hits["SNP"])

                # Use vectorized operations to check membership
                data["is_hit"] = (data["chromosome"].astype(str) + " " + data["SNP"]).isin(hit_set)

                # Convert p-values to -log10 scale
                data["log_pval"] = copy.deepcopy(data["p_value"])#-np.log10(data["p_value"])
                
                # Adjust positions for plotting
                cumulative_pos = 0
                chrom_offsets = {}
                for chrom in sorted(data["chromosome"].unique()):
                    chrom_data = data[data["chromosome"] == chrom]
                    min_pos = chrom_data["position"].min()
                    data.loc[data["chromosome"] == chrom, "adjusted_pos"] = chrom_data["position"] - min_pos + cumulative_pos
                    chrom_offsets[chrom] = cumulative_pos + (chrom_data["position"].max() - min_pos) / 2
                    cumulative_pos = data["adjusted_pos"].max() + 2e7  # Offset chromosomes

                # Define a color palette (you can tweak this or use seaborn's color palette)
                chromosomes = sorted(data["chromosome"].unique())  # Get unique chromosomes
                #c = matplotlib.colormaps.get_cmap("tab20")
                #vals = [(0.1212, 0.467, 0.706, 1.0),(1.0, 0.498, 0.055, 1.0)]
                vals = [(0.0, 0.0, 0.0, 1.0),(0.33, 0.33, 0.33, 1.0)]
                c = matplotlib.colors.ListedColormap(vals)
                colors = [c(i) for i in [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]]
                plt.figure(figsize=(10, 4))

                # Loop through chromosomes and plot separately to maintain color differentiation
                for i, chrom in enumerate(chromosomes):
                    chrom_data = data[data["chromosome"] == chrom]
                    plt.plot(chrom_data["adjusted_pos"], chrom_data["log_pval"], '.', markersize=2, color=colors[i])  # Use different colors

                # Add SNP labels for significant hits
                significant_hits = data[(data["is_hit"]) & (data["log_pval"] > gwas_level)]
                #for _, row in significant_hits.iterrows():
                #    plt.text(row["adjusted_pos"], row["log_pval"], row["SNP"], fontsize=6, ha='right', rotation=45, color='black')

                # Add GWAS and Bonferroni threshold lines
                plt.axhline(y=gwas_level, color="gray", linestyle="--", linewidth=0.5)
                #plt.axhline(y=bf_level, color="black", linestyle="--", linewidth=0.5)
                chrom_offsets["X"] = chrom_offsets.pop("23")
                plt.xticks(ticks=list(chrom_offsets.values()), labels=list(chrom_offsets.keys()), fontsize=8)
                plt.xlabel("Chromosome")
                plt.ylabel("-log10(p-value)")
                plt.title(f"Manhattan Plot: {pheno}")
                # print(rex)
                plt.savefig(f"{out_file}-l.png", dpi=300, bbox_inches="tight")
                plt.close()
                # print(rex)

                # Q-Q Plot
                # expected = -np.log10(np.linspace(1/len(data), 1, len(data)))
                observed = np.sort(data["log_pval"].values)[::-1]
                expected = -np.log10(stats.uniform.ppf(np.arange(1, len(data)+1) / (len(data)+1)))

                plt.figure(figsize=(4, 4))
                plt.scatter(expected, observed, s=3)
                plt.plot([0, max(expected)], [0, max(expected)], color="gray", linestyle="--")
                plt.axhline(y=gwas_level, color="gray", linestyle="--", linewidth=0.5)
                plt.axhline(y=bf_level, color="black", linestyle="--", linewidth=0.5)
                plt.xlabel("Expected -log10(p-value)")
                plt.ylabel("Observed -log10(p-value)")
                plt.title("Q-Q Plot")
                plt.savefig(f"{out_file}-r.png", dpi=300, bbox_inches="tight")
                plt.close()

                # Combine plots
                from PIL import Image

                manhattan_img = Image.open(f"{out_file}-l.png")
                qq_img = Image.open(f"{out_file}-r.png")

                combined_img = Image.new("RGB", (manhattan_img.width + qq_img.width, max(manhattan_img.height, qq_img.height)))
                combined_img.paste(manhattan_img, (0, 0))
                combined_img.paste(qq_img, (manhattan_img.width, 0))

                combined_img.save(f"{out_file}-f.png")
                combined_img.convert("RGB").save(f"{out_file}-f.pdf", "PDF")
                print(f"Step X took {time.time() - start_time:.2f} seconds")

                # print(rex)
            
    except:
        print ("didn't finish running")
