import pandas as pd
import scanpy as sc
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.inspection import permutation_importance


immune_data = sc.read_h5ad(r"data\H5AD Objects\Immune_Raw_New_Harmonized.h5ad")
immune_data = immune_data[immune_data.obs["cell_type"] != 'Erythrocyte']

# Filter to only Early and Late groups
immune_data = immune_data[immune_data.obs["treatment_group"].isin(["WoundEarly", "WoundLate"])].copy()


ensembl = pd.read_csv("Ensembl_New.txt")

# Save the original ENSEMBL IDs
original_ids = immune_data.var_names.copy()

# Create dictionaries containing salmon, zebrafish and human gene names
salmon_to_gene = dict(zip(ensembl["Gene stable ID"], ensembl["Gene name"]))
zebra_to_gene = dict(zip(ensembl["Gene stable ID"], ensembl["Zebrafish gene name"]))
human_to_gene = dict(zip(ensembl["Gene stable ID"], ensembl["Human gene name"]))
var_names_series = pd.Series(immune_data.var_names, index=immune_data.var_names)

# Map the dictionaries to the Ensembl IDs stored in .var
immune_data.var["salmon_name"] = original_ids.map(salmon_to_gene)
immune_data.var["zebra_name"] = original_ids.map(zebra_to_gene)
immune_data.var["human_name"] = original_ids.map(human_to_gene)
immune_data.var["final_name"] = immune_data.var["salmon_name"].fillna(immune_data.var["zebra_name"]).fillna(immune_data.var["human_name"]).fillna(var_names_series)
# Create a final mapping dictionary
gene_id_to_name = dict(zip(original_ids, immune_data.var["final_name"]))


# Select variable features (DEGs) and plot them
sc.pp.highly_variable_genes(immune_data, n_top_genes=3000)
immune_data = immune_data[:, immune_data.var.highly_variable]


# Random forest process
X = immune_data.X
y = immune_data.obs["treatment_group"].values

# Convert to dense array if sparse
if isinstance(X, np.ndarray) is False:
    X = X.toarray()

label_encoder = LabelEncoder()
y_encoded = label_encoder.fit_transform(y)

# Train/Test Split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y_encoded)

# Train Random Forest Classifier
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_train, y_train)

# Predict on test set
y_pred = rf.predict(X_test)


result = permutation_importance(
    rf, 
    X_test, 
    y_test, 
    n_repeats=5, 
    random_state=42, 
    n_jobs=-1
)

# Organize results in a DataFrame
importance_df = pd.DataFrame({
    "Gene": immune_data.var_names,  # List of gene names from your scRNA-seq data
    "Importance": result.importances_mean
})

importance_df.sort_values(by="Importance", ascending=False, inplace=True)
importance_df = importance_df.head(200)

desired_order = ["Ensembl ID", "Gene Name", "Importance"]
importance_df["Gene Name"] = importance_df["Gene"].map(gene_id_to_name)
importance_df.rename(columns={'Gene': 'Ensembl ID'}, inplace=True)
importance_df = importance_df[desired_order]
importance_df.to_csv("Ranked genes by TG Early-Late only.csv", index=False)

