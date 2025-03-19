# %%
import pandas as pd
import requests
from bs4 import BeautifulSoup
import re
import time
#%%
def get_taxon_id(species_name):
    """Retrieve taxon ID using the TimeTree API if missing."""
    url = f"https://www.timetree.org/ajax/names/{species_name}/Homo sapiens"
    response = requests.get(url)
    soup = BeautifulSoup(response.text, "html.parser")

    taxon_option = soup.find("select", {"id": "pairwise-resolve-taxon-a1"})
    if taxon_option:
        return int(taxon_option.find("option")["value"])  # Convert to integer
    
    print(f"❌ Warning: Could not resolve taxon ID for {species_name}")
    return None  # Failsafe: Return None if unresolved

def extract_divergence_time(taxon_a_id, taxon_b_id):
    """Fetch divergence time using taxon IDs."""
    try:
        # Skip if taxon IDs are missing
        if not taxon_a_id or not taxon_b_id:
            return None, None, None  

        url = f"https://www.timetree.org/ajax/pairwise/{int(taxon_a_id)}/{int(taxon_b_id)}"
        response = requests.get(url)
        soup = BeautifulSoup(response.text, "html.parser")

        # Extract median time
        median_text = soup.find("text", string=re.compile(r"Median Time:"))
        median_time = None
        if median_text:
            median_time_tag = median_text.find_next_sibling("text")
            median_time = float(median_time_tag.text.replace(" MYA", "")) if median_time_tag else None

        # Extract confidence interval (CI)
        ci_min, ci_max = None, None
        ci_text = soup.find("text", string=re.compile(r"CI:"))
        if ci_text:
            ci_match = re.search(r"CI: \((\d+\.\d+) - (\d+\.\d+) MYA\)", ci_text.text)
            if ci_match:
                ci_min, ci_max = float(ci_match.group(1)), float(ci_match.group(2))

        return median_time, ci_min, ci_max

    except Exception as e:
        print(f"❌ Error fetching data for taxon {taxon_a_id} vs {taxon_b_id}: {e}")
        return None, None, None  # Failsafe: Return None in case of failure
#%%
# === Load the table ===
df = pd.read_csv("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species_info447_pre.txt", sep="\t", index_col=0)

# 🔹 Fix column name with non-breaking space
df.columns = df.columns.str.replace("\xa0", " ")  # Rename 'taxon\xa0id' → 'taxon id'
#%%
# Target species taxon ID (Homo sapiens)
human_taxon_id = 9606  # Convert to integer directly

# Ensure taxon_id column exists and convert it to integers
df["taxon id"] = df["taxon id"].apply(lambda x: int(float(x)) if pd.notna(x) else None)

# Create new columns for results
df["median_time"] = None
df["ci_min"] = None
df["ci_max"] = None

# === Process Each Species ===
for index, row in df.iterrows():
    species_name = row["scientific name"].strip()
    taxon_id = row["taxon id"]

    # Resolve taxon ID if missing
    if pd.isna(taxon_id) or taxon_id is None:
        taxon_id = get_taxon_id(species_name)
        df.at[index, "taxon id"] = taxon_id  # Store resolved taxon ID

    if taxon_id:  # Only fetch divergence time if taxon ID is available
        print(f"🔍 Fetching data for {species_name} (Taxon ID: {taxon_id}) vs Homo sapiens...")
        median_time, ci_min, ci_max = extract_divergence_time(human_taxon_id, int(taxon_id))

        # Store results in DataFrame
        df.at[index, "median_time"] = median_time
        df.at[index, "ci_min"] = ci_min
        df.at[index, "ci_max"] = ci_max

        time.sleep(1)  # Prevent server overload

# === Save the updated table ===
df.to_csv("updated_species_table_taxon_id.csv", index=False)
print("✅ Updated table saved as 'updated_species_table.csv'!")
# %%
df.to_csv("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species447_info_b.txt", sep='\t')
# %%
df_447 = pd.read_csv("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species447_info_b.txt", sep="\t", index_col=0)
df_241 = pd.read_csv("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species241_info.tsv", sep="\t", index_col=0)

assembly_info = pd.read_csv("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/PRJEB67744_AssemblyDetails.txt", sep="\t", index_col=False)
#%%
# Function to extract the name outside parentheses
def extract_name_outside_parentheses(name):
    match = re.match(r"^[^(]+", name)  # Extract everything before the first '('
    if match:
        return match.group(0).strip()  # Strip any extra whitespace
    return name.strip()
# Function to extract the name inside parentheses (if it exists)
def extract_name_in_parentheses(name):
    match = re.search(r"\((.*?)\)", name)  # Extract the name inside parentheses
    if match:
        return match.group(1).strip()
    return None
# %%
# 🔹 Fix column names if they contain non-breaking spaces
df_447.columns = df_447.columns.str.replace("\xa0", " ")
df_241.columns = df_241.columns.str.replace("\xa0", " ")

# Create a dictionary mapping scientific name → Accession from 241-species table
species_to_accession = assembly_info.set_index("Taxonomy")["# Assembly"].to_dict()

# Track progress
missing_count = 0
updated_count = 0
not_found_species = []
# Loop through 447-species table and update missing Accessions
for index, row in df_447.iterrows():
    if pd.isna(row["Accession"]) or row["Accession"] == "":
        scientific_name = row["scientific name"].strip()  # First, try name outside parentheses
        missing_count += 1  # Count missing accessions

        # Try matching with the name outside the parentheses
        scientific_name_outside = extract_name_outside_parentheses(scientific_name)
        print(scientific_name_outside)
        if scientific_name_outside in species_to_accession:
            df_447.at[index, "Accession"] = species_to_accession[scientific_name_outside]  # Update Accession
            updated_count += 1
        else:
            # If no match, try the name inside the parentheses (if it exists)
            scientific_name_in_parentheses = extract_name_in_parentheses(scientific_name)
            if scientific_name_in_parentheses and scientific_name_in_parentheses in species_to_accession:
                df_447.at[index, "Accession"] = species_to_accession[scientific_name_in_parentheses]  # Update Accession
                updated_count += 1
            else:
                not_found_species.append(scientific_name)  # Log species not found in 241-species table

# Save the updated table
df_447.to_csv("species_table_447_fixed.csv", index=False)

# Print status report
print("✅ Process completed!")
print(f"🔍 Total missing accessions found: {missing_count}")
print(f"✅ Successfully updated: {updated_count}")
print(f"⚠️ Not found in 241-species table: {len(not_found_species)} species")

# Show a few missing ones for manual checking
if not_found_species:
    print(f"⚠️ Example species not found: {not_found_species[:5]} (showing up to 5)")
    
print("✅ Updated table saved as 'species_table_447_fixed.csv'!")
# %%
df_447.to_csv("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/zoonomia_divergence_ref_table/species447_info_b_partial_fix_a.txt", sep='\t')
# %%
# 🔹 Fix column names if they contain non-breaking spaces
df_447.columns = df_447.columns.str.replace("\xa0", " ")
df_241.columns = df_241.columns.str.replace("\xa0", " ")

# Create a dictionary mapping scientific name → Accession from 241-species table
species_to_accession = df_241.set_index("Species")["Accession"].to_dict()

# Track progress
missing_count = 0
updated_count = 0
not_found_species = []
# Loop through 447-species table and update missing Accessions
for index, row in df_447.iterrows():
    if pd.isna(row["Accession"]) or row["Accession"] == "":
        scientific_name = row["scientific name"].strip()  # First, try name outside parentheses
        missing_count += 1  # Count missing accessions

        # Try matching with the name outside the parentheses
        scientific_name_outside = extract_name_outside_parentheses(scientific_name)
        if scientific_name_outside in species_to_accession:
            df_447.at[index, "Accession"] = species_to_accession[scientific_name_outside]  # Update Accession
            updated_count += 1
        else:
            # If no match, try the name inside the parentheses (if it exists)
            scientific_name_in_parentheses = extract_name_in_parentheses(scientific_name)
            if scientific_name_in_parentheses and scientific_name_in_parentheses in species_to_accession:
                df_447.at[index, "Accession"] = species_to_accession[scientific_name_in_parentheses]  # Update Accession
                updated_count += 1
            else:
                not_found_species.append(scientific_name)  # Log species not found in 241-species table

# Save the updated table
df_447.to_csv("species_table_447_fixed.csv", index=False)

# Print status report
print("✅ Process completed!")
print(f"🔍 Total missing accessions found: {missing_count}")
print(f"✅ Successfully updated: {updated_count}")
print(f"⚠️ Not found in 241-species table: {len(not_found_species)} species")

# Show a few missing ones for manual checking
if not_found_species:
    print(f"⚠️ Example species not found: {not_found_species[:5]} (showing up to 5)")
    
print("✅ Updated table saved as 'species_table_447_fixed.csv'!")
#%%


# %%
