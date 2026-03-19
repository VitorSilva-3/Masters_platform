
import pandas as pd
import os
import time
import logging
from Bio import Entrez, SeqIO
from config import AppConfig

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)


class DataManager:
    """Manages the local dataset: building it from NCBI and loading it for the app."""

    def __init__(self, file_path: str = "data/enzymes_data.csv"):
        self.file_path = file_path

    def load_data(self) -> pd.DataFrame:
        if os.path.exists(self.file_path):
            try:
                df = pd.read_csv(self.file_path)
                if df.empty: return pd.DataFrame()
                return df
            except: return pd.DataFrame()
        return pd.DataFrame()

    def build_dataset_offline(self):
        logger.info(f"Starting dataset build. Target Taxa: {len(AppConfig.TARGET_TAXA)} groups.")
        Entrez.email = AppConfig.EMAIL
        
        existing_df = self.load_data()
        existing_ids = set()
        if not existing_df.empty and 'Source ID (NCBI)' in existing_df.columns:
            existing_ids = set(existing_df['Source ID (NCBI)'].dropna().astype(str))
            logger.info(f"Found existing dataset with {len(existing_df)} records.")
            logger.info(f"Loaded {len(existing_ids)} unique NCBI Source IDs to skip.")
        else:
            logger.info("No existing dataset found or missing Source ID column. Full download will proceed.")

        taxa_query_parts = [f'"{t}"[Organism]' for t in AppConfig.TARGET_TAXA]
        taxa_query_full = " OR ".join(taxa_query_parts)
        
        target_taxa_lower = [t.lower() for t in AppConfig.TARGET_TAXA]
        
        new_records = []
        
        for enzyme_name, info in AppConfig.ENZYMES.items():
            logger.info(f"Fetching data for: {enzyme_name} (EC: {info.ec_number})...")
            
            term_query = f'("{enzyme_name}"[Protein Name] OR "{enzyme_name}"[Title] OR "{info.ec_number}"[EC/RN Number])'
            full_query = f'({term_query}) AND ({taxa_query_full})'
            
            try:
                handle = Entrez.esearch(db="protein", term=full_query, retmax=5000)
                search_results = Entrez.read(handle)
                handle.close()
                
                all_found_ids = search_results.get("IdList", [])
                
                if not all_found_ids:
                    logger.warning(f"No records found for {enzyme_name}.")
                    continue
                
                ids_to_fetch = [ncbi_id for ncbi_id in all_found_ids if ncbi_id not in existing_ids]
                
                skipped_count = len(all_found_ids) - len(ids_to_fetch)
                logger.info(f"Total found: {len(all_found_ids)}. Skipping {skipped_count} already saved. Fetching {len(ids_to_fetch)} new records.")

                if not ids_to_fetch:
                    continue
                
                batch_size = 100
                for i in range(0, len(ids_to_fetch), batch_size):
                    batch_ids = ids_to_fetch[i:i+batch_size]
                    
                    try:
                        handle = Entrez.efetch(db="protein", id=batch_ids, rettype="gb", retmode="text")
                        records = SeqIO.parse(handle, "genbank")
                        
                        for record in records:
                            desc = record.description.lower()
                            
                            if any(bad_term in desc for bad_term in AppConfig.EXCLUDED_TERMS):
                                continue

                            raw_organism = record.annotations.get("organism", "Unknown")
                            organism_clean = raw_organism.replace("['", "").replace("']", "")
                            
                            taxonomy_list = record.annotations.get("taxonomy", [])
                            taxonomy_lower = [tax.lower() for tax in taxonomy_list]
                            
                            is_valid_taxa = any(target in taxonomy_lower for target in target_taxa_lower)
                            
                            if not is_valid_taxa:
                                logger.debug(f"Skipping '{organism_clean}' - Lineage does not match Target Taxa.")
                                continue 

                            if "uncharacterized" in desc:
                                status = "🔴 Uncharacterized"
                            elif any(term in desc for term in AppConfig.FUTURE_TERMS):
                                status = "🟡 Probable/Hypothetical"
                            else:
                                status = "🟢 Confirmed"

                            new_records.append({
                                "Specie": organism_clean,
                                "Enzyme": enzyme_name,
                                "EC number": info.ec_number,
                                "Target sugar": info.target_sugar,
                                "Description": record.description,
                                "Status": status,
                                "Source ID (NCBI)": record.id 
                            })
                        handle.close()
                    except Exception as e:
                        logger.error(f"Error in batch for {enzyme_name}: {e}")
                        continue
                    
                    time.sleep(1) 

            except Exception as e:
                logger.error(f"Error processing {enzyme_name}: {e}")

        if new_records or not existing_df.empty:
            new_df = pd.DataFrame(new_records)
            
            if not existing_df.empty and not new_df.empty:
                final_df = pd.concat([existing_df, new_df], ignore_index=True)
                logger.info(f"Merged {len(existing_df)} old records with {len(new_df)} new records.")
            elif not new_df.empty:
                final_df = new_df
            else:
                final_df = existing_df
            
            logger.info(f"Records before cleaning (duplicate removal): {len(final_df)}")
            
            final_df = final_df.drop_duplicates(subset=['Specie', 'Enzyme', 'Description'])
            
            os.makedirs(os.path.dirname(self.file_path), exist_ok=True)
            final_df.to_csv(self.file_path, index=False)
            logger.info(f"Success! Saved total of {len(final_df)} records to '{self.file_path}'.")
        else:
            logger.warning("No records were found in the global search and no existing data to save.")

if __name__ == "__main__":
    manager = DataManager()
    manager.build_dataset_offline()