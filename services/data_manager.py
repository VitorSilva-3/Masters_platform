
import pandas as pd
import os
import time
import logging
import http.client
from Bio import Entrez, SeqIO
from config import AppConfig

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)

class DataManager:
    """
    Manages the local datasets(enzymes and transporters). Data comes from NCBI, filtered by target taxa and specific enzyme/transporter criteria.
    """

    def __init__(self, file_path: str = "data/enzymes_data.csv", data_type: str = "enzyme"):
        self.file_path = file_path
        self.data_type = data_type 

    def load_data(self) -> pd.DataFrame:
        if os.path.exists(self.file_path):
            try:
                df = pd.read_csv(self.file_path)
                return df if not df.empty else pd.DataFrame()
            except Exception:
                return pd.DataFrame()
        return pd.DataFrame()

    def _fetch_batch_with_retry(self, batch_ids: list, max_retries: int = 3) -> list:
        for attempt in range(max_retries):
            try:
                handle = Entrez.efetch(db="protein", id=batch_ids, rettype="gb", retmode="text")
                records = list(SeqIO.parse(handle, "genbank"))
                handle.close()
                return records
            except (http.client.IncompleteRead, Exception) as e:
                wait_time = (attempt + 1) * 5
                logger.warning(f"Network error ({e}). Attempt {attempt+1}/{max_retries}. Retrying in {wait_time}s...")
                time.sleep(wait_time)
        
        logger.error(f"Failed to fetch batch {batch_ids[:3]}... after {max_retries} attempts.")
        return []

    def _parse_and_validate_record(self, record, item_name, info, target_taxa_lower) -> dict:
        desc = record.description.lower()
        
        excluded_list = AppConfig.EXCLUDED_TERMS_ENZYMES if self.data_type == "enzyme" else AppConfig.EXCLUDED_TERMS_TRANSPORTERS

        if any(bad_term in desc for bad_term in excluded_list):
            return None

        raw_organism = record.annotations.get("organism", "Unknown")
        organism_clean = raw_organism.replace("['", "").replace("']", "")
        
        taxonomy_list = record.annotations.get("taxonomy", [])
        taxonomy_lower = [tax.lower() for tax in taxonomy_list]
        is_valid_taxa = any(target in taxonomy_lower for target in target_taxa_lower)
        
        if not is_valid_taxa:
            return None

        if "uncharacterized" in desc:
            status = "🔴 Uncharacterized"
        elif any(term in desc for term in AppConfig.FUTURE_TERMS):
            status = "🟡 Probable/Hypothetical"
        else:
            status = "🟢 Confirmed"

        if self.data_type == "enzyme":
            return {
                "Specie": organism_clean,
                "Enzyme": item_name,
                "EC number": info.ec_number,
                "Target sugar": info.target_sugar,
                "Description": record.description,
                "Status": status,
                "Source ID (NCBI)": record.id 
            }
        elif self.data_type == "transporter":
            return {
                "Specie": organism_clean,
                "Transporter": item_name,
                "TC number": info.tc_number,
                "Family": info.family,
                "Target sugar": info.target_sugar,
                "Description": record.description,
                "Status": status,
                "Source ID (NCBI)": record.id 
            }

    def _save_checkpoint(self, new_records: list):
        if not new_records:
            return

        existing_df = self.load_data()
        new_df = pd.DataFrame(new_records)
        
        if not existing_df.empty:
            final_df = pd.concat([existing_df, new_df], ignore_index=True)
        else:
            final_df = new_df
            
        final_df = final_df.drop_duplicates(subset=['Source ID (NCBI)'])
        
        os.makedirs(os.path.dirname(self.file_path), exist_ok=True)
        final_df.to_csv(self.file_path, index=False)
        logger.debug(f"Checkpoint saved: {len(final_df)} total records now in disk.")

    def build_dataset_offline(self):
        logger.info(f"Starting {self.data_type.upper()} dataset build. Target Taxa: {len(AppConfig.TARGET_TAXA)} groups.")
        
        Entrez.email = AppConfig.EMAIL
        if hasattr(AppConfig, 'NCBI_API_KEY'):
            Entrez.api_key = AppConfig.NCBI_API_KEY

        existing_df = self.load_data()
        existing_ids = set()
        if not existing_df.empty and 'Source ID (NCBI)' in existing_df.columns:
            existing_ids = set(existing_df['Source ID (NCBI)'].dropna().astype(str))
            logger.info(f"Skipping {len(existing_ids)} IDs already present in CSV.")

        taxa_query = " OR ".join([f'"{t}"[Organism]' for t in AppConfig.TARGET_TAXA])
        target_taxa_lower = [t.lower() for t in AppConfig.TARGET_TAXA]
        
        target_dict = AppConfig.ENZYMES if self.data_type == "enzyme" else AppConfig.TRANSPORTERS
        
        for item_name, info in target_dict.items():
            
            if self.data_type == "enzyme":
                logger.info(f"Processing: {item_name} (EC: {info.ec_number})")
                query = f'("{item_name}"[Protein Name] OR "{item_name}"[Title] OR "{info.ec_number}"[EC/RN Number]) AND ({taxa_query})'
            else:
                logger.info(f"Processing: {item_name} (TC: {info.tc_number})")
                query = f'("{item_name}"[Protein Name] OR "{item_name}"[Title] OR "{info.tc_number}"[All Fields]) AND ({taxa_query})'
            
            try:
                handle = Entrez.esearch(db="protein", term=query, retmax=5000)
                search_results = Entrez.read(handle)
                handle.close()
                
                all_ids = search_results.get("IdList", [])
                ids_to_fetch = [ncbi_id for ncbi_id in all_ids if ncbi_id not in existing_ids]
                
                if not ids_to_fetch:
                    logger.info(f"No new records for {item_name}.")
                    continue

                logger.info(f"Found {len(all_ids)} total. Fetching {len(ids_to_fetch)} new records...")

                batch_size = 100
                for i in range(0, len(ids_to_fetch), batch_size):
                    batch_ids = ids_to_fetch[i:i+batch_size]
                    records = self._fetch_batch_with_retry(batch_ids)
                    
                    batch_valid_records = []
                    for record in records:
                        data = self._parse_and_validate_record(record, item_name, info, target_taxa_lower)
                        if data:
                            batch_valid_records.append(data)
                    
                    if batch_valid_records:
                        self._save_checkpoint(batch_valid_records)
                        for r in batch_valid_records:
                            existing_ids.add(str(r["Source ID (NCBI)"]))
                    
                    time.sleep(0.5 if Entrez.api_key else 1.0) 

            except Exception as e:
                logger.error(f"Critical error during {item_name} processing: {e}")

        logger.info(f"Build process finished. Final dataset saved at {self.file_path} with {len(existing_ids)} unique records.")

if __name__ == "__main__":
    logger.info("Beginning dataset construction for enzymes.")
    manager_enzymes = DataManager(file_path="data/enzymes_data.csv", data_type="enzyme")
    manager_enzymes.build_dataset_offline()
    
    logger.info("Beginning dataset construction for transporters.")
    manager_transporters = DataManager(file_path="data/transporters_data.csv", data_type="transporter")
    manager_transporters.build_dataset_offline()