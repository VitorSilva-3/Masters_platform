
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
    Manages the local dataset: building it from NCBI and loading it for the app.
    """

    def __init__(self, file_path: str = "data/enzymes_data.csv"):
        self.file_path = file_path

    def load_data(self) -> pd.DataFrame:
        """Loads the existing CSV dataset."""

        if os.path.exists(self.file_path):
            try:
                df = pd.read_csv(self.file_path)
                return df if not df.empty else pd.DataFrame()
            except Exception:
                return pd.DataFrame()
        return pd.DataFrame()

    def _fetch_batch_with_retry(self, batch_ids: list, max_retries: int = 3) -> list:
        """Fetches a batch of records from NCBI."""

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

    def _parse_and_validate_record(self, record, enzyme_name, info, target_taxa_lower) -> dict:
        """Parses a GenBank record and validates it against excluded terms and target taxa.
        Returns a dictionary with data if valid, else None."""

        desc = record.description.lower()
        
        if any(bad_term in desc for bad_term in AppConfig.EXCLUDED_TERMS):
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

        return {
            "Specie": organism_clean,
            "Enzyme": enzyme_name,
            "EC number": info.ec_number,
            "Target sugar": info.target_sugar,
            "Description": record.description,
            "Status": status,
            "Source ID (NCBI)": record.id 
        }

    def _save_checkpoint(self, new_records: list):
        """Merges new records with existing file and saves immediately to prevent data loss."""

        if not new_records:
            return

        existing_df = self.load_data()
        new_df = pd.DataFrame(new_records)
        
        if not existing_df.empty:
            final_df = pd.concat([existing_df, new_df], ignore_index=True)
        else:
            final_df = new_df
            
        final_df = final_df.drop_duplicates(subset=['Specie', 'Enzyme', 'Description'])
        
        os.makedirs(os.path.dirname(self.file_path), exist_ok=True)
        final_df.to_csv(self.file_path, index=False)
        logger.debug(f"Checkpoint saved: {len(final_df)} total records now in disk.")

    def build_dataset_offline(self):
        """Orchestrates the data collection with API Key and progressive saving."""

        logger.info(f"Starting dataset build. Target Taxa: {len(AppConfig.TARGET_TAXA)} groups.")
        
        Entrez.email = AppConfig.EMAIL
        if hasattr(AppConfig, 'NCBI_API_KEY'):
            Entrez.api_key = AppConfig.NCBI_API_KEY
            logger.info("Using NCBI API Key for higher throughput.")

        existing_df = self.load_data()
        existing_ids = set()
        if not existing_df.empty:
            existing_ids = set(existing_df['Source ID (NCBI)'].dropna().astype(str))
            logger.info(f"Skipping {len(existing_ids)} IDs already present in CSV.")

        taxa_query = " OR ".join([f'"{t}"[Organism]' for t in AppConfig.TARGET_TAXA])
        target_taxa_lower = [t.lower() for t in AppConfig.TARGET_TAXA]
        
        for enzyme_name, info in AppConfig.ENZYMES.items():
            logger.info(f"Processing: {enzyme_name} (EC: {info.ec_number})")
            
            query = f'("{enzyme_name}"[Protein Name] OR "{enzyme_name}"[Title] OR "{info.ec_number}"[EC/RN Number]) AND ({taxa_query})'
            
            try:
                handle = Entrez.esearch(db="protein", term=query, retmax=5000)
                search_results = Entrez.read(handle)
                handle.close()
                
                all_ids = search_results.get("IdList", [])
                ids_to_fetch = [ncbi_id for ncbi_id in all_ids if ncbi_id not in existing_ids]
                
                if not ids_to_fetch:
                    logger.info(f"No new records for {enzyme_name}.")
                    continue

                logger.info(f"Found {len(all_ids)} total. Fetching {len(ids_to_fetch)} new records...")

                batch_size = 100
                enzyme_new_records = []
                
                for i in range(0, len(ids_to_fetch), batch_size):
                    batch_ids = ids_to_fetch[i:i+batch_size]
                    
                    records = self._fetch_batch_with_retry(batch_ids)
                    
                    batch_valid_records = []
                    for record in records:
                        data = self._parse_and_validate_record(record, enzyme_name, info, target_taxa_lower)
                        if data:
                            batch_valid_records.append(data)
                    
                    if batch_valid_records:
                        self._save_checkpoint(batch_valid_records)
                        for r in batch_valid_records:
                            existing_ids.add(str(r["Source ID (NCBI)"]))
                    
                    time.sleep(0.5 if Entrez.api_key else 1.0) 

            except Exception as e:
                logger.error(f"Critical error during {enzyme_name} processing: {e}")

        logger.info(f"Build process finished. Final dataset saved at {self.file_path} with {len(existing_ids)} unique records.")

if __name__ == "__main__":
    manager = DataManager()
    manager.build_dataset_offline()