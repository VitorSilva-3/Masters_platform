
import os
import time
import logging
import pandas as pd
from Bio import Entrez
from config import AppConfig

logger = logging.getLogger(__name__)

class FastaService:
    """Class responsible for downloading and batching FASTA sequences from NCBI."""

    def __init__(self, tax_service, email: str, output_dir: str = "data"):
        self.tax_service = tax_service
        self.email = email
        self.output_dir = output_dir
        
        Entrez.email = self.email
        if hasattr(AppConfig, 'NCBI_API_KEY'):
            Entrez.api_key = AppConfig.NCBI_API_KEY
        
        os.makedirs(self.output_dir, exist_ok=True)

    def _split_ids_by_domain(self, df: pd.DataFrame) -> tuple:
        """Classifies species and splits NCBI IDs into eukaryotes and prokaryotes."""

        eukaryote_ids = []
        prokaryote_ids = []

        df_clean = df.dropna(subset=['Source ID (NCBI)'])

        for _, row in df_clean.iterrows():
            specie = row['Specie']
            ncbi_id = row['Source ID (NCBI)']
            
            taxonomy_data = self.tax_service.cache.get(specie, [])
            
            if isinstance(taxonomy_data, list):
                lineage = [node.get('name', '') for node in taxonomy_data]
            else:
                lineage = []
            
            if "Cyanophyceae" in lineage or "Cyanobacteria" in lineage:
                prokaryote_ids.append(ncbi_id)
            else:
                eukaryote_ids.append(ncbi_id)
                
        return eukaryote_ids, prokaryote_ids

    def _get_existing_ids(self, filepath: str) -> set:
        """Reads an existing FASTA file and extracts all NCBI IDs to avoid redundant downloads."""

        existing_ids = set()
        if os.path.exists(filepath):
            with open(filepath, "r") as file:
                for line in file:
                    if line.startswith(">"):
                        full_id = line.split()[0][1:]
                        existing_ids.add(full_id)
                        existing_ids.add(full_id.split('.')[0])
        return existing_ids

    def _fetch_and_save_batches(self, id_list: list, output_file: str, batch_size: int = 200) -> None:
        """Fetches sequences from NCBI in batches and saves them to a FASTA file."""

        if not id_list:
            logger.warning(f"No IDs provided for {output_file}. Skipping download.")
            return

        existing_ids = self._get_existing_ids(output_file)
        new_ids = list(set([ncbi_id for ncbi_id in id_list if ncbi_id not in existing_ids and ncbi_id.split('.')[0] not in existing_ids]))

        if not new_ids:
            logger.info(f"All sequences for '{os.path.basename(output_file)}' are already downloaded. Skipping.")
            return

        total_batches = (len(new_ids) // batch_size) + 1
        logger.info(f"Downloading {len(new_ids)} new sequences in {total_batches} batches to '{output_file}'...")

        with open(output_file, "a") as file:
            for i in range(0, len(new_ids), batch_size):
                batch = new_ids[i:i + batch_size]
                request_ids = ",".join(batch)
                batch_num = (i // batch_size) + 1
                
                logger.info(f"Processing batch {batch_num}/{total_batches} ({len(batch)} sequences)...")
                
                max_retries = 3
                for attempt in range(max_retries):
                    try:
                        handle = Entrez.efetch(db="protein", id=request_ids, rettype="fasta", retmode="text")
                        records = handle.read()
                        
                        if records.strip(): 
                            file.write(records)
                            if not records.endswith('\n'):
                                file.write('\n') 
                                
                        handle.close()
                        break 
                    except Exception as e:  
                        if attempt < max_retries - 1:
                            logger.warning(f"Network error on batch {batch_num}. Retrying in 5s... (Attempt {attempt + 1}/{max_retries})")
                            time.sleep(5)
                        else:
                            logger.error(f"Failed to download batch {batch_num} after {max_retries} attempts. Error: {e}")
                
                time.sleep(0.5)

    def generate_fasta_files(self, df: pd.DataFrame, dataset_name: str) -> None:
        """Generates the FASTA files for a specific dataset dynamically."""

        logger.info(f"Starting FASTA sequence generation process for: {dataset_name.upper()}")
        
        eukaryote_fasta = os.path.join(self.output_dir, f"euk_{dataset_name}.fasta")
        prokaryote_fasta = os.path.join(self.output_dir, f"prok_{dataset_name}.fasta")
        
        euk_ids, pro_ids = self._split_ids_by_domain(df)
        logger.info(f"Classification complete: {len(euk_ids)} eukaryotic, {len(pro_ids)} prokaryotic sequences.")
        
        if euk_ids:
            self._fetch_and_save_batches(euk_ids, eukaryote_fasta)
            
        if pro_ids:
            self._fetch_and_save_batches(pro_ids, prokaryote_fasta)
            
        logger.info(f"FASTA generation for {dataset_name} completed.")

    def process_all_datasets(self, datasets: dict):
        """Orchestrates the FASTA generation for multiple datasets automatically."""

        for ds_name, df in datasets.items():
            if not df.empty:
                self.generate_fasta_files(df, dataset_name=ds_name)
            else:
                logger.warning(f"Dataset '{ds_name}' is empty. Skipping FASTA generation.")