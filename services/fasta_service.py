
import os
import time
import logging
import pandas as pd
from urllib.error import HTTPError
from Bio import Entrez
from config import AppConfig

logger = logging.getLogger(__name__)

class FastaService:
    """Class responsible for downloading and batching FASTA sequences from NCBI."""

    def __init__(self, tax_service, email: str, data_filepath: str = "data/enzymes_data.csv", output_dir: str = "data"):
        self.tax_service = tax_service
        self.email = email
        self.data_filepath = data_filepath
        self.output_dir = output_dir
        
        Entrez.email = self.email
        if hasattr(AppConfig, 'NCBI_API_KEY'):
            Entrez.api_key = AppConfig.NCBI_API_KEY
        
        os.makedirs(self.output_dir, exist_ok=True)
        
        self.eukaryote_fasta = os.path.join(self.output_dir, "euk_enzymes.fasta")
        self.prokaryote_fasta = os.path.join(self.output_dir, "prok_enzymes.fasta")

    def _split_ids_by_domain(self, df: pd.DataFrame) -> tuple:
        """Classifies species and splits NCBI IDs into eukaryotes and prokaryotes."""

        logger.info("Classifying species to separate eukaryotes from prokaryotes...")
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
                
        logger.info(f"Classification complete: {len(eukaryote_ids)} eukaryotic enzymes, {len(prokaryote_ids)} prokaryotic enzymes.")
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
                            logger.warning(f"Network error on batch {batch_num} ({type(e).__name__}). Retrying in 5 seconds... (Attempt {attempt + 1}/{max_retries})")
                            time.sleep(5)
                        else:
                            logger.error(f"Failed to download batch {batch_num} after {max_retries} attempts. Error: {e}")
                
                time.sleep(0.5)

    def generate_fasta_files(self, df: pd.DataFrame) -> None:
        """Main method to orchestrate the FASTA generation process."""

        logger.info("Starting FASTA sequence generation process...")
        
        euk_ids, pro_ids = self._split_ids_by_domain(df)
        
        if euk_ids:
            logger.info("Processing eukaryotic enzyme sequences...")
            self._fetch_and_save_batches(euk_ids, self.eukaryote_fasta)
            
        if pro_ids:
            logger.info("Processing prokaryotic enzyme sequences...")
            self._fetch_and_save_batches(pro_ids, self.prokaryote_fasta)
            
        logger.info("FASTA sequence generation successfully completed.")