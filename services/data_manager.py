
import pandas as pd
import os
import time
import logging
from Bio import Entrez, SeqIO
from config import AppConfig, FUTURE_TERMS

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)

EXCLUDED_TERMS = [
    "inhibitor", "regulator", "activator", "repressor", 
    "receptor", "transcription factor", "fingers", 
    "binding protein", "domain-containing",
    "dna", "rna", "trna", "mrna", "ribosomal", "ribosome",
    "recombinase", "integrase", "transposase", "nuclease",
    "polymerase", "helicase", "chromosome", "plasmid",
    "protein kinase", "histidine kinase", "tyrosine kinase", 
    "serine/threonine", "signal transduction",
    "synthase", "biosynthesis", "assembly"
]

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
        
        taxa_query_parts = [f'"{t}"[Organism]' for t in AppConfig.TARGET_TAXA]
        taxa_query_full = " OR ".join(taxa_query_parts)
        
        all_records = []
        
        for enzyme_name, info in AppConfig.ENZYMES.items():
            logger.info(f"Fetching data for: {enzyme_name} (EC: {info.ec_number})...")
            
            term_query = f'("{enzyme_name}"[Protein Name] OR "{enzyme_name}"[Title] OR "{info.ec_number}"[EC/RN Number])'
            full_query = f'({term_query}) AND ({taxa_query_full})'
            
            try:
                handle = Entrez.esearch(db="protein", term=full_query, retmax=2000)
                search_results = Entrez.read(handle)
                handle.close()
                
                id_list = search_results["IdList"]
                if not id_list:
                    logger.warning(f"No records found for {enzyme_name}.")
                    continue
                
                batch_size = 100
                for i in range(0, len(id_list), batch_size):
                    batch_ids = id_list[i:i+batch_size]
                    
                    try:
                        handle = Entrez.efetch(db="protein", id=batch_ids, rettype="gb", retmode="text")
                        records = SeqIO.parse(handle, "genbank")
                        
                        for record in records:
                            desc = record.description.lower()
                            
                            if any(bad_term in desc for bad_term in EXCLUDED_TERMS):
                                continue

                            raw_organism = record.annotations.get("organism", "Unknown")
                            organism_clean = raw_organism.replace("['", "").replace("']", "")

                            if "uncharacterized" in desc:
                                status = "🔴 Uncharacterized"
                            elif any(term in desc for term in FUTURE_TERMS):
                                status = "🟡 Probable/Hypothetical"
                            else:
                                status = "🟢 Confirmed"

                            all_records.append({
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

        if all_records:
            df = pd.DataFrame(all_records)
            
            logger.info(f"Records before cleaning (duplicate removal): {len(df)}")
            
            df = df.drop_duplicates(subset=['Specie', 'Enzyme', 'Description'])
            
            os.makedirs(os.path.dirname(self.file_path), exist_ok=True)
            df.to_csv(self.file_path, index=False)
            logger.info(f"Success! Saved {len(df)} records to '{self.file_path}'.")
        else:
            logger.warning("No records were found in the global search.")

if __name__ == "__main__":
    manager = DataManager()
    manager.build_dataset_offline()