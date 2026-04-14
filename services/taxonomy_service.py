
import json
import os
import logging
from Bio import Entrez
from typing import Dict, Any, List, Union
from config import AppConfig  

logger = logging.getLogger(__name__)

class TaxonomyService:
    """Service to fetch organism lineage from NCBI Taxonomy"""

    def __init__(self, email: str, cache_file: str = "data/taxonomy_cache.json"):
        """Initialize with user email and load local taxonomy cache."""

        self.email = email
        Entrez.email = self.email
        self.cache_file = cache_file
        self.cache = self._load_cache()

        if hasattr(AppConfig, 'NCBI_API_KEY'):
            Entrez.api_key = AppConfig.NCBI_API_KEY

    def _load_cache(self) -> Dict[str, Any]:
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception as e:
                logger.error(f"[TaxonomyService] Error loading cache: {e}")
        return {}

    def _save_cache(self):
        try:
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump(self.cache, f, indent=4)
        except Exception as e:
            logger.error(f"[TaxonomyService] Error saving cache: {e}")

    def fetch_taxonomy_lineage(self, organism_name: str) -> Union[List[Dict[str, str]], str]:
        """
        Fetches the full taxonomy lineage for a specific organism.
        """
        
        if organism_name in self.cache:
            return self.cache[organism_name]

        try:
            handle = Entrez.esearch(db="taxonomy", term=f"{organism_name}[Scientific Name]")
            record = Entrez.read(handle)
            handle.close()

            id_list = record.get("IdList", [])

            if not id_list:
                result = "Taxonomy ID not found."
            else:
                handle = Entrez.efetch(db="taxonomy", id=",".join(id_list), retmode="xml")
                details_list = Entrez.read(handle)
                handle.close()
                
                chosen_details = None

                if len(details_list) == 1:
                    chosen_details = details_list[0]
                else:
                    target_taxa_lower = [t.lower() for t in AppConfig.TARGET_TAXA]
                    
                    for details in details_list:
                        lineage_str = details.get("Lineage", "").lower()
                        
                        if any(taxa in lineage_str for taxa in target_taxa_lower):
                            chosen_details = details
                            logger.info(f"Resolved homonym for '{organism_name}': Chose correct TaxID {details.get('TaxId')}")
                            break
                    
                    if not chosen_details:
                        chosen_details = details_list[0]

                if chosen_details:
                    lineage_ex = chosen_details.get("LineageEx", [])
                    if lineage_ex:
                        result = []
                        for node in lineage_ex:
                            result.append({
                                "rank": node.get("Rank", "no rank"),
                                "name": node.get("ScientificName", "Unknown")
                            })
                    else:
                        result = "Lineage not available"
                else:
                    result = "No details found."

            self.cache[organism_name] = result
            self._save_cache()
            
            return result
            
        except Exception as e:
            logger.error(f"[TaxonomyService] Error fetching taxonomy for {organism_name}: {e}")
            return f"Error fetching taxonomy: {str(e)}"