
import json
import os
import logging
from Bio import Entrez
from typing import Dict, Any, List, Union

logger = logging.getLogger(__name__)

class TaxonomyService:
    """Service to fetch organism lineage from NCBI Taxonomy"""

    def __init__(self, email: str, cache_file: str = "data/taxonomy_cache.json"):
        """Initialize with user email and load local taxonomy cache."""

        self.email = email
        Entrez.email = self.email
        self.cache_file = cache_file
        self.cache = self._load_cache()

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
        Returns a list of dictionaries with 'rank' and 'name' or an error string.
        """
        
        if organism_name in self.cache:
            return self.cache[organism_name]

        try:
            handle = Entrez.esearch(db="taxonomy", term=f"{organism_name}[Scientific Name]")
            record = Entrez.read(handle)
            handle.close()

            if not record["IdList"]:
                result = "Taxonomy ID not found."
            else:
                tax_id = record["IdList"][0]
                handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
                details = Entrez.read(handle)
                handle.close()
                
                if details:
                    lineage_ex = details[0].get("LineageEx", [])
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