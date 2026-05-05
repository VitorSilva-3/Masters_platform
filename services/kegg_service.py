
import logging
from Bio.KEGG import REST
from typing import Dict, Any
from utils import load_json_cache, save_json_cache

logger = logging.getLogger(__name__)

class KeggService:
    """Service to interact with the KEGG Enzyme Database with local caching."""

    def __init__(self, cache_file: str = "data/kegg_cache.json"):
        self.cache_file = cache_file
        self.cache = load_json_cache(self.cache_file, service_name = "KeggService")

    def fetch_enzyme_details(self, enzyme_name: str, ec_number: str) -> Dict[str, Any]:
        """Fetches detailed information from KEGG for a specific enzyme."""

        if ec_number in self.cache:
            return self.cache[ec_number]

        entry_id = f"ec:{ec_number}"
        try:
            raw_data = REST.kegg_get(entry_id).read()
            info = self._parse_kegg_text(raw_data, ec_number)
            
            if info and info.get('name'): 
                self.cache[ec_number] = info
                save_json_cache(self.cache_file, self.cache, service_name = "KeggService")
                
            return info
            
        except Exception as e:
            logger.error(f"[KeggService] Error fetching {entry_id}: {e}") 
            return {}

    def _parse_kegg_text(self, raw_text: str, ec_number: str) -> Dict[str, Any]:
        """Parses the raw text response from KEGG into a dictionary."""

        info = {
            'ec_number': ec_number,
            'name': None,
            'definition': None,
            'reaction': None,
            'pathways': []
        }

        for line in raw_text.splitlines():
            if line.startswith("NAME"):
                info['name'] = line.split("NAME")[1].strip().rstrip(";")
            elif line.startswith("DEFINITION"):
                info['definition'] = line.split("DEFINITION")[1].strip()
            elif line.startswith("REACTION"):
                info['reaction'] = line.split("REACTION")[1].strip()
            elif line.startswith("PATHWAY"):
                parts = line.split()
                if len(parts) > 2:
                    code = parts[1]
                    desc = " ".join(parts[2:])
                    info['pathways'].append({'code': code, 'description': desc})

        return info
    
    def save_cache(self):
        """Saves the cache to disk."""
        save_json_cache(self.cache_file, self.cache, service_name = "KeggService")