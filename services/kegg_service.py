
import json
import os
from Bio.KEGG import REST
from typing import Dict, Any

class KeggService:
    """Service to interact with the KEGG Enzyme Database with local caching."""

    def __init__(self, cache_file: str = "data/kegg_cache.json"):
        self.cache_file = cache_file
        self.cache = self._load_cache()

    def _load_cache(self) -> Dict[str, Any]:
        """Loads the cache from a JSON file if it exists, otherwise returns an empty dictionary."""

        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception as e:
                print(f"[KeggService] Error loading cache: {e}")
        return {}

    def _save_cache(self):
        """Saves the current cache dictionary to a JSON file."""

        try:
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump(self.cache, f, indent=4)
        except Exception as e:
            print(f"[KeggService] Error saving cache: {e}")

    def fetch_enzyme_details(self, enzyme_name: str, ec_number: str) -> Dict[str, Any]:
        """
        Fetches detailed information from KEGG for a specific enzyme.
        """

        if ec_number in self.cache:
            return self.cache[ec_number]

        entry_id = f"ec:{ec_number}"
        try:
            raw_data = REST.kegg_get(entry_id).read()
            info = self._parse_kegg_text(raw_data, ec_number)
            
            if info and info.get('name'): 
                self.cache[ec_number] = info
                self._save_cache()
                
            return info
            
        except Exception as e:
            print(f"[KeggService] Error fetching {entry_id}: {e}")
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