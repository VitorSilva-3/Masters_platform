
import requests
import json
import os
import logging
from typing import Dict, Any, List

logger = logging.getLogger(__name__)

class UniprotService:
    """Service to interact with the UniProtKB Database via API with local caching."""

    BASE_URL = "https://rest.uniprot.org/uniprotkb/search"

    def __init__(self, cache_file: str = "data/uniprot_cache.json"):
        """Initialize and load local UniProt cache."""
        self.cache_file = cache_file
        self.cache = self._load_cache()

    def _load_cache(self) -> Dict[str, Any]:
        """Loads the JSON cache from disk."""
        
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception as e:
                logger.error(f"[UniprotService] Error loading cache: {e}")
        return {}

    def _save_cache(self):
        """Saves the current cache dictionary to disk."""

        try:
            os.makedirs(os.path.dirname(self.cache_file), exist_ok=True)
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump(self.cache, f, indent=4)
        except Exception as e:
            logger.error(f"[UniprotService] Error saving cache: {e}")

    def fetch_enzyme_data(self, enzyme_name: str, ec_number: str, max_entries: int = 30) -> Dict[str, Any]:
        """
        Searches UniProt for general protein properties matching the EC number.
        """

        if ec_number in self.cache:
            return self.cache[ec_number]

        params = {
            'query': f'ec:{ec_number}',
            'format': 'json',
            'size': max_entries,
            'fields': 'accession,protein_name,cc_function,cc_catalytic_activity,cc_pathway,cc_cofactor,cc_subunit'
        }
        
        try:
            response = requests.get(self.BASE_URL, params=params, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            results = data.get('results', [])
            if not results:
                self.cache[ec_number] = {}
                self._save_cache()
                return {}
            
            summary = self._summarize_results(enzyme_name, ec_number, results)
            self.cache[ec_number] = summary
            self._save_cache()
            return summary
            
        except Exception as e:
            logger.error(f"[UniprotService] Error fetching general data for {enzyme_name} (EC {ec_number}): {e}")
            return {}

    def fetch_species_link(self, ec_number: str, species_name: str) -> Dict[str, str]:
        """
        Searches UniProt specifically for an entry matching the EC number AND Species.
        """

        cache_key = f"{ec_number}_{species_name}_link"
        if cache_key in self.cache:
            return self.cache[cache_key]

        params = {
            'query': f'ec:{ec_number} AND organism_name:"{species_name}"',
            'format': 'json',
            'size': 1, 
            'fields': 'accession'
        }
        
        try:
            response = requests.get(self.BASE_URL, params=params, timeout=10)
            if response.status_code == 200:
                results = response.json().get('results', [])
                if results:
                    accession = results[0].get('primaryAccession')
                    if accession:
                        result = {
                            "status": "Found",
                            "accession": accession,
                            "url": f"https://www.uniprot.org/uniprotkb/{accession}/entry"
                        }
                        self.cache[cache_key] = result
                        self._save_cache()
                        return result
        except Exception as e:
            logger.error(f"[UniprotService] Error fetching link for {species_name}: {e}")

        fallback_url = f"https://www.uniprot.org/uniprotkb?query=ec:{ec_number}+AND+organism_name:\"{species_name}\""
        result = {
            "status": "Not Found",
            "accession": "-",
            "url": fallback_url
        }
        self.cache[cache_key] = result
        self._save_cache()
        return result

    def _summarize_results(self, name: str, ec: str, results: List[dict]) -> Dict[str, Any]:
        """Aggregates general data."""

        functions = set()
        catalytic_activities = []
        pathways = set()
        cofactors = set()
        subunits = set()

        for entry in results:
            for comment in entry.get('comments', []):
                ctype = comment.get('commentType')
                
                if ctype == 'FUNCTION':
                    for txt in comment.get('texts', []):
                        functions.add(txt.get('value'))
                elif ctype == 'CATALYTIC ACTIVITY':
                    reaction = comment.get('reaction', {}).get('name')
                    if reaction: catalytic_activities.append(reaction)
                elif ctype == 'PATHWAY':
                    for txt in comment.get('texts', []):
                        pathways.add(txt.get('value'))
                elif ctype == 'COFACTOR':
                    for txt in comment.get('texts', []):
                        cofactors.add(txt.get('value'))
                elif ctype == 'SUBUNIT':
                    for txt in comment.get('texts', []):
                        subunits.add(txt.get('value'))

        first_desc = results[0].get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
        first_accession = results[0].get('primaryAccession')
        
        link = f"https://www.uniprot.org/uniprotkb/{first_accession}/entry" if first_accession else f"https://www.uniprot.org/uniprotkb?query=ec:{ec}"

        return {
            'enzyme_name': name,
            'ec_number': ec,
            'uniprot_link': link,  
            'general_info': {
                'protein_name': first_desc,
                'functions': list(functions)[:3], 
                'catalytic_activities': list(set(catalytic_activities))[:3],
                'pathways': list(pathways)[:3],
                'cofactors': list(cofactors)[:3],
                'subunit_structure': list(subunits)[:2]
            },
            'total_entries_analyzed': len(results)
        }