
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

    def fetch_protein_data(self, protein_name: str, identifier: str, id_type: str = "EC", max_entries: int = 30) -> Dict[str, Any]:
        """Searches UniProt for general protein properties matching the EC or TC number."""

        clean_key_name = protein_name.replace(' ', '_').replace('/', '_').replace('-', '_')
        cache_key = f"{id_type}_{identifier}_{clean_key_name}"
        
        if cache_key in self.cache:
            return self.cache[cache_key]

        query_name = protein_name.replace('/', ' ').replace('-', ' ')

        if id_type == "EC":
            query_str = f'ec:{identifier}'
        elif id_type == "TC":
            query_str = f'xref:tcdb-{identifier} AND "{query_name}"'
        else:
            query_str = f'{identifier} AND "{query_name}"'

        def _do_request(q_str):
            """Auxiliary function to avoid repeating the requests block."""

            params = {
                'query': q_str,
                'format': 'json',
                'size': max_entries,
                'fields': 'accession,protein_name,cc_function,cc_catalytic_activity,cc_pathway,cc_cofactor,cc_subunit'
            }
            resp = requests.get(self.BASE_URL, params=params, timeout=10)
            resp.raise_for_status()
            return resp.json().get('results', [])

        try:
            results = []
            
            try:
                results = _do_request(query_str)
            except requests.exceptions.HTTPError as req_err:
                logger.warning(f"Specific query rejected by UniProt for '{protein_name}': {req_err}")

            if not results and id_type == "TC":
                logger.info(f"No exact match (or error) for '{protein_name}'. Falling back to general TC family {identifier} data.")
                fallback_query = f'xref:tcdb-{identifier}'
                results = _do_request(fallback_query)

            if not results:
                self.cache[cache_key] = {}
                self._save_cache()
                return {}
            
            summary = self._summarize_results(protein_name, identifier, id_type, results)
            self.cache[cache_key] = summary
            self._save_cache()
            return summary
            
        except Exception as e:
            logger.error(f"[UniprotService] Error fetching general data for {protein_name} ({id_type} {identifier}): {e}")
            return {}

    def _summarize_results(self, name: str, identifier: str, id_type: str, results: List[dict]) -> Dict[str, Any]:
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
        
        if first_accession:
            link = f"https://www.uniprot.org/uniprotkb/{first_accession}/entry"
        else:
            if id_type == "EC":
                link = f"https://www.uniprot.org/uniprotkb?query=ec:{identifier}"
            else:
                link = f"https://www.uniprot.org/uniprotkb?query=xref:tcdb-{identifier}"

        return {
            'protein_name': name,
            'identifier': identifier,
            'id_type': id_type,
            'uniprot_link': link,  
            'general_info': {
                'protein_name_uniprot': first_desc,
                'functions': list(functions)[:3], 
                'catalytic_activities': list(set(catalytic_activities))[:3],
                'pathways': list(pathways)[:3],
                'cofactors': list(cofactors)[:3],
                'subunit_structure': list(subunits)[:2]
            },
            'total_entries_analyzed': len(results)
        }