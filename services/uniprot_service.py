
import requests
from typing import Dict, Any, List

class UniprotService:
    """Service to interact with the UniProtKB Database via API."""

    BASE_URL = "https://rest.uniprot.org/uniprotkb/search"

    def fetch_enzyme_data(self, enzyme_name: str, ec_number: str, max_entries: int = 30) -> Dict[str, Any]:
        """
        Searches UniProt for proteins matching the EC number.
        Returns aggregated stats and general info.
        """
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
                return {}
            
            # Process the raw results to get a summary
            return self._summarize_results(enzyme_name, ec_number, results)
            
        except Exception as e:
            print(f"[UniprotService] Error for {enzyme_name}: {e}")
            return {}

    def _summarize_results(self, name: str, ec: str, results: List[dict]) -> Dict[str, Any]:
        """Aggregates data from multiple entries to give a general enzyme profile."""
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

        # Get protein name from the first result as a representative
        first_desc = results[0].get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')

        return {
            'enzyme_name': name,
            'ec_number': ec,
            'general_info': {
                'protein_name': first_desc,
                'functions': list(functions)[:3], # Limit to top 3
                'catalytic_activities': list(set(catalytic_activities))[:3],
                'pathways': list(pathways)[:3],
                'cofactors': list(cofactors)[:3],
                'subunit_structure': list(subunits)[:2]
            },
            'total_entries_analyzed': len(results)
        }