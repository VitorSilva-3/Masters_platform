
from Bio.KEGG import REST
from typing import Dict, Any

class KeggService:
    """Service to interact with the KEGG Enzyme Database."""

    def fetch_enzyme_details(self, enzyme_name: str, ec_number: str) -> Dict[str, Any]:
        """
        Fetches detailed information from KEGG for a specific enzyme.
        """
        entry_id = f"ec:{ec_number}"
        try:
            # Fetch data from KEGG API
            raw_data = REST.kegg_get(entry_id).read()
            return self._parse_kegg_text(raw_data, ec_number)
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
                # Clean up the name string
                info['name'] = line.split("NAME")[1].strip().rstrip(";")
            elif line.startswith("DEFINITION"):
                info['definition'] = line.split("DEFINITION")[1].strip()
            elif line.startswith("REACTION"):
                info['reaction'] = line.split("REACTION")[1].strip()
            elif line.startswith("PATHWAY"):
                # Extract pathway code and description
                parts = line.split()
                if len(parts) > 2:
                    code = parts[1]
                    desc = " ".join(parts[2:])
                    info['pathways'].append({'code': code, 'description': desc})

        return info