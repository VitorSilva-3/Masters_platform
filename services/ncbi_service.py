
from Bio import Entrez, Medline
from typing import List, Dict, Any

class NCBIService:
    """Service responsible for interacting with NCBI databases (PubMed, Taxonomy)."""

    def __init__(self, email: str):
        """Initialize with user email for API tracking."""
        self.email = email
        Entrez.email = self.email

    def fetch_taxonomy_lineage(self, organism_name: str) -> str:
        """
        Fetches the full taxonomy lineage for a specific organism.
        Returns a string with the lineage or an error message.
        """
        try:
            # Search for the organism ID
            handle = Entrez.esearch(db="taxonomy", term=f"{organism_name}[Scientific Name]")
            record = Entrez.read(handle)
            handle.close()

            if not record["IdList"]:
                return "Taxonomy ID not found."

            # Fetch details using the ID
            tax_id = record["IdList"][0]
            handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            details = Entrez.read(handle)
            handle.close()
            
            # Extract lineage string
            if details:
                return details[0].get("Lineage", "Lineage not available")
            return "No details found."
            
        except Exception as e:
            return f"Error fetching taxonomy: {str(e)}"

    def search_articles(self, organism: str, enzyme: str, ec_number: str, max_results: int = 5) -> List[Dict[str, Any]]:
        """
        Searches PubMed for articles linking an organism to an enzyme/EC number.
        """
        # Constructing a specific query
        query = f'"{organism}"[Organism] AND ("{enzyme}"[Title/Abstract] OR "{ec_number}"[Title/Abstract])'
        
        try:
            # Search for Article IDs (PMIDs)
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            pmids = record.get("IdList", [])
            if not pmids:
                return []

            # Fetch Article Details
            handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="text")
            articles = []
            
            # Parse Medline format
            for article in Medline.parse(handle):
                articles.append({
                    'title': article.get('TI', 'No title available'),
                    'authors': article.get('AU', []),
                    'journal': article.get('JT', 'Unknown Journal'),
                    'year': article.get('DP', 'Unknown Year'),
                    'pmid': article.get('PMID', 'N/A'),
                    'abstract': article.get('AB', 'No abstract available.')
                })
            handle.close()
            return articles
            
        except Exception as e:
            print(f"Error searching PubMed for {organism} + {enzyme}: {e}")
            return []