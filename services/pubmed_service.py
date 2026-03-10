
import json
import os
from Bio import Entrez, Medline
from typing import List, Dict, Any

class PubMedService:
    """Service to search NCBI PubMed for enzyme-related articles with local caching."""

    def __init__(self, email: str, cache_file: str = "data/pubmed_cache.json"):
        """Initialize with user email and load local PubMed cache."""
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
                print(f"[PubMedService] Error loading cache: {e}")
        return {}

    def _save_cache(self):
        try:
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump(self.cache, f, indent=4)
        except Exception as e:
            print(f"[PubMedService] Error saving cache: {e}")

    def search_articles(self, organism: str, enzyme: str, ec_number: str, max_results: int = 5) -> List[Dict[str, Any]]:
        """
        Searches PubMed for articles linking an organism to an enzyme.
        """
          
        cache_key = f"{organism}|{enzyme}|{ec_number}"
        
        if cache_key in self.cache:
            return self.cache[cache_key]

        query = f'"{organism}"[Organism] AND ("{enzyme}"[Title/Abstract] OR "{ec_number}"[Title/Abstract])'
        
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            pmids = record.get("IdList", [])
            if not pmids:
                self.cache[cache_key] = []
                self._save_cache()
                return []

            handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="text")
            articles = []
            
            for article in Medline.parse(handle):
                articles.append({
                    'title': article.get('TI', 'No title available'),
                    'authors': article.get('AU', []),
                    'journal': article.get('JT', 'Unknown Journal'),
                    'year': article.get('DP', 'Unknown Year'),
                    'pmid': article.get('PMID', 'N/A'),
                    'abstract': article.get('AB', 'No abstract available.'),
                    'doi': article.get('LID', '').split(' [doi]')[0] if 'LID' in article else None
                })
            handle.close()
            
            self.cache[cache_key] = articles
            self._save_cache()
            
            return articles
            
        except Exception as e:
            print(f"Error searching PubMed for {organism} + {enzyme}: {e}")
            return []