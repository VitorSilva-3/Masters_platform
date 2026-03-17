
import json
import os
import logging
from Bio import Entrez, Medline
from typing import List, Dict, Any

logger = logging.getLogger(__name__)

class PubMedService:
    """Service to search NCBI PubMed for articles."""

    def __init__(self, email: str, cache_file: str = "data/literature_cache.json"):
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
                logger.error(f"[PubMedService] Error loading cache: {e}")
        return {}

    def save_cache(self):
        """Saves the current cache dictionary to a JSON file."""
        try:
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump(self.cache, f, indent=4)
        except Exception as e:
            logger.error(f"[PubMedService] Error saving cache: {e}")

    def search_articles(self, organism: str, enzyme: str, ec_number: str, keywords: List[str] = None, max_results: int = 15) -> List[Dict[str, Any]]:
        """Searches PubMed for articles."""

        cache_key = f"{organism}|{enzyme}|{ec_number}"
        
        if cache_key in self.cache:
            return self.cache[cache_key]

        if keywords:
            keywords_block = " OR ".join([f'"{kw}"[Title/Abstract]' for kw in keywords])
            query = (
                f'("{organism}"[Organism] OR "{organism}"[Title/Abstract]) AND '
                f'("{enzyme}"[Title/Abstract] OR "{ec_number}"[Title/Abstract] OR {keywords_block})'
            )
        else:
            query = (
                f'("{organism}"[Organism] OR "{organism}"[Title/Abstract]) AND '
                f'("{enzyme}"[Title/Abstract] OR "{ec_number}"[Title/Abstract])'
            )
        
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            pmids = record.get("IdList", [])
            if not pmids:
                return [] 

            handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="text")
            articles = []
            
            for article in Medline.parse(handle):
                authors_list = article.get('AU', [])
                authors_str = ", ".join(authors_list) if authors_list else "Unknown Authors"

                date_pub = article.get('DP', 'Unknown Year')
                year = date_pub.split(' ')[0] if date_pub != 'Unknown Year' else 'Unknown Year'

                doi = None
                if 'LID' in article:
                    for lid_part in article['LID'].split(' '):
                        if '10.' in lid_part:  
                            doi = lid_part
                            break

                articles.append({
                    'title': article.get('TI', 'No title available'),
                    'authors': authors_str,
                    'journal': article.get('JT', 'Unknown Journal'),
                    'year': year,
                    'pmid': article.get('PMID', 'N/A'),
                    'abstract': article.get('AB', 'No abstract available.'),
                    'doi': doi,
                })
            handle.close()
            
            return articles
            
        except Exception as e:
            logger.error(f"[PubMedService] Error searching PubMed for {organism} + {enzyme}: {e}")
            return []