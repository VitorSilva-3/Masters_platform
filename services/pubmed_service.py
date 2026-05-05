
import logging
from Bio import Entrez, Medline
from typing import List, Dict, Any
from utils import load_json_cache, save_json_cache, setup_ncbi_entrez

logger = logging.getLogger(__name__)

class PubMedService:
    """Service to search NCBI PubMed for articles."""

    def __init__(self, email: str, cache_file: str = "data/literature_cache.json"):
        self.email = email
        setup_ncbi_entrez(self.email)
        self.cache_file = cache_file
        self.cache = load_json_cache(self.cache_file, service_name = "PubMedService")

    def search_articles(self, organism: str, protein_name: str, identifier: str, keywords: List[str] = None, max_results: int = 15) -> List[Dict[str, Any]]:
        """Searches PubMed for articles."""

        cache_key = f"{organism}|{protein_name}|{identifier}"
        
        if cache_key in self.cache:
            return self.cache[cache_key]

        query_name = protein_name.replace('/', ' ').replace('-', ' ')

        if keywords:
            keywords_block = " OR ".join([f'"{kw}"[Title/Abstract]' for kw in keywords])
            query = (
                f'("{organism}"[Organism] OR "{organism}"[Title/Abstract]) AND '
                f'("{query_name}"[Title/Abstract] OR "{identifier}"[Title/Abstract] OR {keywords_block})'
            )
        else:
            query = (
                f'("{organism}"[Organism] OR "{organism}"[Title/Abstract]) AND '
                f'("{query_name}"[Title/Abstract] OR "{identifier}"[Title/Abstract])'
            )
        
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            pmids = record.get("IdList", [])
            if not pmids:
                self.cache[cache_key] = []
                save_json_cache(self.cache_file, self.cache, service_name = "PubMedService")
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
            
            self.cache[cache_key] = articles
            save_json_cache(self.cache_file, self.cache, service_name = "PubMedService")
            
            return articles
            
        except Exception as e:
            logger.error(f"[PubMedService] Error searching PubMed for {organism} + {protein_name}: {e}")
            return []
        
    def save_cache(self):
        """Saves the cache to disk."""
        save_json_cache(self.cache_file, self.cache, service_name = "PubMedService")