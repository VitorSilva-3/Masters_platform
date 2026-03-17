

# -----------------------------
# Currently waiting for the Semantic Scholar response regarding API access
# The code is working with the email-based User-Agent header
# -----------------------------


import requests
import time
import logging
from typing import List, Dict, Any

logger = logging.getLogger(__name__)

class SemanticScholarService:
    """Service to search Semantic Scholar for enzyme-related articles."""

    def __init__(self, email: str):
        """Initialize with user email for API requests."""

        self.headers = {
            "User-Agent": f"MicroValue_Thesis_Project/1.0 (mailto:{email})"
        }
        self.base_url = "https://api.semanticscholar.org/graph/v1/paper/search"

    def search_articles(self, species: str, enzyme: str, keywords: List[str] = None, max_results: int = 15) -> List[Dict[str, Any]]:
        """Searches for articles in Semantic Scholar."""

        query_parts = [f'"{species}"', f'"{enzyme}"']
        
        if keywords:
            kw_string = " ".join(keywords)
            query_parts.append(kw_string)
            
        query = " ".join(query_parts).strip()
        
        params = {
            "query": query,
            "limit": max_results * 2, 
            "fields": "title,abstract,year,authors,venue,url,externalIds"
        }

        try:
            max_retries = 3
            for attempt in range(max_retries):
                response = requests.get(self.base_url, params=params, headers=self.headers)
                
                if response.status_code == 429:
                    sleep_time = 10 * (attempt + 1) 
                    logger.warning(f"[SemanticScholarService] Rate limit hit! Sleeping for {sleep_time}s (Attempt {attempt+1}/{max_retries})...")
                    time.sleep(sleep_time)
                    continue 
                
                response.raise_for_status()
                break
            else:
                logger.error(f"[SemanticScholarService] Max retries reached for {species}. Skipping S2.")
                return []

            data = response.json()
            articles = []
            
            if "data" not in data:
                return articles

            for item in data["data"]:
                title = item.get("title", "")
                abstract = item.get("abstract", "")

                if not abstract:
                    continue

                full_text = f"{title} {abstract}".lower()
                
                if species.lower() not in full_text:
                    continue 

                authors_list = item.get("authors", [])
                if authors_list:
                    if len(authors_list) > 2:
                        authors_str = f"{authors_list[0]['name']} et al."
                    else:
                        authors_str = ", ".join([a["name"] for a in authors_list])
                else:
                    authors_str = "Unknown Authors"

                external_ids = item.get("externalIds", {})
                doi = external_ids.get("DOI", "")
                pmid = external_ids.get("PubMed", "")
                
                article_url = item.get("url", "")
                if doi:
                    article_url = f"https://doi.org/{doi}"
                elif pmid:
                    article_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

                articles.append({
                    "title": title,
                    "authors": authors_str,
                    "journal": item.get("venue", "Unknown Journal"),
                    "year": str(item.get("year", "")) if item.get("year") else "N/A",
                    "abstract": abstract,
                    "url": article_url,
                    "pmid": pmid,
                    "doi": doi,
                    "source": "Semantic Scholar"
                })

                if len(articles) >= max_results:
                    break

            return articles

        except Exception as e:
            logger.error(f"[SemanticScholarService] Error searching for {species} + {enzyme}: {e}")
            return []