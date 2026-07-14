
import os
import sys
import json
import time
import random
import requests
from bs4 import BeautifulSoup

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

from config import AppConfig

class FeedipediaService:
    """Service to scrape and extract data from Feedipedia database."""

    def __init__(self):
        self.base_url = "https://www.feedipedia.org/node/{}"
        self.headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36'
        }
        
        self.data_dir = os.path.join(parent_dir, "data")
        self.output_file = os.path.join(self.data_dir, "feedipedia_raw_data.json")
        
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)

    def extract_node_ids(self, data_dict, main_category=None, result_list=None):
        """Recursively extracts node IDs from the nested dictionary structure."""

        if result_list is None:
            result_list = []
            
        for key, value in data_dict.items():
            if isinstance(value, dict):
                if main_category is None:
                    self.extract_node_ids(value, main_category=key, result_list=result_list)
                else:
                    self.extract_node_ids(value, main_category=main_category, result_list=result_list)
            
            elif isinstance(value, int):
                result_list.append({
                    "Main_Category": main_category or "Uncategorized",
                    "Residue_Name": key,
                    "Node_ID": value
                })
                
        return result_list

    def fetch_table_data(self, node_id):
        """Scrapes the 'Main Analysis' and 'Minerals' tables for a given node."""

        url = self.base_url.format(node_id)
        
        try:
            response = requests.get(url, headers=self.headers, timeout=10)
            
            extracted_data = {
                "Main analysis": {},
                "Minerals": {}
            }

            if response.status_code != 200:
                print(f"HTTP Error {response.status_code} for Node {node_id}")
                return None

            soup = BeautifulSoup(response.content, 'html.parser')
            tables = soup.find_all('table')
            
            for table in tables:
                rows = table.find_all('tr')
                table_category = None
                column_names = []
                
                for row in rows:
                    cells = row.find_all(['td', 'th'])
                    cell_texts = [cell.text.strip() for cell in cells]
                    
                    if "Unit" in cell_texts and "Avg" in cell_texts:
                        column_names = cell_texts
                        first_cell = cell_texts[0]
                        
                        if "Main analysis" in first_cell:
                            table_category = "Main analysis"
                        elif "Minerals" in first_cell:
                            table_category = "Minerals"
                        else:
                            table_category = None
                        continue
                        
                    if table_category and column_names and len(cell_texts) > 1:
                        parameter_name = cell_texts[0]
                        parameter_metrics = {}
                        
                        for i in range(1, min(len(column_names), len(cell_texts))):
                            col_name = column_names[i]
                            value = cell_texts[i]
                            
                            if value and value != "*": 
                                parameter_metrics[col_name] = value
                                
                        extracted_data[table_category][parameter_name] = parameter_metrics
                        
            return extracted_data
            
        except Exception as e:
            print(f"Request failed for Node {node_id}: {str(e)}")
            return None

    def run_extraction_pipeline(self):
        """Runs the full extraction process and saves to JSON."""

        print("Starting Feedipedia extraction pipeline...")
        
        targets = self.extract_node_ids(AppConfig.FEEDIPEDIA_RESIDUES)
        total_targets = len(targets)
        print(f"Found {total_targets} target residues in configuration.")
        
        master_dataset = []
        
        for index, target in enumerate(targets, 1):
            name = target["Residue_Name"]
            node_id = target["Node_ID"]
            
            print(f"[{index}/{total_targets}] Fetching: {name} (Node: {node_id})...")
            
            scraped_data = self.fetch_table_data(node_id)
            
            if scraped_data:
                record = {
                    "Main_Category": target["Main_Category"], 
                    "Residue_Name": name,
                    "Node_ID": node_id,
                    "Data": scraped_data
                }
                master_dataset.append(record)
            
            if index < total_targets:
                time.sleep(random.uniform(1.5, 3.5))
                
        with open(self.output_file, 'w', encoding='utf-8') as f:
            json.dump(master_dataset, f, indent=4, ensure_ascii=False)

        relative_path = os.path.relpath(self.output_file)

        print(f"Extraction complete. Data saved to {relative_path}")
            

if __name__ == "__main__":
    service = FeedipediaService()
    service.run_extraction_pipeline()