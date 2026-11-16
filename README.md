# MicroValue Platform 

**App:** [https://microvalue.streamlit.app/](https://microvalue.streamlit.app/)

## About
This platform was developed as part of a Master's thesis in Bioinformatics. It is an interactive web application designed to integrate and analyze data on agro-industrial byproducts and the biochemical repertoire of various microalgae and cyanobacteria. 

The goal is to identify which byproducts can serve as optimal substrates for specific microalgae species based on their enzymatic and transporter profiles.

## Modules
The platform includes the following analytical modules: 
* app: main page with an overview of the platform’s objectives and functionalities. 
* Data: datasets with the microalgal and cyanobacterial species, their targeted sugars, and the corresponding degrading enzymes or transport proteins.  
* Scientific Literature: collection of PubMed publications related to the microalgae and cyanobacteria. 
* Taxonomy: complete taxonomic classification for each organism. 
* KEGG: functional annotations, biochemical reactions, and metabolic pathway maps linked to the enzymes. 
* UniProt: extensive protein metadata, detailing specific catalytic activities, pathways and subunit structures for the enzymes. 
* Subcellular Localization: spatial distribution predictions of enzymes and transport proteins calculated by DeepLoc2.1 and DeepLocPro. 
* Agro-industrial Waste Composition: biochemical profiles of various agro-industrial byproducts sourced from the Feedipedia database. 
* Matchmaker: recommendation engine that calculates a metabolic affinity index to optimally pair specific waste streams with the most suitable species.
* Chatbot: an integrated LLM assistant providing contextual support, data interpretation and guidance on the platform. 

## Built with
* **Frontend:** Streamlit (Python)
* **Data Handling:** Pandas, Biopython, NumPy
* **Database:** MongoDB
* **AI:** Google Gemini API

Developed by Vítor Silva as part of a Master's degree in Bioinformatics at Universidade do Minho.  
