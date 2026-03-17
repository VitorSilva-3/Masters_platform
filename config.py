
import re
from dataclasses import dataclass
from typing import Dict, List 

def normalize_text(text: str) -> str:
    if not text: return ""
    return re.sub(r"[^a-z0-9]", "", str(text).lower())

@dataclass
class EnzymeInfo:
    """Data structure to hold enzyme details."""
    
    ec_number: str
    target_sugar: str

class AppConfig:
    """Central configuration for the application."""
    
    EMAIL: str = "vtsilva3@gmail.com"
    
    TARGET_TAXA: List[str] = [
        # --- Green Microalgae (Chlorophyta) ---
        "Chlorophyceae",          # Chlamydomonas, Scenedesmus, Dunaliella
        "Trebouxiophyceae",       # Chlorella, Auxenochlorella
        "Prasinophyceae",         # Ostreococcus (Picoplankton)
        "Mamiellophyceae",        # Micromonas
        
        # --- Diatoms & Ochrophytes ---
        "Bacillariophyceae",      # Diatoms (Phaeodactylum)
        "Eustigmatophyceae",      # Nannochloropsis
        "Chrysophyceae",          # Golden algae
        "Xanthophyceae",          # Yellow-green algae
        
        # --- Red Microalgae ---
        "Porphyridiophyceae",     # Porphyridium (Exopolysaccharides)
        "Cyanidiophyceae",        # Extremophiles (Galdieria - loves waste/acid)
        "Rhodellophyceae",        # Unicellular red algae
        
        # --- Haptophytes ---
        "Prymnesiophyceae",       # Isochrysis, Pavlova
        
        # --- Key Heterotrophs (Crucial for Waste) ---
        "Euglenida",              # Euglena (Excellent for waste sugars)
        "Dinophyceae",            # Crypthecodinium (DHA production)
        "Cryptophyceae",          # Cryptomonas
        
        # --- Blue-Green Algae ---
        "Cyanobacteria"           # Spirulina, Synechocystis, Anabaena
    ]

    ENZYMES: Dict[str, EnzymeInfo] = {
        "beta-galactosidase":         EnzymeInfo("3.2.1.23", "lactose"),
        "lactase":                    EnzymeInfo("3.2.1.23", "lactose"),
        "beta-fructofuranosidase":    EnzymeInfo("3.2.1.26", "sucrose"),
        "invertase":                  EnzymeInfo("3.2.1.26", "sucrose"),
        "sucrase":                    EnzymeInfo("3.2.1.26", "sucrose"),
        "fructokinase":               EnzymeInfo("2.7.1.4",  "fructose"),
        "maltase":                    EnzymeInfo("3.2.1.20", "maltose"),
        "beta-amylase":               EnzymeInfo("3.2.1.2",  "starch"),
        "alpha-amylase":              EnzymeInfo("3.2.1.1",  "starch"),
        "cellulase":                  EnzymeInfo("3.2.1.4",  "cellulose")
    }

    FUTURE_TERMS: List[str] = [
    "hypothetical", "similar", "putative", 
    "uncharacterized", "probable", "partial"
    ]

    EXCLUDED_TERMS: List[str] = [
    "inhibitor", "regulator", "activator", "repressor", 
    "receptor", "transcription factor", "fingers", 
    "binding protein", "domain-containing",
    "dna", "rna", "trna", "mrna", "ribosomal", "ribosome",
    "recombinase", "integrase", "transposase", "nuclease",
    "polymerase", "helicase", "chromosome", "plasmid",
    "protein kinase", "histidine kinase", "tyrosine kinase", 
    "serine/threonine", "signal transduction",
    "synthase", "biosynthesis", "assembly"
    ]

    PUBMED_KEYWORDS: List[str] = [
        # Concepts
        "circular economy", "biorefinery", "valorization", "valorisation", "upcycling", "bioeconomy",
        # General waste
        "byproduct*", "by-product*", "waste*", "agro-industrial", "agroindustrial", "effluent*", "food waste",
        # Dairy
        "dairy", "whey", "cheese whey", "milk permeate",
        # Brewing and wine
        "brewing", "brewery", "spent grain", "BSG", "vinasse", "grape pomace",
        # Cereals
        "cereal*", "bran", "straw", "stover", "husks", "hulls", "lignocellulosic",
        # Fruits and veggies
        "pomace", "peel*", "sugar beet pulp", "molasses"
    ]

    @staticmethod
    def get_enzyme_names() -> List[str]:
        return list(AppConfig.ENZYMES.keys())