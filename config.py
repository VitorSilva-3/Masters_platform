
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

@dataclass
class TransporterInfo:
    """Data structure to hold transporter details."""

    tc_number: str  
    target_sugar: str
    family: str     

class AppConfig:
    """Central configuration for the application."""
    
    EMAIL: str = "vtsilva3@gmail.com"
    
    TARGET_TAXA: List[str] = [
        # --- Green Microalgae ---
        "Chloropicophyceae",
        "Mamiellophyceae",
        "Nephroselmidophyceae",
        "Picocystophyceae",
        "Pseudoscourfieldiophyceae",
        "Pyramimonadophyceae",
        "Chlorodendrophyceae",
        "Chlorophyceae",
        "Pedinophyceae",       
        "Trebouxiophyceae", 

        # --- Red Microalgae ---
        "Rhodellophyceae",       
        "Stylonematophyceae",    
    
        # --- Yellow-green Microalgae ---
        "Xanthophyceae",

        # --- Golden Microalgae ---
        "Chrysophyceae",

        # --- Diatoms ---
        "Bacillariophyceae",
        "Coscinodiscophyceae",
        "Fragilariophyceae",
        "Mediophyceae",     

        # --- Haptophytes ---
        "Prymnesiophyceae",
        "Pavlovophyceae",
        "Rappephyceae",

        # --- Dinoflagellates ---
        "Dinophyceae",

        # --- Ochrophytes ---
        "Eustigmatophyceae",
        "Dictyochophyceae",
        "Raphidophyceae",
        
        # --- Flagellates ---
        "Euglenida",
        "Cryptophyceae",
        
        # --- Ancestral & Specialized Lineages ---
        "Glaucocystophyceae",
        "Chlorarachniophyceae",
                
        # --- Cyanobacteria ---
        "Cyanophyceae"          
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

    TRANSPORTERS: Dict[str, TransporterInfo] = {
        # --- Transportadores de Lactose / Galactose ---
        "lactose permease":           TransporterInfo("2.A.1.14", "lactose", "MFS"),
        "lac y":                      TransporterInfo("2.A.1.14", "lactose", "MFS"),
        "galactose transporter":      TransporterInfo("2.A.1.1",  "galactose", "MFS"),
        "gal2":                       TransporterInfo("2.A.1.1",  "galactose", "MFS"),
        
        # --- Transportadores de Sacarose ---
        "sucrose transporter":        TransporterInfo("2.A.1.1",  "sucrose", "MFS"),
        "sucrose permease":           TransporterInfo("2.A.1.1",  "sucrose", "MFS"),
        "sut1":                       TransporterInfo("2.A.1.1",  "sucrose", "MFS"),
        
        # --- Transportadores de Hexoses ---
        "hexose transporter":         TransporterInfo("2.A.1.1",  "glucose/fructose", "MFS"),
        "hxt":                        TransporterInfo("2.A.1.1",  "glucose/fructose", "MFS"),
        "fructose transporter":       TransporterInfo("2.A.1.1",  "fructose", "MFS"),
        "glucose transporter":        TransporterInfo("2.A.1.1",  "glucose", "MFS"),
        
        # --- Transportadores de Maltose ---
        "maltose transporter":        TransporterInfo("2.A.1.1",  "maltose", "MFS"),
        "maltose permease":           TransporterInfo("2.A.1.1",  "maltose", "MFS"),
        
        # --- Sistemas ABC (Mais complexos, usam ATP) ---
        "sugar abc transporter":      TransporterInfo("3.A.1.1",  "multiple sugars", "ABC"),
        "multiple sugar transport":   TransporterInfo("3.A.1.1",  "multiple sugars", "ABC")
    }

    FUTURE_TERMS: List[str] = [
    "hypothetical", "similar", "putative", 
    "uncharacterized", "probable", "possible",
    "potential", "like", "related"
    ]

    EXCLUDED_TERMS: List[str] = [
    "inhibitor", "regulator", "activator", "repressor", 
    "receptor", "transcription factor", "fingers", 
    "binding protein", "domain-containing",
    "dna", "rna", "trna", "mrna", "ribosomal", "ribosome",
    "recombinase", "integrase", "transposase", "nuclease",
    "polymerase", "helicase", "chromosome", "plasmid",
    "protein kinase", "histidine kinase", "tyrosine kinase", 
    "serine/threonine", "signal transduction", "partial",
    "synthase", "biosynthesis", "assembly", "pdb", "chain",
    "fragment", "phage", "viral", "virus", "capsid",
    "chaperone", "flagellar", "flagellum", "pilus", "porin",
    "synthetic", "construct", "vector", "fusion", "mutant", "chimeric",
    "operon", "promoter", "sensor", "toxin", "antitoxin", "domain",
    "terminal", "n-terminal", "ferredoxin", "n-terminal", "c-terminal"
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