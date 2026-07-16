
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
        "maltase":                    EnzymeInfo("3.2.1.20", "maltose"),
        "fructokinase":               EnzymeInfo("2.7.1.4",  "fructose"),
        "galactokinase":              EnzymeInfo("2.7.1.6",  "galactose"),
        "glucokinase":                EnzymeInfo("2.7.1.2",  "glucose"),
        "hexokinase":                EnzymeInfo("2.7.1.1",  "glucose"),
        "beta-amylase":               EnzymeInfo("3.2.1.2",  "starch"),
        "alpha-amylase":              EnzymeInfo("3.2.1.1",  "starch"),
        "cellulase":                  EnzymeInfo("3.2.1.4",  "cellulose")
    }

    TRANSPORTERS: Dict[str, TransporterInfo] = {
        # --- Família MFS ---
        "glucose transporter":        TransporterInfo("2.A.1.1",  "glucose", "MFS"),
        "fructose transporter":       TransporterInfo("2.A.1.1",  "fructose", "MFS"),
        "galactose transporter":      TransporterInfo("2.A.1.1",  "galactose", "MFS"),
        "hexose transporter":         TransporterInfo("2.A.1.1",  "glucose/fructose/galactose", "MFS"),
        "sucrose transporter":        TransporterInfo("2.A.1.1",  "sucrose", "MFS"),
        "maltose transporter":        TransporterInfo("2.A.1.1",  "maltose", "MFS"),
        "lactose permease":           TransporterInfo("2.A.1.14", "lactose", "MFS"),
        "lactose transporter":        TransporterInfo("2.A.1.14", "lactose", "MFS"),
        
        # --- Família SWEET ---
        "sweet sugar transporter":    TransporterInfo("2.A.123",  "glucose/sucrose", "SWEET"),
        "sweet transporter":          TransporterInfo("2.A.123",  "glucose/sucrose", "SWEET"),
        
        # --- Sistemas ABC ---
        "maltose abc transporter":    TransporterInfo("3.A.1.1",  "maltose", "ABC"),
        "galactose abc transporter":  TransporterInfo("3.A.1.1",  "galactose", "ABC"),
        
        # --- Sistemas PTS (Cianobactérias) ---
        "pts system glucose":         TransporterInfo("4.A.1",    "glucose", "PTS"),
        "pts system fructose":        TransporterInfo("4.A.2",    "fructose", "PTS"),
        "pts system sucrose":         TransporterInfo("4.A.1",    "sucrose", "PTS"),
        "pts system lactose":         TransporterInfo("4.A.3",    "lactose", "PTS"),

        "sugar abc transporter":      TransporterInfo("3.A.1.1",  "broad specificity", "ABC"),
        "sugar transporter":          TransporterInfo("2.A.1.1",  "broad specificity", "MFS")
    }

    FUTURE_TERMS: List[str] = [
    "hypothetical", "similar", "putative", 
    "uncharacterized", "probable", "possible",
    "potential", "like", "related"
    ]

    EXCLUDED_TERMS_ENZYMES: List[str] = [
        "inhibitor", "regulator", "activator", "repressor", 
        "receptor", "transcription factor", "fingers", 
        "binding protein", "domain-containing",
        "dna", "rna", "trna", "mrna", "ribosomal", "ribosome",
        "recombinase", "integrase", "transposase", "nuclease",
        "polymerase", "helicase", "chromosome", "plasmid",
        "protein kinase", "histidine kinase", "tyrosine kinase", 
        "serine/threonine", "signal transduction", "partial",
        "assembly", "pdb", "chain", "fragment", "phage", "viral", "virus", "capsid",
        "chaperone", "flagellar", "flagellum", "pilus", "porin",
        "synthetic", "construct", "vector", "fusion", "mutant", "chimeric",
        "operon", "promoter", "sensor", "toxin", "antitoxin", "domain",
        "terminal", "n-terminal", "ferredoxin", "c-terminal"
    ]

    EXCLUDED_TERMS_TRANSPORTERS: List[str] = [
        "inhibitor", "regulator", "activator", "repressor", 
        "transcription factor", "fingers",
        "dna", "rna", "trna", "mrna", "ribosomal", "ribosome",
        "recombinase", "integrase", "transposase", "nuclease",
        "polymerase", "helicase", "chromosome", "plasmid",
        "protein kinase", "histidine kinase", "tyrosine kinase", 
        "serine/threonine", "partial", "synthase", "biosynthesis",
        "assembly", "pdb", "chain", "fragment", "phage", "viral", "virus", "capsid",
        "chaperone", "flagellar", "flagellum", "pilus", 
        "synthetic", "construct", "vector", "fusion", "mutant", "chimeric",
        "operon", "promoter", "toxin", "antitoxin", "terminal", "n-terminal", "c-terminal"
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

    FEEDIPEDIA_RESIDUES: Dict[str, Dict[str, Dict[str, int]]] = {

        "Cereal grains and by-products": {

            "Barley distillers grains (ethanol)": 19499,
            "Brewers grains, dehydrated": 11893,
            "Maize bran": 12280,
            "Hominy feed": 12290,
            "Maize cobs": 12875, 
            "Maize stover, dry": 12874,
            "Corn thin stillage": 12852,
            "Corn gluten feed": 12287,
            "Barley rootlets, dehydrated": 11829,
            "Oat hulls": 12387,
            "Oat mill feed": 15136,
            "Rice bran, defatted, fibre > 20%": 11641, 
            "Rice hulls": 11643,
            "Wheat middlings": 12756 
        },

        "Legume seeds and by-products": {

            "Carob pod meal, without seeds": 11933,
            "Carob germ, dehydrated": 11932,
            "Chickpea bran (chuni)": 11971,
            "Chickpea pod husks": 17367,
            "Chickpea straw": 11976,
            "Common bean straw": 12006,
            "Faba bean aerial part, straw": 19703,
            "Guar meal crop residues": 11703,
            "Guar meal": 11697,
            "Guar meal, high protein": 11698,
            "Guar meal, raw": 11700,
            "Guar meal forage, dry": 11702,
            "Lentil screenings": 12249,
            "Lentil bran": 12246,
            "Lentil pod husks": 12247,
            "Lentil straw": 12250,
            "Lupin seed hulls": 24366,
            "Peanut hulls": 12155,
            "Peanut skins": 12166,
            "Peanut crop residues, dry": 12167
        },

        "Oil plants and by-products": {

            "Almond hulls": 11756,
            "Almond shells": 26577,
            "Cocoa hulls": 11692,
            "Cocoa pod husks": 15454,
            "Cotton straw": 12019,
            "Cottonseed hulls": 12020,
            "Oil palm fronds, dried": 15389,
            "Olive oil cake, exhausted, with stones": 12887,
            "Olive oil cake, exhausted, without stones": 12888,
            "Olive oil vegetation water": 12889,
            "Palm oil mill effluent, decanted": 15497,
            "Palm press fibre, low oil": 15394,
            "Rapeseed hulls": 12501,
            "Soybean hulls": 12623,
            "Soybean straw": 19916,
            "Sunflower stover (stalks and heads)": 12662,
            "Sunflower hulls": 12661,
            "Sunflower screenings": 15601
        },

        "Fruits and by-products": {

            "Apple pomace, dehydrated": 22380,
            "Citrus molasses": 11986,
            "Citrus pulp, dried": 11987,
            "Coffee husks": 11612,
            "Coffee pulp, dehydrated": 11613,
            "Instant coffee byproduct": 11615,
            "Grape pomace, dehydrated": 12133,
            "Guava waste, dried": 20399,
            "Mango juice extraction by-product, dried": 14827,
            "Papaya pomace, dried": 11554,
            "Pineapple canning byproduct, dehydrated": 12453,
            "Pumpkin hulls": 12481,
            "Tomato pomace, dehydrated": 12705,
            "Watermelon hulls": 12743
        },

        "Roots, tubers and by-products": {

            "Beet pulp, dehydrated": 11838,
            "Beet molasses": 12340,
            "Cassava pomace, dehydrated": 11945,
            "Cassava peels, dry": 12805,
            "Potato pulp, dehydrated": 26045,
            "Steamed potato peels (liquid potato feed)": 26046
        },

        "Sugar processing by-products": {

            "Sugarcane molasses": 12341,
            "Sugarcane bagasse, dehydrated": 11673,
            "Sugarcane bagasse pith": 11677,
            "Sugarcane filter-press mud": 12653
        },

        "Dairy by-products": {

            "Whey, sweet, dehydrated, skimmed": 12893,
            "Whey, acid, dehydrated, skimmed": 12894
        }
    }

    @staticmethod
    def get_enzyme_names() -> List[str]:
        return list(AppConfig.ENZYMES.keys())