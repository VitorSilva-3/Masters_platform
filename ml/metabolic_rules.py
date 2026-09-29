"""
Biological rules for the content-based recommendation system.
Maps nutrients to the required Enzyme EC numbers and Transporter targets.
"""

NUTRIENT_TARGETS = [
    "Lactose",
    "Starch (polarimetry)",
    "Starch (enzymatic)",
    "Crude fibre",
    "NDF",
    "Neutral detergent fibre",
    "ADF",
    "Acid detergent fibre",
    "Total sugars",
    "Fructose"
]

METABOLIC_PATHWAYS = {
    "Lactose_Pathway": {
        "feedipedia_triggers": ["Lactose"],
        "enzymes": ["3.2.1.23", "2.7.1.6", "2.7.1.1", "2.7.1.2"], 
        "transporters": ["lactose", "galactose", "glucose"],
        "requires_extracellular": False
    },
    "Starch_Pathway": {
        "feedipedia_triggers": ["Starch (polarimetry)", "Starch (enzymatic)"],
        "enzymes": ["3.2.1.1", "3.2.1.2", "3.2.1.20", "2.7.1.1", "2.7.1.2"], 
        "transporters": ["maltose", "glucose"],
        "requires_extracellular": True 
    },
    "Cellulose_Pathway": {
        "feedipedia_triggers": ["Crude fibre", "NDF", "ADF", "Neutral detergent fibre", "Acid detergent fibre"],
        "enzymes": ["3.2.1.4", "3.2.1.21", "2.7.1.1", "2.7.1.2"], 
        "transporters": ["cellobiose", "glucose"],
        "requires_extracellular": True
    },
    "Fructose_Pathway": {
        "feedipedia_triggers": ["Fructose"],
        "enzymes": ["2.7.1.4"], 
        "transporters": ["fructose"],
        "requires_extracellular": False
    },
    "Mixed_Sugars_Pathway": {
        "feedipedia_triggers": ["Total sugars"],
        "enzymes": ["3.2.1.26", "2.7.1.1", "2.7.1.2", "2.7.1.4"], 
        "transporters": ["sucrose", "fructose", "glucose"],
        "requires_extracellular": False
    }
}

PURE_SUGARS_PATHWAYS = {
    "glucose": {
        "enzymes": ["2.7.1.1", "2.7.1.2"], 
        "transporters": ["glucose"],
        "requires_extracellular": False
    },
    "fructose": {
        "enzymes": ["2.7.1.4"], 
        "transporters": ["fructose"],
        "requires_extracellular": False
    },
    "galactose": {
        "enzymes": ["2.7.1.6"], 
        "transporters": ["galactose"],
        "requires_extracellular": False
    },
    "sucrose": {
        "enzymes": ["3.2.1.26", "2.7.1.1", "2.7.1.2", "2.7.1.4"], 
        "transporters": ["sucrose", "fructose", "glucose"],
        "requires_extracellular": False 
    },
    "lactose": {
        "enzymes": ["3.2.1.23", "2.7.1.6", "2.7.1.1", "2.7.1.2"], 
        "transporters": ["lactose", "galactose", "glucose"],
        "requires_extracellular": False
    },
    "maltose": {
        "enzymes": ["3.2.1.20", "2.7.1.1", "2.7.1.2"], 
        "transporters": ["maltose", "glucose"],
        "requires_extracellular": False
    },
    "starch": {
        "enzymes": ["3.2.1.1", "3.2.1.2", "3.2.1.20", "2.7.1.1", "2.7.1.2"], 
        "transporters": ["maltose", "glucose"],
        "requires_extracellular": True
    },
    "cellulose": {
        "enzymes": ["3.2.1.4", "3.2.1.21", "2.7.1.1", "2.7.1.2"], 
        "transporters": ["cellobiose", "glucose"],
        "requires_extracellular": True
    }
}