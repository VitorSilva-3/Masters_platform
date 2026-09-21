
"""
Biological rules for the content-based recommendation system.
Maps nutrients to the required Enzyme EC numbers and Transporter targets.
"""

# The specific nutrients extracted from Feedipedia analysis
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

# DICTIONARY 1: FOR COMPLEX AGRO-INDUSTRIAL WASTE (FEEDIPEDIA)
METABOLIC_PATHWAYS = {
    "Lactose_Pathway": {
        "feedipedia_triggers": ["Lactose"],
        "enzymes": ["3.2.1.23", "2.7.1.6", "2.7.1.1", "2.7.1.2"], # lactase, galactokinase, hexokinase, glucokinase
        "transporters": ["lactose", "glucose", "galactose", "glucose/fructose/galactose", "glucose/sucrose", "broad specificity"],
        "requires_extracellular": False
    },
    "Starch_Pathway": {
        "feedipedia_triggers": ["Starch (polarimetry)", "Starch (enzymatic)"],
        "enzymes": ["3.2.1.1", "3.2.1.2", "3.2.1.20", "2.7.1.1", "2.7.1.2"], # alpha-amylase, beta-amylase, maltase, hexokinase, glucokinase
        "transporters": ["maltose", "glucose", "glucose/fructose/galactose", "glucose/sucrose", "broad specificity"],
        "requires_extracellular": True 
    },
    "Cellulose_Pathway": {
        "feedipedia_triggers": ["Crude fibre", "NDF", "ADF", "Neutral detergent fibre", "Acid detergent fibre"],
        "enzymes": ["3.2.1.4", "3.2.1.21", "2.7.1.1", "2.7.1.2"], # cellulase, beta-glucosidase, hexokinase, glucokinase
        "transporters": ["glucose", "cellobiose", "glucose/fructose/galactose", "glucose/sucrose", "broad specificity"],
        "requires_extracellular": True
    },
    "Fructose_Pathway": {
        "feedipedia_triggers": ["Fructose"],
        "enzymes": ["2.7.1.4"], # fructokinase
        "transporters": ["fructose", "glucose/fructose/galactose", "broad specificity"],
        "requires_extracellular": False
    },
    "Mixed_Sugars_Pathway": {
        # Catch-all for undefined Feedipedia "Total sugars" (usually sucrose + glucose + fructose)
        "feedipedia_triggers": ["Total sugars"],
        "enzymes": ["3.2.1.26", "2.7.1.1", "2.7.1.2", "2.7.1.4"], # invertase, sucrase, beta-fructofuranosidase, hexokinase, glucokinase, fructokinase
        "transporters": ["sucrose", "fructose", "glucose", "glucose/fructose/galactose", "glucose/sucrose", "broad specificity"],
        "requires_extracellular": False
    }
}

# DICTIONARY 2: FOR PURE SUBSTRATE SELECTION
PURE_SUGARS_PATHWAYS = {
    "glucose": {
        "enzymes": ["2.7.1.1", "2.7.1.2"], # hexokinase, glucokinase
        "transporters": ["glucose", "glucose/fructose/galactose", "glucose/sucrose", "broad specificity"],
        "requires_extracellular": False
    },
    "fructose": {
        "enzymes": ["2.7.1.4"], # fructokinase
        "transporters": ["fructose", "glucose/fructose/galactose", "broad specificity"],
        "requires_extracellular": False
    },
    "galactose": {
        "enzymes": ["2.7.1.6"], # galactokinase
        "transporters": ["galactose", "glucose/fructose/galactose", "broad specificity"],
        "requires_extracellular": False
    },
    "sucrose": {
        "enzymes": ["3.2.1.26", "2.7.1.1", "2.7.1.2", "2.7.1.4"], # invertase, sucrase, beta-fructofuranosidase, hexokinase, glucokinase, fructokinase
        "transporters": ["sucrose", "glucose/sucrose", "broad specificity"],
        "requires_extracellular": False # Can be intra or extra depending on the strain
    },
    "lactose": {
        "enzymes": ["3.2.1.23", "2.7.1.6", "2.7.1.1", "2.7.1.2"], # lactase, beta-galactosidase, galactokinase, hexokinase, glucokinase
        "transporters": ["lactose", "glucose", "galactose", "glucose/sucrose", "glucose/fructose/galactose", "broad specificity"],
        "requires_extracellular": False
    },
    "maltose": {
        "enzymes": ["3.2.1.20", "2.7.1.1", "2.7.1.2"], # maltase, hexokinase, glucokinase
        "transporters": ["maltose", "glucose", "glucose/sucrose", "glucose/fructose/galactose", "broad specificity"],
        "requires_extracellular": False
    },
    "starch": {
        "enzymes": ["3.2.1.1", "3.2.1.2", "3.2.1.20", "2.7.1.1", "2.7.1.2"], # alpha-amylase, beta-amylase, maltase, hexokinase, glucokinase
        "transporters": ["maltose", "glucose", "glucose/sucrose", "glucose/fructose/galactose", "broad specificity"],
        "requires_extracellular": True
    },
    "cellulose": {
        "enzymes": ["3.2.1.4", "3.2.1.21", "2.7.1.1", "2.7.1.2"], # cellulase, beta-glucosidase, hexokinase, glucokinase
        "transporters": ["glucose", "cellobiose", "glucose/sucrose", "glucose/fructose/galactose", "broad specificity"],
        "requires_extracellular": True
    }
}