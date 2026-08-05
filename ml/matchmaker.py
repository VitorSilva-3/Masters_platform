
import pandas as pd
import numpy as np
from ml.metabolic_rules import METABOLIC_PATHWAYS, PURE_SUGARS_PATHWAYS

class BiotransformationMatchmaker:
    """
    A class to recommend microalgae and cyanobacteria strains based on the metabolic potential
    to biotransform specific agro-industrial residues or pure sugars.
    """

    def __init__(self, df_feedipedia: pd.DataFrame, df_strains_enzymes: pd.DataFrame, df_strains_transporters: pd.DataFrame):
        self.df_feedipedia = df_feedipedia
        self.df_enzymes = df_strains_enzymes
        self.df_transporters = df_strains_transporters
        self.pathways = METABOLIC_PATHWAYS
        self.pure_sugars_pathways = PURE_SUGARS_PATHWAYS
        
        self.enz_strain_col = self._find_column(self.df_enzymes, ['strain', 'organism', 'name'])
        self.trans_strain_col = self._find_column(self.df_transporters, ['strain', 'organism', 'name'])

    def _find_column(self, df: pd.DataFrame, possible_names: list) -> str:
        if df.empty: return None
        lower_cols = {str(col).lower().strip(): col for col in df.columns}
        for name in possible_names:
            if name.lower() in lower_cols:
                return lower_cols[name.lower()]
        return df.columns[0] if len(df.columns) > 0 else None

    def _vectorize_residue(self, residue_row: pd.Series) -> np.ndarray:
        vector = []
        for pathway_name, rules in self.pathways.items():
            pathway_value = 0.0
            for nutrient in rules["feedipedia_triggers"]:
                for col in residue_row.index:
                    if str(col).strip().lower() == nutrient.lower():
                        val = residue_row[col]
                        if pd.notna(val) and float(val) > 0:
                            pathway_value += float(val)
            
            vector.append(min(pathway_value, 1.0))
            
        return np.array(vector).reshape(1, -1)

    def _calculate_pathway_score(self, enzyme_fraction: float, has_transporter: bool, is_simple_sugar: bool) -> float:
        if is_simple_sugar:
            weight_enz, weight_trans = 0.4, 0.6
        else:
            weight_enz, weight_trans = 0.7, 0.3
            
        score = 0.0
        score += (weight_enz * enzyme_fraction)
        if has_transporter: 
            score += weight_trans
            
        return score

    def _vectorize_strain(self, strain_id: str, pathways_dict: dict, active_pathways: list) -> tuple:
        vector = []
        strain_enz = self.df_enzymes[self.df_enzymes[self.enz_strain_col] == strain_id] if self.enz_strain_col else pd.DataFrame()
        strain_trans = self.df_transporters[self.df_transporters[self.trans_strain_col] == strain_id] if self.trans_strain_col else pd.DataFrame()

        ec_col = self._find_column(strain_enz, ['ec_number', 'ec number', 'ec'])
        loc_col = self._find_column(strain_enz, ['localization','localizations', 'subcellular', 'location', 'deeploc'])
        enz_name_col = self._find_column(strain_enz, ['enzyme', 'protein name', 'description'])

        target_col = self._find_column(strain_trans, ['target sugar', 'target', 'substrate', 'nutrient'])
        trans_name_col = self._find_column(strain_trans, ['transporter', 'protein name', 'description'])

        strain_enz_list = []
        if not strain_enz.empty and ec_col:
            for _, row in strain_enz.iterrows():
                strain_enz_list.append({
                    "ec": str(row[ec_col]).lower().strip(),
                    "loc": str(row[loc_col]).lower() if loc_col else "",
                    "name": str(row[enz_name_col]) if enz_name_col else str(row[ec_col])
                })

        strain_trans_list = []
        if not strain_trans.empty and target_col:
            for _, row in strain_trans.iterrows():
                strain_trans_list.append({
                    "target": str(row[target_col]).lower().strip(),
                    "name": str(row[trans_name_col]) if trans_name_col else str(row[target_col])
                })

        matched_enzymes = set()
        matched_transporters = set()

        for pathway_name, rules in pathways_dict.items():
            enzymes_required = rules["enzymes"]
            enzymes_found = 0
            
            current_pathway_enzymes = []
            for required_ec in enzymes_required:
                enzyme_matched = False
                for enz in strain_enz_list:
                    if str(required_ec).lower().strip() in enz["ec"]:
                        if str(required_ec).startswith("3.") and rules.get("requires_extracellular", False):
                            if "extracellular" in enz["loc"] or "secreted" in enz["loc"]:
                                enzyme_matched = True
                                current_pathway_enzymes.append(f"{enz['name']} (Ext)")
                                break
                        else:
                            enzyme_matched = True
                            current_pathway_enzymes.append(enz['name'])
                            break
                if enzyme_matched:
                    enzymes_found += 1
            
            enzyme_fraction = enzymes_found / len(enzymes_required) if len(enzymes_required) > 0 else 0.0

            current_pathway_transporters = []
            has_trans = False
            for rule_t in rules["transporters"]:
                rule_t_clean = str(rule_t).lower().strip()
                for t_dict in strain_trans_list:
                    if rule_t_clean in t_dict["target"]:
                        has_trans = True
                        current_pathway_transporters.append(t_dict["name"])
                        break 
            
            if pathway_name in active_pathways:
                matched_enzymes.update(current_pathway_enzymes)
                matched_transporters.update(current_pathway_transporters)

            is_simple_sugar = ("Simple_Sugars" in pathway_name or "Fructose" in pathway_name or pathway_name in ["glucose", "fructose", "galactose", "sucrose"])
            vector.append(self._calculate_pathway_score(enzyme_fraction, has_trans, is_simple_sugar))

        return np.array(vector).reshape(1, -1), list(matched_enzymes), list(matched_transporters)

    def recommend_strains(self, residue_name: str) -> pd.DataFrame:
        if residue_name not in self.df_feedipedia.index:
            raise ValueError(f"Residue '{residue_name}' not found.")

        residue_row = self.df_feedipedia.loc[residue_name]
        residue_vector = self._vectorize_residue(residue_row)

        if np.sum(residue_vector) == 0:
            return pd.DataFrame(columns=["Strain", "Metabolic affinity index", "Mapped enzymes", "Mapped transporters"])

        pathway_names = list(self.pathways.keys())
        active_pathways = [pathway_names[i] for i, val in enumerate(residue_vector[0]) if val > 0]

        results = []
        unique_strains = self.df_enzymes[self.enz_strain_col].unique() if self.enz_strain_col else []

        for strain in unique_strains:
            strain_vector, found_ecs, found_trans = self._vectorize_strain(strain, self.pathways, active_pathways)
            
            # Novo Cálculo: Multiplicamos o vetor da estirpe pelo vetor do resíduo.
            # Se o resíduo tem muita celulose (ex: 0.8), vai valorizar muito o score da estirpe para a celulose.
            # Se a estirpe não tem a via completa, o score (magnitude) é penalizado, não importa o ângulo.
            weighted_scores = np.multiply(residue_vector[0], strain_vector[0])
            
            # Normalizamos o score final para uma escala de 0 a 10
            # A soma máxima possível de weighted_scores ocorre se a estirpe tiver nota 1.0 em todas 
            # as vias que o resíduo exige.
            max_possible_score = np.sum(residue_vector[0])
            if max_possible_score > 0:
                final_score = (np.sum(weighted_scores) / max_possible_score) * 10
            else:
                final_score = 0.0

            results.append({
                "Strain": strain,
                "Metabolic affinity index": round(final_score, 2),
                "Mapped enzymes": ", ".join(found_ecs) if found_ecs else "-",
                "Mapped transporters": ", ".join(found_trans) if found_trans else "-"
            })

        if not results:
            return pd.DataFrame(columns=["Strain", "Metabolic affinity index", "Mapped enzymes", "Mapped transporters"])

        df_results = pd.DataFrame(results).sort_values(by="Metabolic affinity index", ascending=False)
        return df_results
        
    def recommend_strains_for_pure_sugar(self, sugar_name: str) -> pd.DataFrame:
        sugar_name_lower = sugar_name.lower()
        if sugar_name_lower not in self.pure_sugars_pathways:
            raise ValueError(f"Sugar '{sugar_name}' is not mapped to any pure metabolic pathway.")

        active_pathways = [sugar_name_lower]
        
        results = []
        unique_strains = self.df_enzymes[self.enz_strain_col].unique() if self.enz_strain_col else []

        for strain in unique_strains:
            strain_vector, found_ecs, found_trans = self._vectorize_strain(strain, self.pure_sugars_pathways, active_pathways)
            
            # Para açúcares puros, o vetor do açúcar é 1.0. 
            # Logo, o Índice é diretamente o score da estirpe para essa via * 10.
            # Aqui localizamos o índice da via no vetor para extrair o score.
            pathway_idx = list(self.pure_sugars_pathways.keys()).index(sugar_name_lower)
            strain_score = strain_vector[0][pathway_idx]
            
            final_score = strain_score * 10

            results.append({
                "Strain": strain,
                "Metabolic affinity index": round(final_score, 2),
                "Mapped enzymes": ", ".join(found_ecs) if found_ecs else "-",
                "Mapped transporters": ", ".join(found_trans) if found_trans else "-"
            })

        if not results:
            return pd.DataFrame(columns=["Strain", "Metabolic affinity index", "Mapped enzymes", "Mapped transporters"])

        df_results = pd.DataFrame(results).sort_values(by="Metabolic affinity index", ascending=False)
        return df_results