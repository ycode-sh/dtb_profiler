import logging
from dataclasses import dataclass
from typing import Dict, List, Optional
import pandas as pd
from pathlib import Path
from pysam import VariantRecord

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class DrugData:
    """Container for drug resistance information."""
    drug_name: str
    resistance_gene: str
    resistance_mutation: str
    var_ann: str
    mutation_effect: str
    confidence_grading: str  # Changed to lowercase
    compared_to_version_1: str
    genomic_position: str = ""
    additional_grading_criteria_applied: str = ""  # Changed to lowercase
    comment: str = ""  # Changed to lowercase

class DrugResistanceEnricher:
    """Enriches records with drug resistance information."""
    
    def __init__(self, drug_data_dir: Path, logger=logger):  # Allow logger injection
        self.logger = logger or logging.getLogger(__name__)  # Store logger as instance variable
        self.var_ann_index: Dict[str, List[DrugData]] = {}
        self._load_drug_data(drug_data_dir)

    def _load_drug_data(self, drug_data_dir: Path) -> None:
        """Load and index all drug data files."""
        required_columns = {
            'drug_name', 'resistance_gene', 'resistance_mutation',
            'var_ann', 'mutation_effect', 'confidence_grading',
            'compared_to_version_1'
        }

        for tsv_file in drug_data_dir.glob('*.tsv'):
            df = pd.read_csv(tsv_file, sep='\t')
            
            
            # Validate required columns
            for col in required_columns:
                if col not in df.columns:
                    raise ValueError(f"Missing required column {col} in {tsv_file}")
                if df[col].isnull().any():
                    raise ValueError(f"Found empty values in required column {col} in {tsv_file}")

            # Index by var_ann
            for _, row in df.iterrows():
                drug_data = DrugData(**row.to_dict())
                if drug_data.var_ann not in self.var_ann_index:
                    self.var_ann_index[drug_data.var_ann] = []
                self.var_ann_index[drug_data.var_ann].append(drug_data)

    def enrich_records(self, enriched_records: List[VariantRecord]) -> List[VariantRecord]:
        """Further enrich records with drug resistance information."""
        drecen_records = []
        drug_names = {data.drug_name for data_list in self.var_ann_index.values() 
                    for data in data_list}

        for record in enriched_records:
            if not hasattr(record, 'var_ann'):
                self.logger.warning(f"Record at {record.chrom}:{record.pos} missing var_ann")  # Use instance logger
                continue

            # Initialize all drug attributes as None
            for drug in drug_names:
                setattr(record, drug.lower(), None)

            # Find matching drug data
            if record.var_ann in self.var_ann_index:
                matched = False
                for drug_data in self.var_ann_index[record.var_ann]:
                    drug_attr = drug_data.drug_name.lower()
                    current_value = getattr(record, drug_attr)
                    
                    # Create or update drug attribute
                    if current_value is None:
                        current_value = {}
                    current_value[record.var_ann] = {
                        'resistance_gene': drug_data.resistance_gene,
                        'resistance_mutation': drug_data.resistance_mutation,
                        'mutation_effect': drug_data.mutation_effect,
                        'confidence_grading': drug_data.confidence_grading,
                        'compared_to_version_1': drug_data.compared_to_version_1,
                        'genomic_position': drug_data.genomic_position,
                        'additional_grading_criteria_applied': drug_data.additional_grading_criteria_applied,
                        'comment': drug_data.comment
                    }
                    setattr(record, drug_attr, current_value)
                    matched = True

                if matched:
                    drecen_records.append(record)
                    self.logger.info(f"Enriched record at {record.chrom}:{record.pos} with drug resistance data")

        self.logger.info(f"Created {len(drecen_records)} DRECEN records")
        return drecen_records