from dataclasses import dataclass
from typing import Dict, Tuple, Optional
import pandas as pd
import warnings

@dataclass
class VariantInfo:
    """
    Container for variant information from resistance database.
    
    Attributes:
        chrom (str): Chromosome identifier
        pos (int): Position on chromosome
        ref (str): Reference allele
        alt (str): Alternate allele
        var_ann (str): Variant annotation from resistance database
    """
    chrom: str
    pos: int
    ref: str
    alt: str
    var_ann: str

class ResistanceDB:
    """
    Database for efficient lookup of variants associated with drug resistance.
    
    Provides fast access to resistance variant information and caches results
    for repeated queries. Uses two-step position matching for efficiency.
    """
    
    def __init__(self, tsv_path: str):
        """
        Initialize resistance database from TSV file.
        
        Args:
            tsv_path (str): Path to TSV file containing resistance data
                        Expected columns: chrom, pos, ref, alt, var_ann
        
        Raises:
            ValueError: If TSV file is missing required columns
        """
        self.variants: Dict[Tuple[str, int, str, str], VariantInfo] = {}
        self._cached_matches: Dict[str, VariantInfo] = {}
        self._load_tsv(tsv_path)

    def _load_tsv(self, tsv_path: str) -> None:
        """
        Load and index TSV data for efficient querying.
        
        Args:
            tsv_path (str): Path to TSV file
            
        Creates a dictionary with compound keys (chrom, pos, ref, alt) for O(1) lookups.
        """
        df = pd.read_csv(tsv_path, 
                        sep='\t',
                        dtype={
                            'chrom': str,
                            'pos': int,
                            'ref': str,
                            'alt': str,
                            'var_ann': str
                        })
        
        required_columns = {'chrom', 'pos', 'ref', 'alt', 'var_ann'}
        if not required_columns.issubset(df.columns):
            missing = required_columns - set(df.columns)
            raise ValueError(f"TSV file missing required columns: {missing}")
        
        for _, row in df.iterrows():
            key = (row['chrom'], row['pos'], row['ref'], row['alt'])
            self.variants[key] = VariantInfo(**row.to_dict())

    def get_matching_variants(self, vcf_records: list) -> Dict[str, VariantInfo]:
        """
        Get variants that match the resistance database.
        
        Args:
            vcf_records (list): List of VCF records (expected to be normalized/bi-allelic)
            
        Returns:
            Dict[str, VariantInfo]: Dictionary of matched variants with their annotations
            
        Note:
            Issues warning if multi-allelic variants are found, suggesting normalization
        """
        # Return cached results if available
        if self._cached_matches:
            return self._cached_matches

        results = {}
        
        # Create sets for fast position matching
        positions = {(record.chrom, record.pos) for record in vcf_records}
        potential_positions = {(v.chrom, v.pos) for v in self.variants.values()}
        matching_positions = positions & potential_positions

        # Process matching positions
        for record in vcf_records:
            # Warn if record is multi-allelic
            if len(record.alts) > 1:
                warnings.warn(
                    f"Multi-allelic variant found at {record.chrom}:{record.pos}. "
                    "Records should be normalized to bi-allelic variants first.",
                    UserWarning
                )

            if (record.chrom, record.pos) in matching_positions:
                for alt in record.alts:
                    key = (record.chrom, record.pos, record.ref, alt)
                    if key in self.variants:
                        var_key = f"{record.chrom}_{record.pos}_{record.ref}_{alt}"
                        results[var_key] = self.variants[key]

        # Cache results for future queries
        self._cached_matches = results
        return results

    def clear_cache(self) -> None:
        """
        Clear the sample cache to free memory.
        
        Args:
            sample_id (str, optional): Specific sample to clear from cache.
                                    If None, clears entire cache.
        """
        self._cached_matches.clear()