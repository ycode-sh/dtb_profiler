import logging
from typing import List, Dict, Optional
from pysam import VariantRecord
from .resistance_variant_matcher import VariantInfo

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class VariantRecordEnricher:
    """Enriches VCF records with resistance annotations."""
    
    @staticmethod
    def _validate_bi_allelic(record: VariantRecord) -> Optional[str]:
        """
        Validate that a record is bi-allelic.
        
        Args:
            record: VCF record to validate
            
        Returns:
            str: Single alternate allele if valid, None if invalid
        """
        if not record.alts:
            logger.warning(f"Record at {record.chrom}:{record.pos} has no alternate alleles")
            return None
        
        if len(record.alts) > 1:
            logger.error(
                f"Multi-allelic variant found at {record.chrom}:{record.pos}. "
                "Expected normalized bi-allelic record."
            )
            return None
            
        return record.alts[0]
    
    @staticmethod
    def enrich_records(vcf_records: List[VariantRecord], 
                    resistance_matches: Dict[str, VariantInfo]) -> List[VariantRecord]:
        """
        Enrich VCF records with resistance annotations.
        
        Args:
            vcf_records: List of normalized VCF records
            resistance_matches: Dictionary of resistance matches
            
        Returns:
            List of enriched VCF records that match resistance variants
        """
        enriched_records = []
        
        for record in vcf_records:
            alt = VariantRecordEnricher._validate_bi_allelic(record)
            if not alt:
                continue
                
            key = f"{record.chrom}_{record.pos}_{record.ref}_{alt}"
            
            if key in resistance_matches:
                variant_info = resistance_matches[key]
                
                # Verify fields match
                if (record.chrom == variant_info.chrom and
                    record.pos == variant_info.pos and
                    record.ref == variant_info.ref and
                    alt == variant_info.alt):
                    
                    # Add var_ann as new attribute
                    record.var_ann = variant_info.var_ann
                    enriched_records.append(record)
                    
                    logger.info(
                        f"Enriched record at {key} with annotation: {variant_info.var_ann}"
                    )
                else:
                    logger.warning(
                        f"Field mismatch for {key}. Record: "
                        f"{record.chrom}:{record.pos} {record.ref}>{alt}"
                    )
        
        logger.info(f"Enriched {len(enriched_records)} records with resistance annotations")
        return enriched_records

        