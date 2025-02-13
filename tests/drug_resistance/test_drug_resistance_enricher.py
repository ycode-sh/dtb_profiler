import pytest
from pathlib import Path
import pandas as pd
from unittest.mock import Mock, patch
import tempfile
from profiler.utils.drug_resistance_enricher import DrugResistanceEnricher, DrugData

MUTATION_EFFECTS = {
    'feature_ablation', 'frameshift', 'inframe_deletion', 'inframe_insertion',
    'initiator_codon_variant', 'lof', 'missense_variant', 'start_lost',
    'stop_gained', 'stop_lost', 'stop_retained_variant', 'synonymous_variant',
    'upstream_gene_variant', 'non_coding_transcript_exon_variant'
}

CONFIDENCE_GRADES = {
    'Assoc_w_R', 'Assoc_w_RI', 'Uncertain_significance',
    'Not_assoc_w_RI', 'Not_assoc_w_R'
}

COMPARISON_VALUES = {
    'DOWN from AWR to AWRI', 'DOWN from AwR to Uncertain',
    'DOWN from AWRI to Uncertain', 'DOWN from NotAWR to NotAWRI',
    'DOWN from NotAwR to Uncertain', 'DOWN from NotAWRI to Uncertain',
    'New AWR', 'New AWRI', 'New NotAWR', 'New NotAWRI', 'New Uncertain',
    'No change', 'Now listed', 'SWITCH from AwRI to NotAWRI',
    'UP from AwRI to AwR', 'UP from NotAWRI to NotAwR',
    'UP from Uncertain to AwR', 'UP from Uncertain to AwRI',
    'UP from Uncertain to NotAwR', 'UP from Uncertain to NotAWRI'
}

@pytest.fixture
def mock_drug_files():
    """Create temporary drug TSV files with comprehensive test data."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create isoniazid.tsv with various mutation effects and confidence grades
        isoniazid_data = pd.DataFrame({
            'drug_name': ['Isoniazid'] * 4,
            'resistance_gene': ['katG', 'inhA', 'fabG1', 'ahpC'],
            'resistance_mutation': ['p.Ser315Thr', 'c.-15C>T', 'deletion', 'insertion'],
            'var_ann': [
                'katG_p.Ser315Thr', 'inhA_c.-15C>T',
                'fabG1_deletion', 'ahpC_insertion'
            ],
            'mutation_effect': [
                'missense_variant', 'upstream_gene_variant',
                'feature_ablation', 'frameshift'
            ],
            'confidence_grading': [
                'Assoc_w_R', 'Assoc_w_RI',
                'Uncertain_significance', 'Not_assoc_w_RI'
            ],
            'compared_to_version_1': [
                'Now listed', 'UP from Uncertain to AwR',
                'New Uncertain', 'DOWN from AWR to AWRI'
            ],
            'genomic_position': ['2155168', '1674481', '1673425', '2726105'],
            'additional_grading_criteria_applied': ['Yes', 'No', '', 'Yes'],
            'comment': ['High confidence', 'Promoter mutation', '', 'Complex change']
        })
        isoniazid_data.to_csv(f"{temp_dir}/isoniazid.tsv", sep='\t', index=False)

        # Create rifampicin.tsv with different combinations
        rifampicin_data = pd.DataFrame({
            'drug_name': ['Rifampicin'] * 3,
            'resistance_gene': ['rpoB'] * 3,
            'resistance_mutation': [
                'p.Ser450Leu', 'p.His445Tyr', 'c.1349_1350delGC'
            ],
            'var_ann': [
                'rpoB_p.Ser450Leu', 'rpoB_p.His445Tyr', 'rpoB_c.1349_1350delGC'
            ],
            'mutation_effect': [
                'missense_variant', 'stop_gained', 'inframe_deletion'
            ],
            'confidence_grading': [
                'Assoc_w_R', 'Uncertain_significance', 'Not_assoc_w_R'
            ],
            'compared_to_version_1': [
                'No change', 'New AWR', 'DOWN from NotAWR to NotAWRI'
            ],
            'genomic_position': ['761155', '761140', '761110'],
            'additional_grading_criteria_applied': ['Yes', '', 'No'],
            'comment': ['RRDR mutation', '', 'Complex deletion']
        })
        rifampicin_data.to_csv(f"{temp_dir}/rifampicin.tsv", sep='\t', index=False)

        yield Path(temp_dir)

@pytest.fixture
def mock_enriched_records():
    """Create mock enriched VCF records."""
    record1 = Mock()
    record1.chrom = "NC_000962.3"
    record1.pos = 2155168
    record1.ref = "C"
    record1.alts = ("T",)
    record1.var_ann = "katG_p.Ser315Thr"

    record2 = Mock()
    record2.chrom = "NC_000962.3"
    record2.pos = 761155
    record2.ref = "C"
    record2.alts = ("T",)
    record2.var_ann = "rpoB_p.Ser450Leu"

    return [record1, record2]

class TestDrugResistanceEnricher:
    """ PRE-DRECEN VALIDATION TESTS 
    
    a. test_init_and_load: Validates DrugResistanceEnricher initialization
    b. test_initial_records_have_var_ann: Ensures input records are properly annotated
    c. test_missing_required_columns_in_drug_files: Validates TSV file structure
    
    """

    def test_init_and_load(self, mock_drug_files):
        """Test initialization and drug data loading."""
        enricher = DrugResistanceEnricher(mock_drug_files)
        
        # Check var_ann index
        assert "katG_p.Ser315Thr" in enricher.var_ann_index
        assert "inhA_c.-15C>T" in enricher.var_ann_index
        assert "rpoB_p.Ser450Leu" in enricher.var_ann_index

        # Check drug data content
        isoniazid_data = enricher.var_ann_index["katG_p.Ser315Thr"][0]
        assert isoniazid_data.drug_name == "Isoniazid"
        assert isoniazid_data.confidence_grading == "Assoc_w_R"

    
    def test_initial_records_have_var_ann(self, mock_drug_files, mock_enriched_records):
        """Test that initial records have var_ann attribute before DRECEN enrichment."""
    
        # Verify all initial records have var_ann
        for record in mock_enriched_records:
            assert hasattr(record, 'var_ann'), f"Record at {record.chrom}:{record.pos} missing var_ann"
            assert record.var_ann is not None, f"Record at {record.chrom}:{record.pos} has None var_ann"
            assert isinstance(record.var_ann, str), f"Record at {record.chrom}:{record.pos} var_ann must be string"
    
    def test_missing_required_columns_in_drug_files(self, mock_drug_files):
        """Test handling of missing required columns."""
        bad_data = pd.DataFrame({
            'drug_name': ['Isoniazid'],
            'var_ann': ['katG_p.Ser315Thr']
        })
        bad_file = mock_drug_files / "bad.tsv"
        bad_data.to_csv(bad_file, sep='\t', index=False)

        with pytest.raises(ValueError) as exc_info:
            DrugResistanceEnricher(mock_drug_files)
        assert "Missing required column" in str(exc_info.value)

    """ DRUG DATA VALIDATION TESTS

    a. test_mutation_effect_values: Checks mutation effect values
    b. test_confidence_grading_values: Validates confidence grades
    c. test_comparison_values: Ensures version comparison values are valid

    """

    def test_mutation_effect_values(self, mock_drug_files):
        """Test validation of mutation effect values."""
        enricher = DrugResistanceEnricher(mock_drug_files)
        for data_list in enricher.var_ann_index.values():
            for data in data_list:
                assert data.mutation_effect in MUTATION_EFFECTS

    def test_confidence_grading_values(self, mock_drug_files):
        """Test validation of confidence grading values."""
        enricher = DrugResistanceEnricher(mock_drug_files)
        for data_list in enricher.var_ann_index.values():
            for data in data_list:
                assert data.confidence_grading in CONFIDENCE_GRADES

    def test_comparison_values(self, mock_drug_files):
        """Test validation of compared_to_version_1 values."""
        enricher = DrugResistanceEnricher(mock_drug_files)
        for data_list in enricher.var_ann_index.values():
            for data in data_list:
                assert data.compared_to_version_1 in COMPARISON_VALUES

    """" DRECEN ENRICHMENT TESTS
    test_basic_record_enrichment: Basic enrichment functionality
    test_multiple_variants_per_drug: Multiple variants handling
    test_drecen_record_attributes: VCF attribute preservation
    test_drug_attribute_completeness: Drug attribute handling
    """

    def test_basic_record_enrichment(self, mock_drug_files, mock_enriched_records):
        """Test basic record enrichment with drug resistance data."""
        enricher = DrugResistanceEnricher(mock_drug_files)
        drecen_records = enricher.enrich_records(mock_enriched_records)

        assert len(drecen_records) == 2

        # Check first record (Isoniazid resistance)
        record1 = drecen_records[0]
        assert hasattr(record1, 'isoniazid')
        assert record1.isoniazid is not None
        assert record1.rifampicin is None
        
        isoniazid_data = record1.isoniazid['katG_p.Ser315Thr']
        assert isoniazid_data['resistance_gene'] == 'katG'
        assert isoniazid_data['confidence_grading'] == 'Assoc_w_R'

        # Check second record (Rifampicin resistance)
        record2 = drecen_records[1]
        assert hasattr(record2, 'rifampicin')
        assert record2.rifampicin is not None
        assert record2.isoniazid is None

        rifampicin_data = record2.rifampicin['rpoB_p.Ser450Leu']
        assert rifampicin_data['resistance_gene'] == 'rpoB'
        assert rifampicin_data['confidence_grading'] == 'Assoc_w_R'

    def test_multiple_variants_per_drug(self, mock_drug_files):
        """Test handling of multiple variants for the same drug."""
        enricher = DrugResistanceEnricher(mock_drug_files)
        
        # Create record with multiple variants for isoniazid
        record1 = Mock(
            chrom="NC_000962.3",
            pos=2155168,
            ref="C",
            alt="T",
            var_ann="katG_p.Ser315Thr"
        )
        record2 = Mock(
            chrom="NC_000962.3",
            pos=1674481,
            ref="C",
            alt="T",
            var_ann="inhA_c.-15C>T"
        )
        
        # Test enrichment
        drecen_records = enricher.enrich_records([record1, record2])
        
        # Verify results
        assert len(drecen_records) == 2
        for record in drecen_records:
            assert record.isoniazid is not None
            assert len(record.isoniazid) >= 1  # At least one variant
        
        # Check specific variants
        isoniazid_variants = set()
        for record in drecen_records:
            isoniazid_variants.update(record.isoniazid.keys())
        
        assert "katG_p.Ser315Thr" in isoniazid_variants
        assert "inhA_c.-15C>T" in isoniazid_variants
    
    def test_drecen_record_attributes(self, mock_drug_files, mock_enriched_records):
        """Test preservation of original VCF attributes in DRECEN records."""
        # Add additional attributes to mock record
        record = mock_enriched_records[0]
        record.qual = 30
        record.filter = {"PASS"}
        record.info = {"DP": 100}
        record.format = {"GT": "1/1", "AD": "0,30"}
        
        enricher = DrugResistanceEnricher(mock_drug_files)
        drecen_records = enricher.enrich_records([record])
        
        # Verify original attributes preserved
        assert len(drecen_records) == 1
        drecen_record = drecen_records[0]
        
        # Check VCF attributes
        assert drecen_record.chrom == record.chrom
        assert drecen_record.pos == record.pos
        assert drecen_record.ref == record.ref
        assert drecen_record.alts == record.alts
        assert drecen_record.qual == record.qual
        assert drecen_record.filter == record.filter
        assert drecen_record.info == record.info
        assert drecen_record.format == record.format
        
        # Check drug resistance attributes
        assert hasattr(drecen_record, "isoniazid")
        assert drecen_record.isoniazid is not None
        assert "katG_p.Ser315Thr" in drecen_record.isoniazid

    def test_drug_attribute_completeness(self, mock_drug_files, mock_enriched_records):
        """Test that all drug attributes are present, even if None."""
        enricher = DrugResistanceEnricher(mock_drug_files)
        drecen_records = enricher.enrich_records(mock_enriched_records)
        
        # Get all possible drug names
        drug_names = {data.drug_name.lower() 
                    for data_list in enricher.var_ann_index.values() 
                    for data in data_list}
        
        # Check each record has all drug attributes
        for record in drecen_records:
            for drug in drug_names:
                assert hasattr(record, drug)
                # At least one drug should not be None
                assert any(getattr(record, d) is not None for d in drug_names)
    

    
    
    


    

