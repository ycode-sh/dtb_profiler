import pytest
from unittest.mock import Mock, patch
import logging
from profiler.utils.variant_record_enricher import VariantRecordEnricher
from profiler.utils.resistance_variant_matcher import VariantInfo

@pytest.fixture
def mock_variant_records():
    """Create mock VCF records for testing."""
    # Valid bi-allelic record
    record1 = Mock()
    record1.chrom = "NC_000962.3"
    record1.pos = 70
    record1.ref = "C"
    record1.alts = ("T",)

    # Multi-allelic record (invalid)
    record2 = Mock()
    record2.chrom = "NC_000962.3"
    record2.pos = 574702
    record2.ref = "A"
    record2.alts = ("T", "G")

    # No alternate alleles (invalid)
    record3 = Mock()
    record3.chrom = "NC_000962.3"
    record3.pos = 1253460
    record3.ref = "C"
    record3.alts = ()

    # Complex variant
    record4 = Mock()
    record4.chrom = "NC_000962.3"
    record4.pos = 1406744
    record4.ref = "CGG"
    record4.alts = ("ACT",)

    return [record1, record2, record3, record4]

@pytest.fixture
def mock_resistance_matches():
    """Create mock resistance matches."""
    return {
        "NC_000962.3_70_C_T": VariantInfo(
            chrom="NC_000962.3",
            pos=70,
            ref="C",
            alt="T",
            var_ann="dnaA_p.Pro24Ser"
        ),
        "NC_000962.3_1406744_CGG_ACT": VariantInfo(
            chrom="NC_000962.3",
            pos=1406744,
            ref="CGG",
            alt="ACT",
            var_ann="Rv1258c_p.Pro199Ser"
        )
    }

class TestVariantRecordEnricher:
    @patch('profiler.utils.variant_record_enricher.logger')
    def test_validate_bi_allelic(self, mock_logger, mock_variant_records):
        """Test validation of bi-allelic records."""
        # Valid record
        assert VariantRecordEnricher._validate_bi_allelic(mock_variant_records[0]) == "T"
        
        # Multi-allelic record
        result = VariantRecordEnricher._validate_bi_allelic(mock_variant_records[1])
        assert result is None
        mock_logger.error.assert_called_once()

        # No alternates
        result = VariantRecordEnricher._validate_bi_allelic(mock_variant_records[2])
        assert result is None
        mock_logger.warning.assert_called_once()

    def test_enrich_records(self, mock_variant_records, mock_resistance_matches):
        """Test enrichment of records with resistance annotations."""
        enriched = VariantRecordEnricher.enrich_records(
            mock_variant_records,
            mock_resistance_matches
        )

        # Should only enrich valid bi-allelic records with matches
        assert len(enriched) == 2
        
        # Check first enriched record
        assert hasattr(enriched[0], 'var_ann')
        assert enriched[0].var_ann == "dnaA_p.Pro24Ser"
        assert enriched[0].chrom == "NC_000962.3"
        assert enriched[0].pos == 70
        
        # Check complex variant
        assert hasattr(enriched[1], 'var_ann')
        assert enriched[1].var_ann == "Rv1258c_p.Pro199Ser"
        assert enriched[1].ref == "CGG"
        assert enriched[1].alts[0] == "ACT"

    @patch('profiler.utils.variant_record_enricher.logger')
    def test_field_mismatch(self, mock_logger, mock_variant_records):
        """Test handling of field mismatches."""
        # Use only the valid bi-allelic record
        valid_record = mock_variant_records[0]  # Using just the first record
        
        mismatched_matches = {
            "NC_000962.3_70_C_T": VariantInfo(
                chrom="NC_000962.3",
                pos=71,  # Mismatched position
                ref="C",
                alt="T",
                var_ann="dnaA_p.Pro24Ser"
            )
        }

        enriched = VariantRecordEnricher.enrich_records(
            [valid_record],  # Pass only the single record
            mismatched_matches
        )
        
        assert len(enriched) == 0
        mock_logger.warning.assert_called_once_with(
            "Field mismatch for NC_000962.3_70_C_T. Record: NC_000962.3:70 C>T"
        )

    def test_no_matches(self, mock_variant_records):
        """Test behavior with no matching variants."""
        enriched = VariantRecordEnricher.enrich_records(
            mock_variant_records,
            {}  # Empty matches
        )
        assert len(enriched) == 0

    def test_enriched_record_attributes(self, mock_variant_records, mock_resistance_matches):
        """Test that enriched records maintain original attributes."""
        # Add some additional attributes to mock record
        mock_variant_records[0].qual = 30
        mock_variant_records[0].filter = {"PASS"}
        mock_variant_records[0].info = {"DP": 100}

        enriched = VariantRecordEnricher.enrich_records(
            mock_variant_records,
            mock_resistance_matches
        )

        assert len(enriched) > 0
        assert enriched[0].qual == 30
        assert enriched[0].filter == {"PASS"}
        assert enriched[0].info["DP"] == 100
        assert hasattr(enriched[0], 'var_ann')

