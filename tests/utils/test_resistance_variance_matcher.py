import pytest
import pandas as pd
from pathlib import Path
import tempfile
from unittest.mock import Mock
import warnings
from profiler.utils.resistance_variant_matcher import ResistanceDB, VariantInfo

@pytest.fixture
def sample_tsv_file():
    """Create a temporary TSV file with test resistance data."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        f.write("chrom\tpos\tref\talt\tvar_ann\n")
        f.write("NC_000962.3\t70\tC\tT\tdnaA_p.Pro24Ser\n")
        f.write("NC_000962.3\t574702\tA\tT\tmshA_c.-646A>T\n")
        f.write("NC_000962.3\t1253460\tC\tG\tRv1129c_p.Val359Leu\n")
        f.write("NC_000962.3\t1406744\tCGG\tACT\tRv1258c_p.Pro199Ser\n")
    yield f.name
    Path(f.name).unlink()

@pytest.fixture
def mock_vcf_records():
    """Create mock VCF records for testing."""
    record1 = Mock()
    record1.chrom = "NC_000962.3"
    record1.pos = 70
    record1.ref = "C"
    record1.alts = ("T",)

    record2 = Mock()
    record2.chrom = "NC_000962.3"
    record2.pos = 574702
    record2.ref = "A"
    record2.alts = ("T",)

    record3 = Mock()
    record3.chrom = "NC_000962.3"
    record3.pos = 1406744
    record3.ref = "CGG"
    record3.alts = ("ACT",)

    return [record1, record2, record3]

class TestResistanceDB:
    def test_init_and_load(self, sample_tsv_file):
        """Test initialization and TSV loading."""
        db = ResistanceDB(sample_tsv_file)
        assert len(db.variants) == 4
        assert isinstance(db.variants, dict)
        
        # Check specific variant
        key = ("NC_000962.3", 70, "C", "T")
        assert key in db.variants
        assert isinstance(db.variants[key], VariantInfo)
        assert db.variants[key].var_ann == "dnaA_p.Pro24Ser"

        key = ("NC_000962.3", 1406744, "CGG", "ACT")
        assert key in db.variants
        assert isinstance(db.variants[key], VariantInfo)
        assert db.variants[key].var_ann == "Rv1258c_p.Pro199Ser"

    def test_missing_columns(self):
        """Test handling of TSV with missing required columns."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv') as f:
            f.write("chrom\tpos\tref\n")  # Missing columns
            f.write("chr1\t1000\tA\n")
            f.flush()
            
            with pytest.raises(ValueError) as exc_info:
                ResistanceDB(f.name)
            assert "missing required columns" in str(exc_info.value)

    def test_get_sample_variants(self, sample_tsv_file, mock_vcf_records):
        """Test variant matching functionality."""
        db = ResistanceDB(sample_tsv_file)
        matches = db.get_matching_variants(mock_vcf_records)
        
        assert len(matches) == 3
        assert "NC_000962.3_70_C_T" in matches
        assert "NC_000962.3_1406744_CGG_ACT" in matches
        
        # Check specific match details
        variant = matches["NC_000962.3_1406744_CGG_ACT"]
        assert variant.chrom == "NC_000962.3"
        assert variant.pos == 1406744
        assert variant.var_ann == "Rv1258c_p.Pro199Ser"

    def test_cache_functionality(self, sample_tsv_file, mock_vcf_records):
        """Test caching of results."""
        db = ResistanceDB(sample_tsv_file)
        
        # First query
        matches1 = db.get_matching_variants(mock_vcf_records)
        assert len(db._cached_matches) > 0
        
        # Second query should use cache
        matches2 = db.get_matching_variants(mock_vcf_records)
        assert matches1 == matches2  # Should be identical
        
        # Clear ca
        db.clear_cache()
        assert len(db._cached_matches) == 0
        

    def test_multi_allelic_warning(self, sample_tsv_file):
        """Test warning for multi-allelic variants."""
        db = ResistanceDB(sample_tsv_file)
        
        # Create mock record with multiple alternates
        multi_record = Mock()
        multi_record.chrom = "NC_000962.3"
        multi_record.pos = 70
        multi_record.ref = "C"
        multi_record.alts = ("T", "G")
        
        with pytest.warns(UserWarning) as warning_info:
            db.get_matching_variants([multi_record])
        
        assert len(warning_info) == 1
        assert "Multi-allelic variant found" in str(warning_info[0].message)

    def test_no_matches(self, sample_tsv_file):
        """Test behavior when no variants match."""
        db = ResistanceDB(sample_tsv_file)
        
        # Create mock record that won't match
        no_match_record = Mock()
        no_match_record.chrom = "chr3"
        no_match_record.pos = 5000
        no_match_record.ref = "A"
        no_match_record.alts = ("T",)
        
        matches = db.get_matching_variants([no_match_record])
        assert len(matches) == 0

    def test_variant_info_dataclass(self):
        """Test VariantInfo dataclass functionality."""
        variant = VariantInfo(
            chrom="NC_000962.3",
            pos=70,
            ref="C",
            alt="T",
            var_ann="dnaA_p.Pro24Ser"
        )
        
        assert variant.chrom == "NC_000962.3"
        assert variant.pos == 70
        assert variant.ref == "C"
        assert variant.alt == "T"
        assert variant.var_ann == "dnaA_p.Pro24Ser"