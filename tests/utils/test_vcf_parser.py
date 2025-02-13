import pytest
from unittest.mock import patch, mock_open, MagicMock
from profiler.utils.vcf_parser import VCFNormalizer
import tempfile
import pysam
import os

@pytest.fixture
def sample_vcf_files():
    return ["tests/test_data/single_sample_bcf/sample_116.vcf", 
            "tests/test_data/multi_sample_bcf/merged_samples_variants.vcf",
            "tests/test_data/minos_single_sample/sample_01_minos.vcf.gz",
            "tests/test_data/minos_single_sample/sample_287_minos.vcf.gz",
            "tests/test_data/gatk_multi_sample/gatk.vcf"]

@pytest.fixture
def empty_vcf_file():
    return "tests/test_data/empty.vcf"

@pytest.fixture
def unexpected_format_vcf_file():
    return "tests/test_data/unexpected_format.vcf"

@pytest.fixture
def non_vcf_file():
    return "tests/test_data/not_a_vcf.tsv"

class TestVCFNormalizer:
    @pytest.fixture
    def create_test_vcf(self):
        """Create a temporary VCF file and its bgzipped version."""
        def _create_vcf(content):
            with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as temp_vcf:
                temp_vcf.write(content.encode())
                temp_vcf_name = temp_vcf.name
            
            # Create bgzipped version
            out_fname = temp_vcf_name + '.gz'
            pysam.tabix_compress(temp_vcf_name, out_fname, force=True)
            pysam.tabix_index(out_fname, preset="vcf", force=True)
            
            return temp_vcf_name, out_fname
        return _create_vcf

    def test_split_multi_allelic_sites(self, create_test_vcf):
        """Test splitting of multi-allelic sites with genotype ordering."""
        vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
1\t100\t.\tA\tC,G\t.\t.\t.\tGT\t0/1
1\t200\t.\tA\tC,G\t.\t.\t.\tGT\t1/2
1\t300\t.\tA\tC,G\t.\t.\t.\tGT\t2/1
"""
        try:
            temp_vcf, gz_vcf = create_test_vcf(vcf_content)
            normalizer = VCFNormalizer([gz_vcf])
            normalized_vcf = normalizer.normalize()

            # Test pos=100: REF/ALT1
            assert normalized_vcf[0].samples["sample1"]["GT"] == (0, 1)  # A/C
            assert normalized_vcf[1].samples["sample1"]["GT"] == (0, 0)  # A/A for G

            # Test pos=200: ALT1/ALT2
            assert normalized_vcf[2].samples["sample1"]["GT"] == (0, 1)  # A/C
            assert normalized_vcf[3].samples["sample1"]["GT"] == (0, 1)  # A/G

            # Test pos=300: ALT2/ALT1
            assert normalized_vcf[4].samples["sample1"]["GT"] == (0, 1)  # A/C
            assert normalized_vcf[5].samples["sample1"]["GT"] == (0, 1)  # A/G

        finally:
            for f in [temp_vcf, gz_vcf, gz_vcf + '.tbi']:
                if os.path.exists(f):
                    os.unlink(f)

    def test_info_field_handling(self, create_test_vcf):
        """Test INFO fields with different Number specifications."""
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
1\t100\t.\tA\tT,C\t.\t.\tAF=0.1,0.2;DP=100;AD=50,30,20\tGT\t0/1
"""
        try:
            temp_vcf, gz_vcf = create_test_vcf(vcf_content)
            normalizer = VCFNormalizer([gz_vcf])
            normalized_vcf = normalizer.normalize()

            # Test Number=A fields with floating point comparison
            assert pytest.approx(normalized_vcf[0].info["AF"][0], rel=1e-6) == 0.1
            assert pytest.approx(normalized_vcf[1].info["AF"][0], rel=1e-6) == 0.2
            
            # Test Number=1 fields (single value)
            assert normalized_vcf[0].info["DP"] == 100
            assert normalized_vcf[1].info["DP"] == 100
            
            # Test Number=R fields (REF + ALT values)
            assert normalized_vcf[0].info["AD"] == (50, 30)
            assert normalized_vcf[1].info["AD"] == (50, 20)

        finally:
            for f in [temp_vcf, gz_vcf, gz_vcf + '.tbi']:
                if os.path.exists(f):
                    os.unlink(f)

    def test_format_field_handling(self, create_test_vcf):
        """Test FORMAT fields with different Number specifications."""
        vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
1\t100\t.\tA\tT,C\t.\t.\t.\tGT:AD:GQ\t0/1:50,30,20:99
"""
        try:
            temp_vcf, gz_vcf = create_test_vcf(vcf_content)
            normalizer = VCFNormalizer([gz_vcf])
            normalized_vcf = normalizer.normalize()

            # Test FORMAT fields
            sample1_data = normalized_vcf[0].samples["sample1"]
            assert sample1_data["GT"] == (0, 1)
            assert sample1_data["AD"] == (50, 30)
            assert sample1_data["GQ"] == 99

        finally:
            for f in [temp_vcf, gz_vcf, gz_vcf + '.tbi']:
                if os.path.exists(f):
                    os.unlink(f)

    def test_invalid_files(self):
        """Test handling of invalid VCF files."""
        # Test empty file
        with tempfile.NamedTemporaryFile(suffix='.vcf') as empty_file:
            normalizer = VCFNormalizer([empty_file.name])
            with pytest.raises(ValueError) as excinfo:
                normalizer.normalize()
            error_msg = str(excinfo.value)
            assert any(msg in error_msg for msg in [
                "Invalid VCF file format",
                "is it VCF/BCF format?",
                "no CHROM line found"
            ])

        # Test non-VCF file
        with tempfile.NamedTemporaryFile(suffix='.txt') as wrong_file:
            wrong_file.write(b"This is not a VCF file")
            wrong_file.flush()
            normalizer = VCFNormalizer([wrong_file.name])
            with pytest.raises(ValueError) as excinfo:
                normalizer.normalize()
            error_msg = str(excinfo.value)
            assert any(msg in error_msg for msg in [
                "Invalid VCF file format",
                "is it VCF/BCF format?",
                "no CHROM line found"
            ])

    def test_empty_vcf_file(self, empty_vcf_file):
        normalizer = VCFNormalizer([empty_vcf_file])
        with pytest.raises(Exception) as excinfo:
            normalizer.validate_and_read_vcf_files()
        assert "Invalid VCF file" in str(excinfo.value)

    def test_unexpected_format_vcf_file(self, unexpected_format_vcf_file):
        normalizer = VCFNormalizer([unexpected_format_vcf_file])
        with pytest.raises(Exception) as excinfo:
            normalizer.validate_and_read_vcf_files()
        assert "Invalid VCF file format" in str(excinfo.value)

    def test_non_vcf_file(self, non_vcf_file):
        normalizer = VCFNormalizer([non_vcf_file])
        with pytest.raises(Exception) as excinfo:
            normalizer.validate_and_read_vcf_files()
        assert "Invalid VCF file format" in str(excinfo.value)