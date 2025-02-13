import pysam
import gzip

"""
VCFNormalizer
============

A module for normalizing VCF (Variant Call Format) files.

Features
--------
* Splits multi-allelic variants into bi-allelic records
* Merges multiple VCF files while preserving sample information
* Validates VCF format compliance
* Handles both single and multi-sample VCFs

Usage
-----
```python
normalizer = VCFNormalizer(["sample1.vcf.gz", "sample2.vcf"])
normalized_records = normalizer.normalize()

Detailed feature

1. Multi-allelic Site Handling
Splits multi-allelic variants into bi-allelic records
Preserves genotype relationships
Maintains INFO and FORMAT field integrity
2. File Format Support
Handles both compressed (.gz) and uncompressed VCF files
Supports both single-sample and multi-sample VCFs
Validates VCF header and field formats
3. Sample Merging
Combines records from multiple single-sample VCFs
Preserves sample-specific information
Handles overlapping variants across samples

Main Methods
1. normalize()

def normalize(self):
    Main entry point for VCF normalization.
    
    Returns:
        list: Normalized VCF records with split multi-allelic sites
    
    Process:
    1. Categorizes input files (single/multi-sample)
    2. Splits multi-allelic sites
    3. Merges samples
    4. Validates final records
    

2. _split_multi_allelic_record()
def _split_multi_allelic_record(self, record):
    
    Splits multi-allelic sites into individual records.
    
    Args:
        record: VCF record with multiple alternate alleles
    
    Returns:
        list: Split bi-allelic records
        
    Handles:
    - INFO fields (Number=A, R, 1)
    - FORMAT fields
    - Genotype conversion

    
Data Handling
INFO Fields
Number=A: One value per alternate allele
Number=R: One value per allele (including reference)
Number=1: Single value shared across alleles
FORMAT Fields
GT: Genotype field (0/1 format after splitting)
AD: Allelic depth (split according to alleles)
Other FORMAT fields preserved according to their Number specification

Input Requirements
VCF Format
Must include mandatory headers (CHROM, POS, REF, ALT, FORMAT)
Properly formatted VCF structure
Valid genotype encodings
File Types
Standard VCF (.vcf)
Compressed VCF (.vcf.gz)
Properly indexed for compressed files (.tbi)
Output Format
Normalized Records
Bi-allelic variants only
Standardized genotype representation (0/1)
Preserved sample information
Validated field formats
Error Handling
Validates VCF format compliance
Checks field data types
Handles missing or malformed headers
Reports specific validation errors
Dependencies
pysam: VCF file handling
gzip: Compressed file support
Performance Considerations
Memory usage scales with number of variants and samples
Processes files sequentially
Creates temporary copies during splitting

"""
class VCFNormalizer:
    def __init__(self, vcf_files):
        self.vcf_files = vcf_files
        self.compulsory_fields = {"CHROM", "POS", "REF", "ALT", "FORMAT"}
        self.optional_fields = {"QUAL", "FILTER", "INFO"}
        self.merged_vcf_data = []
        self.is_multi_sample = False

    def validate_and_read_vcf_files(self):
        # Validate and read VCF files
        single_sample_vcfs = []
        for vcf_file in self.vcf_files:
            try:
                # Check if the file is gzipped
                if vcf_file.endswith('.gz'):
                    open_func = gzip.open
                    mode = 'rt'
                else:
                    open_func = open
                    mode = 'r'
                
                with open_func(vcf_file, mode) as f:
                    header_line = None
                    for line in f:
                        if line.startswith("#CHROM"):
                            header_line = line.strip().split('\t')
                            break
                    if header_line is None:
                        raise ValueError(f"Missing header line in VCF file: {vcf_file}")
                    self._validate_header(header_line)
            except Exception as e:
                raise ValueError(f"Invalid VCF file format: {vcf_file}") from e
            
            reader = pysam.VariantFile(vcf_file)
            if len(reader.header.samples) == 1:
                single_sample_vcfs.append(vcf_file)
            else:
                self.is_multi_sample = True
                for record in reader:
                    self.merged_vcf_data.append(record)
        
        if single_sample_vcfs:
            self._merge_single_sample_vcfs(single_sample_vcfs)

    def _validate_header(self, header_line):
        # Remove the '#' character from the first element if it exists
        if header_line[0].startswith("#"):
            header_line[0] = header_line[0][1:]
        
        # Validate the header line
        header_fields = set(header_line)
        
        # Fixed column headers that are always expected
        fixed_headers = {"CHROM", "POS", "ID", "REF", "ALT", "FORMAT"}
        
        # Validate compulsory fields
        for field in self.compulsory_fields:
            if field not in header_fields:
                raise ValueError(f"Missing compulsory field: {field}")
        
        # Identify sample names by excluding fixed headers and optional fields
        sample_names = header_fields - fixed_headers - self.optional_fields
        
        # Validate unknown fields, excluding sample names
        for field in header_fields:
            if field not in fixed_headers and field not in self.optional_fields and field not in sample_names:
                raise ValueError(f"Unknown field: {field}")

    def validate_each_field(self):
        # Validate the data in each field
        for record in self.merged_vcf_data:
            print(record)
            if not isinstance(record.chrom, str):
                raise ValueError("Invalid CHROM field")
            if not isinstance(record.pos, int):
                raise ValueError("Invalid POS field")
            if not isinstance(record.ref, str):
                raise ValueError("Invalid REF field")
            if not all(isinstance(alt, str) for alt in record.alts):
                raise ValueError("Invalid ALT field")
            #if "FORMAT" in record and not isinstance(record.format, str):
            #    raise ValueError("Invalid FORMAT field")


    def _split_multi_allelic_record(self, record):
        """Split a single multi-allelic record into multiple bi-allelic records."""
        if len(record.alts) <= 1:
            return [record]

        split_records = []
        for idx, alt in enumerate(record.alts):
            new_record = record.copy()
            new_record.alts = (alt,)  
        
            # Handle INFO fields per alternate allele
            for info_key in new_record.info.keys():
                info_val = record.info[info_key]
                if isinstance(info_val, tuple):
                    if len(info_val) == len(record.alts):  # Number=A
                        new_record.info[info_key] = (info_val[idx],)
                    elif len(info_val) == len(record.alts) + 1:  # Number=R
                        new_record.info[info_key] = (info_val[0], info_val[idx + 1])
        
            # Handle FORMAT fields per sample
            for sample in new_record.samples:
                self._adjust_sample_format_fields(new_record, record, sample, idx)
        
            split_records.append(new_record)
    
        return split_records

    def _adjust_genotype(self, gt_value, alt_idx):
        """
    Adjust genotype values when splitting multi-allelic sites.
    
    Args:
        gt_value (tuple): Original genotype values
        alt_idx (int): Index of current alternate allele (0-based)
    
    Returns:
        tuple: Adjusted genotype values (0/1 format)
    
    Examples:
        For REF=A, ALT=C,G:
        Original GT=1/2 (C/G):
            - For C (alt_idx=0): becomes 0/1 (A/C)
            - For G (alt_idx=1): becomes 0/1 (A/G)
    """
        target_alt = alt_idx + 1
        new_gt = []
    
        # Preserve phasing information
        is_phased = hasattr(gt_value, 'phased') and gt_value.phased
    
        #Always convert to reference-first format (0/1)
        for allele in gt_value:
            if allele is None:
                new_gt.append(None)
            elif allele == target_alt:
                new_gt.append(1)  # Current alt becomes 1
            else:
                new_gt.append(0)  # Reference and other alts become 0
    
        # Always ensure reference allele comes first
        if len(new_gt) == 2 and new_gt[0] == 1 and new_gt[1] == 0:
            new_gt = [0, 1]
        # Create tuple with phasing information if available
        if is_phased:
            result = pysam.VariantRecordSample.phased_gt(*new_gt)
        else:
            result = tuple(new_gt)
    
        return result
    


    def _adjust_sample_format_fields(self, new_record, original_record, sample, alt_idx):
        """Adjust FORMAT fields for a given sample when splitting multi-allelic sites."""
        sample_data = new_record.samples[sample]
        for format_key in sample_data.keys():
            format_val = original_record.samples[sample][format_key]
            if isinstance(format_val, tuple):
                if format_key == 'GT':
                    sample_data[format_key] = self._adjust_genotype(format_val, alt_idx)
                elif len(format_val) == len(original_record.alts):  # Number=A
                    sample_data[format_key] = (format_val[alt_idx],)
                elif len(format_val) == len(original_record.alts) + 1:  # Number=R
                    sample_data[format_key] = (format_val[0], format_val[alt_idx + 1])

    def _merge_single_sample_vcfs(self, single_sample_vcfs):
        """Merge pre-split VCF records."""
        merged_records = {}
        
        for vcf_file in single_sample_vcfs:
            reader = pysam.VariantFile(vcf_file)
            for record in reader:
                # Split the record first if it's multi-allelic
                split_records = self._split_multi_allelic_record(record)
                
                for split_record in split_records:
                    # Key now includes only single alt
                    key = (split_record.chrom, split_record.pos, split_record.ref, split_record.alts[0])
                    
                    if key not in merged_records:
                        merged_records[key] = split_record
                    else:
                        # Merge samples
                        if not hasattr(merged_records[key], 'samples'):
                            merged_records[key].samples = {}
                        merged_records[key].samples.update(split_record.samples)
        
        return list(merged_records.values())

    def normalize(self):
        """Main method to normalize the VCF files."""
        single_sample_vcfs = []
        multi_sample_records = []

        # First pass: categorize files and handle multi-sample VCFs
        for vcf_file in self.vcf_files:
            reader = pysam.VariantFile(vcf_file)
            if len(reader.header.samples) == 1:
                single_sample_vcfs.append(vcf_file)
            else:
                self.is_multi_sample = True
                for record in reader:
                    split_records = self._split_multi_allelic_record(record)
                    multi_sample_records.extend(split_records)

        # Handle single sample VCFs
        if single_sample_vcfs:
            self.merged_vcf_data = self._merge_single_sample_vcfs(single_sample_vcfs)
        
        # Combine with multi-sample records if any
        if multi_sample_records:
            self.merged_vcf_data.extend(multi_sample_records)

        self.validate_each_field()
        return self.merged_vcf_data