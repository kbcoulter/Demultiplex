#!/usr/bin/env python

from demultiplex import reverse_comp, collect_record, convert_phred, mean_qual_score
from bioinfo import convert_phred

def test_reverse_comp():
    assert reverse_comp("ATCG") == "CGAT", "Reverse complement failed"
    assert reverse_comp("atcg") == "CGAT", "Reverse complement failed"
    assert reverse_comp("atcg ") == "CGAT", "Reverse complement failed"
    assert reverse_comp("\nAtCg ") == "CGAT", "Reverse complement failed"
    assert reverse_comp("NATCG") == "CGATN", "Reverse complement failed"

def test_convert_phred_basic():
    assert convert_phred("!") == 0, "Convert Phred failed"
    assert convert_phred("I") == 40, "Convert Phred failed"
    assert convert_phred() == None, "Convert Phred failed"

def test_mean_qual_basic():
    assert mean_qual_score("!I") == 20, "Mean qual failed"
    assert mean_qual_score("IIII") == 40, "Mean qual failed"
    assert mean_qual_score("IIIIIIII") == 40, "Mean qual failed"
    # Empty Case Covered, as we are using convert_phred

def test_collect_record():
    with open("/TEST-input_FASTQ/S1_TEST_R2_001.fastq", "rt") as test1, open("/TEST-input_FASTQ/S1_TEST_R1_001.fastq", "rt") as test2:
        result1 = collect_record("/TEST-input_FASTQ/S1_TEST_R2_001.fastq", True)
        result2 = collect_record("/TEST-input_FASTQ/S1_TEST_R1_001.fastq", False, "dog")
        expected1 = "GTAGCGTA"
        expected2 = [
            "@seq1 dog",
            "CATCATCATCAT",
            "+",
            "************"]
        assert result1 == expected1, f"Got\t{result1}"
        assert result2 == expected2, f"Got\t{result2}"

if __name__ == "__main__":
    test_reverse_comp()
    test_convert_phred_basic()
    test_mean_qual_basic()
    test_collect_record()
    print("All Functions in demultiplex.py work as expected!")