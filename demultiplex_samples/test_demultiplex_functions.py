#!/usr/bin/env python

from demultiplex import reverse_comp, collect_record, mean_qual_score
import bioinfo

def test_reverse_comp():
    assert reverse_comp("ATCG") == "CGAT", "Reverse complement failed"
    assert reverse_comp("atcg") == "CGAT", "Reverse complement failed"
    assert reverse_comp("atcg ") == "CGAT", "Reverse complement failed"
    assert reverse_comp("\nAtCg ") == "CGAT", "Reverse complement failed"
    assert reverse_comp("NATCG") == "CGATN", "Reverse complement failed"

def test_convert_phred_basic():
    assert bioinfo.convert_phred("!") == 0, "Convert Phred failed"
    assert bioinfo.convert_phred("I") == 40, "Convert Phred failed"
    assert bioinfo.convert_phred("%") == 4, "Convert Phred failed"
    assert bioinfo.convert_phred("5") == 20, "Convert Phred failed"
    #assert bioinfo.convert_phred("") == "", "Convert Phred failed"

def test_mean_qual_basic():
    assert mean_qual_score("!I") == 20, "Mean qual failed"
    assert mean_qual_score("IIII") == 40, "Mean qual failed"
    assert mean_qual_score("IIIIIIII") == 40, "Mean qual failed"
    assert mean_qual_score("HJ") == 40, "Mean qual failed"
    assert mean_qual_score("%") == 4, "Mean qual failed"


    # Empty Case Covered, as we are using convert_phred

def test_collect_record():
    with open("../TEST-input_FASTQ/S1_TEST_R2_001.fastq", "rt") as test1, open("../TEST-input_FASTQ/S1_TEST_R1_001.fastq", "rt") as test2:
        result1 = collect_record(test1, True)
        result2 = collect_record(test2, False, "dog")
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