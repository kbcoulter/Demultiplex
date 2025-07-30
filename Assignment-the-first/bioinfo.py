#!/usr/bin/env python

# Author: Kaden Coulter kcoulter@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
Functions with tests to ensure they are functioning properly are included below.
This file lives in its own bioinfo repo/module but may be copied elsewhere.'''

__version__ = "0.4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning
                            # 2025 July 17

DNA_bases = "ATGCatgc"
RNA_bases = "AUGCaugc"

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return (ord(letter) - 33)

def qual_score(phred_score: str) -> float:
    '''Accepts a phred_score string as input, returns the average quality score of the entire phred string'''
    score_holder = 0
    len_score = len(phred_score)
    for char in (phred_score):
        score_holder += (convert_phred(char))
    return (score_holder/len_score)

def validate_base_seq(seq: str, RNAflag: bool=False) -> bool:
    '''This function takes a string and an optional bool to denote RNA or DNA (DNA Standard with False). Returns True if string is composed
    of only A's, T's (or U's if RNAflag), G's, C's. False otherwise. Case insensitive.'''
    seq = seq.upper()
    return len(seq) == seq.count("A") + seq.count("U" if RNAflag else "T") + seq.count("G") + seq.count("C")

def gc_content(DNA:str) -> float:
    '''Takes in a DNA (or RNA) string and returns GC content of a DNA or RNA sequence as a float decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - Are you sure you used a valid sequence?"
    DNA = DNA.upper()
    Gs = DNA.count("G")
    Cs = DNA.count("C")
    return (Gs+Cs)/len(DNA)

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    if len(lst) % 2 == 0:
        return (lst[int(len(lst)/2)] + lst[int(len(lst)/2)-1]) / 2
    else:
        return (lst[int((len(lst)-1)/2)])

def oneline_fasta(file: str):
    '''Accepts a fasta file, returns a fasta file where all sequences are single lines named oneline_<file>.fa'''
    output_file = (f"oneline_{file}.fa")
    with open(output_file, "w") as new:
        header_tracker = 0
        with open(file, "r") as old:
            for line in old:
                if line[0] == ">" and header_tracker == 0:
                    header_tracker += 1
                    new.write(f"{line.strip()}\n")
                elif line[0] == ">":
                    new.write(f"\n{line.strip()}\n")
                else:
                    new.write(line.strip())
        new.write("\n")

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    assert calc_median([1,2,100]) == 2, "calc_median function does not work for odd length list"
    assert calc_median([1,2]) == 1.5, "calc_median function does not work for even length list"
    assert calc_median([1,1,1,1,1,1,1,1,1,5000]) == 1
    assert calc_median([1,2,3,4,5,6,7,13]) == 4.5
    print("Your calc_median function is working! Hooray!")

    assert qual_score("EEE") == 36
    assert qual_score("#I") == 21
    assert qual_score("EJ") == 38.5
    print("Your qual_score function is working! Great.")

    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("ATKXM") == False, "Validate base seq does not work with non CATGU chars"
    print("Your validate_base_seq passed DNA and RNA tests. Nice!")

    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATCGAT") == 0.5
    assert gc_content("G") == 1
    assert gc_content("C") == 1
    print("This GC content seems correct, gc_content looks like it works. Nice!")