"""
Code to automatize miRNA design for expression in plants based on
Zhao, J. et al(2020) Potato Virus X-Based microRNA Silencing (VbMS) In Potato.


The code follows the instructions to design short tandem target mimics (STTM)
for directed degradation of specific miRNAs.
"""
from termcolor import colored


"""
Constants for all miRNAs (all 5' - 3')
"""
# constants are designated in capital letters and won't be altered along the code
SPACER = "GTTGTTGTTGTTATGGTCTAATTTAAATATGGTCTAAAGAAGAAGAAT"
PARTIAL_5_SPACER = "GTTGTTGTTGTTATGGT"
PARTIAL_3_SPACER_REVERSE_COMPLEMENT = "ATTCTTCTTCTTTAGACCAT"


LIC1 = "CgACgACAAgACCgT" # used if the plasmid has LIC cassette
LIC2 = "gAggAgAagAgCCgT" # used if the plasmid has LIC cassette


POSITION = 10
LOOP = "UAG"


#enzyme_dict = {"SmaI":"CCCGGG", "HindIII":"AAGCTT", "BamHI":"GGATCC", "EcoRI":"GAATTC", "BsiWI":"CGTACG", "SacI":"GAGCTC", "SacII": "CCGCGG", "XbaI":"TCTAGA"}


print(colored("Welcome to STTM design program", "green"))

# definition of main function (will be called at the end of the code)
def main():
   print(colored("The following program will output the forward and reverse primers that should be used according to ", "green"))
   print (colored("Zhao, J. et al(2020) Potato Virus X-Based microRNA Silencing (VbMS) In Potato", "green"))
   print(colored("to clone the miRNA of your choice in their plasmid.", "green"))
   print("")
  
   miRNA = ask_miRNA()


   # count the number of nucleotides in the miRNA
   count = len(miRNA)
   print(colored("The number of nucleotides in your miRNA is", "yellow"),colored(count, "blue"))


   TAG_miRNA = add_bulge(miRNA, count)


   DNA_TAG_miRNA = RNA_into_DNA(TAG_miRNA)


   complement_DNA_TAG_miRNA = reverse_strand(DNA_TAG_miRNA)


   create_primers(DNA_TAG_miRNA, complement_DNA_TAG_miRNA)




"""
SECONDARY FUNCTIONS
"""
# prompt for miRNA sequence
def ask_miRNA():
   miRNA = input(colored("Please, insert miRNA sequence: ","red"))
   return miRNA


# add central budge
def add_bulge(miRNA, count):
   #Add "CTA" between 10th and 11th position (or between 11th and 12th position if it's longer)
   # depending on the miRNA lenth, the position is the 10th or 11th
  
   list1 = list(miRNA) # 1) convert string into list
   
    # check position 10 and 11 to determine which bulge to add
    if list1[10] == "G":
        if list1[11] == "C":
           list1.insert(10, LOOP) # 2) insert "CTA" between 10th and 11th nucleotides
           TAG_miRNA = ''.join(list1) # 3) convert list into string again
        elif list1[11] == "G":
           list1.insert(10, "UAC") # 2) insert "GTA" between 10th and 11th nucleotides
           TAG_miRNA = ''.join(list1) # 3) convert list into string again
        else:
           list1.insert(10, "GCGA") # 2) insert "GTGC" between 10th and 11th nucleotides
           TAG_miRNA = ''.join(list1) # 3) convert list into string again
    elif list[10] == "C":
        if list1[11] == "G"
           list1.insert(10, "UAC") # 2) insert "GTA" between 10th and 11th nucleotides
           TAG_miRNA = ''.join(list1) # 3) convert list into string again
        elif list[11] == "C":
            list1.insert(10, LOOP) # 2) insert "GTA" between 10th and 11th nucleotides
            TAG_miRNA = ''.join(list1) # 3) convert list into string again
        else:
            list1.insert(10, "GCGA") # 2) insert "GTA" between 10th and 11th nucleotides
            TAG_miRNA = ''.join(list1) # 3) convert list into string again
    else:
        list1.insert(10, "GCGA") # 2) insert "GTGC" between 10th and 11th nucleotides
        TAG_miRNA = ''.join(list1) # 3) convert list into string again
   
   return TAG_miRNA


# create reverse complement of the miRNA sequence with the budge
def RNA_into_DNA(TAG_miRNA):
   DNA_TAG_miRNA = "" # create variable to avoid UnboundError


   for nucleotide in TAG_miRNA:
       if nucleotide == 'G':
           nucleotide ='C'
       elif nucleotide == 'C':
           nucleotide = 'G'
       elif nucleotide == 'A':
           nucleotide = 'T'
       elif nucleotide == 'U' or nucleotide == 'T':  # some sequencing results provide some miRNAs with U and others with T
           nucleotide = 'A'
       else:
           print(colored("Sorry, there must have been a problem. Please, check your sequence", "red"))


       DNA_TAG_miRNA = DNA_TAG_miRNA + nucleotide      # adding complementary nt


   DNA_TAG_miRNA = DNA_TAG_miRNA[::-1]         # inverting the sequence


   return DNA_TAG_miRNA # this DNA is complementary to the RNA (will be the TM molecule)


# create reverse complement from above DNA to obtain dsDNA
def reverse_strand(DNA_TAG_miRNA):
   complement_DNA_TAG_miRNA = "" # create variable to avoid UnboundError
  
   for nucleotide in DNA_TAG_miRNA:
       if nucleotide == 'G':
           nucleotide ='C'
       elif nucleotide == 'C':
           nucleotide = 'G'
       elif nucleotide == 'A':
           nucleotide = 'T'
       elif nucleotide == 'T':
           nucleotide = 'A'
       else:
           print(colored("Sorry, there must have been a problem. Please, check your sequence", "red"))


       complement_DNA_TAG_miRNA = complement_DNA_TAG_miRNA + nucleotide  # adding complementary nt


   complement_DNA_TAG_miRNA = complement_DNA_TAG_miRNA[::-1]   # inverting the sequence


   return complement_DNA_TAG_miRNA  # this DNA is almost equal to the RNA (won't be the TM molecule)


# create primers for PCR creation of TM molecules
def create_primers(DNA_TAG_miRNA, complement_DNA_TAG_miRNA):
   # spacer sequence
   print(colored("The spacer sequence has to be","yellow"),colored(SPACER,"blue"))
   print("")
  
   check_LIC = input(colored("Does your plasmid contain LIC cassette for STTM insertion? Write yes or no: ", "red"))
   if check_LIC.lower() == "yes": # make case-insensitive and check if there's LIC cassette
       # FOR primer
       # LIC1 linker + forward miRNA seq + partial 5' spacer seq
       FOR_primer = LIC1 + DNA_TAG_miRNA + PARTIAL_5_SPACER


       # REV primer
       # LIC2 linker + reverse complement miRNA seq + partial reverse complement 3' spacer seq
       REV_primer = LIC2 + complement_DNA_TAG_miRNA + PARTIAL_3_SPACER_REVERSE_COMPLEMENT


   else:
       check_restriction = input(colored("Do you want to add any restriction sites? Write yes or no: ", "red"))
       if check_restriction.lower() == "no": # make case-insensitive and check if restriction site is not added
           # FOR primer
           # forward miRNA seq + partial 5' spacer seq
           FOR_primer = DNA_TAG_miRNA + PARTIAL_5_SPACER


           # REV primer
           # LIC2 linker + reverse complement miRNA seq + partial reverse complement 3' spacer seq
           REV_primer = complement_DNA_TAG_miRNA + PARTIAL_3_SPACER_REVERSE_COMPLEMENT


       else:
           # FOR primer
           FOR_restriction = input(colored("Write the restriction site sequence for the forward primer: ", "red"))
           # forward miRNA seq + partial 5' spacer seq
           FOR_primer = FOR_restriction.lower() + DNA_TAG_miRNA + PARTIAL_5_SPACER


           # REV primer
           REV_restriction = input(colored("Write the restriction site sequence for the reverse primer: ", "red"))
           # LIC2 linker + reverse complement miRNA seq + partial reverse complement 3' spacer seq
           REV_primer = REV_restriction.lower() + complement_DNA_TAG_miRNA + PARTIAL_3_SPACER_REVERSE_COMPLEMENT


  
   # Print primers
   print("")      # blank line
   print(colored("FOR primer: ", "yellow"), colored(FOR_primer, "blue"))
   print(colored("REV primer: ", "yellow"), colored(REV_primer, "blue"))
   print("")      # blank line


   if check_LIC.lower() == "yes":
       print(colored("Take into account that your primers are meant for LIC cassette insertion", "green"))
   elif check_restriction.lower() == "no":
       print(colored("Take into account that your primers do not have LIC or restriction enzymes", "green"))
   else:
       print(colored("Take into account that your primers are meant for insertion with restriction enzymes", "green"))


main()
