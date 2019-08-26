// Matthew Gomex
// COP 3223
// 8/20/2019
// Generates a strand of DNA based on user input

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

#define NUCLEOTIDES 4

void DNAStrndGen(char bases[], int size);       // This is responsible for the generation and analysis of the DNA sequence

int main() 

{
    srand(time(NULL));      // Seeds the random number generator so that the program will produce different sequences for each run

    char bases[] = {"AGCT"};
    int size, i;

    printf("How long do you want your DNA sequence to be?\n");

    scanf("%d", &size);     // Collects the size from the user so that it can be used across most of the functions below

    DNAStrndGen(bases, size);

    printf("\n");

}

void DNAStrndGen(char bases[], int size)

{
    int i; 
    char DNAStrnd[size];

    printf("Randomly Generated DNA Sequence: ");

    for(i=0; i<size; i++){

            DNAStrnd[i] = bases[rand() % 4];        // This produces the DNA sequence by randomly choosing which nucleotide to produce,
                                                    // zero is included in the generation to ensure that 'A' can be produced.
            printf("%c", DNAStrnd[i]);

        }
    
    printf("\n");

    int baseCount[4] = {0};             // Sets the nucleotide counter to zero so that it can count effectively

        for(i=0; i<size; i++){          // The amount of times the for loop runs is equal to the size of the sequence / user input size
                                        // This is done so that each base is counted. Using an array was the best way to count the nucleotides
            if(DNAStrnd[i] == 'A'){     // because it allows for all the data to be stored in one variable rather than dealing with multiple.

                baseCount[0] = baseCount[0] + 1;

            }

            if(DNAStrnd[i] == 'G'){

                baseCount[1] = baseCount[1] + 1;

            }
        
            if(DNAStrnd[i] == 'C'){

                baseCount[2] = baseCount[2] + 1;

            }

            if(DNAStrnd[i] == 'T'){

                baseCount[3] = baseCount[3] + 1;

            }
        }

    printf("Adenine: %d\n", baseCount[0]);

    printf("Guanine: %d\n", baseCount[1]);

    printf("Cytosine: %d\n", baseCount[2]);

    printf("Thymine: %d\n", baseCount[3]);

        char DNAStrnd2[size];       // Here the DNA sequence is copied except this time each base now is accompanied by its complement
                                    // using a series of basic if statements. 
    for(i=0; i<size; i++){

        if(DNAStrnd[i] == 'A'){

            DNAStrnd2[i] = 'T';

        }

        if(DNAStrnd[i] == 'T'){

            DNAStrnd2[i] = 'A';

        }

        if(DNAStrnd[i] == 'C'){

            DNAStrnd2[i] = 'G';

        }

        if(DNAStrnd[i] == 'G'){

            DNAStrnd2[i] = 'C';

        }
        
    }

    int start, end = 0;             // In order to reverse the complementary sequence, the new reverse complementary strand is declared
                                    // with the same size as DNAStrnd2 (the 2nd sequence that was produced) then each element/base of the reverse comp.            
    end = size - 1;                 // is set equal to the last element/base of DNAStrnd2. The counter is reduced by one to ensure the loop continues 
    char CompDNA[size];             // until the final element of the reverse comp. strand is equal to the first base of DNAStrnd2. 

        for(start = 0; start<size; start++){

            CompDNA[start] = DNAStrnd2[end];
            end--;
        }
    
    CompDNA[start] = '\0';
    
    printf("Reverse Complementary Strand: %s\n", CompDNA);

    printf("mRNA Strand: ");

    char mRNA[size];                // The mRNA strand is created by simply replacing the thymine base, 'T', with uracil 'U'. The for 
                                    // loop verifies that T that is in the CompDNA array is replaced with a U when creating the mRNA sequence
        for(i=0; i<size; i++){

            if(CompDNA[i] == 'T'){

                mRNA[i] = 'U';
            }

            if(CompDNA[i] == 'A'){

                mRNA[i] = 'A';
            }

            if(CompDNA[i] == 'C'){

                mRNA[i] = 'C';
            }

            if(CompDNA[i] == 'G'){

                mRNA[i] = 'G';
            }

        printf("%c", mRNA[i]);

    }
        printf("\n");

        int PotAminoAcidCount; 
        int k;
        char startCodon[] = {"AUG"};
        char stopCodon[] = {"UAG"};
        struct AminoAcid {

            char type[10];
            char Codon[3];
            int freq = 0;

        };

        struct AminoAcid Histidine;
            strcpy( Histidine.type, "Essential");
            strcpy( Histidine.Codon, "CAC");
        struct AminoAcid Isoleucine;
            strcpy( Isoleucine.type, "Essential");
            strcpy( Isoleucine.Codon, "AUU");
        struct AminoAcid Leucine;
            strcpy( Leucine.type, "Essential");
            strcpy( Leucine.Codon, "CUA");
        struct AminoAcid Lysine;
            strcpy( Lysine.type, "Essential");
            strcpy( Lysine.Codon, "AAG");
        struct AminoAcid Methionine;
            strcpy( Methionine.type, "Essential");
        struct AminoAcid Phenylalanine;
            strcpy( Phenylalanine.type, "Essential");
            strcpy( Phenylalanine.Codon, "UUC");
        struct AminoAcid Threonine;
            strcpy( Threonine.type, "Essential");
            strcpy( Threonine.Codon, "ACA");
        struct AminoAcid Tryptophan;
            strcpy( Tryptophan.type, "Essential");
            strcpy( Tryptophan.Codon, "UAU");
        struct AminoAcid Valine;
            strcpy( Valine.type, "Essential");
            strcpy( Valine.Codon, "GUA");

        PotAminoAcidCount = (size/3); 

        printf("%d potential amino acids can be produced\n", PotAminoAcidCount);

        printf("\n");

            for(i=0; i<size; i+=3){

                        if(mRNA[i] == startCodon[0] && mRNA[i+1] == startCodon[1] && mRNA[i+2] == startCodon[2]){

                            printf("Amino Acid Produced: Methionine\n");
                            Methionine.freq +=1;
                        
                        }
                        if(mRNA[i] == Histidine.Codon[0] && mRNA[i+1] == Histidine.Codon[1] && mRNA[i+2] == Histidine.Codon[2]){

                            printf("Amino Acid Produced: Histidine\n");
                            Histidine.freq +=1;
                        }
                        if(mRNA[i] == Isoleucine.Codon[0] && mRNA[i+1] == Isoleucine.Codon[1] && mRNA[i+2] == Isoleucine.Codon[2]){

                            printf("Amino Acid Produced: Isoleucine\n");
                            Isoleucine.freq +=1;
                        }
                        if(mRNA[i] == Leucine.Codon[0] && mRNA[i+1] == Leucine.Codon[1] && mRNA[i+2] == Leucine.Codon[2]){

                            printf("Amino Acid Produced: Leucine\n");
                            Leucine.freq +=1;
                        }
                        if(mRNA[i] == Phenylalanine.Codon[0] && mRNA[i+1] == Phenylalanine.Codon[1] && mRNA[i+2] == Phenylalanine.Codon[2]){

                            printf("Amino Acid Produced: Phenylalanine\n");
                            Phenylalanine.freq +=1;
                        }
                        if(mRNA[i] == Threonine.Codon[0] && mRNA[i+1] == Threonine.Codon[1] && mRNA[i+2] == Threonine.Codon[2]){

                            printf("Amino Acid Produced: Threonine\n");
                            Threonine.freq +=1;
                        }
                        if(mRNA[i] == Tryptophan.Codon[0] && mRNA[i+1] == Tryptophan.Codon[1] && mRNA[i+2] == Tryptophan.Codon[2]){

                            printf("Amino Acid Produced: Tryptophan\n");
                            Tryptophan.freq +=1;
                        }
                        if(mRNA[i] == Valine.Codon[0] && mRNA[i+1] == Valine.Codon[1] && mRNA[i+2] == Valine.Codon[2]){

                            printf("Amino Acid Produced: Valine\n");
                            Valine.freq +=1;
                        }
            }
        
        printf("Methionine Count: %d\n", Methionine.freq);
        printf("Histidine Count: %d\n", Histidine.freq);
        printf("Isoleucine Count: %d\n", Isoleucine.freq); 
        printf("Leucine Count: %d\n", Leucine.freq);
        printf("Phenylalanine Count: %d\n", Phenylalanine.freq);
        printf("Threonine Count: %d\n", Threonine.freq);
        printf("Tryptophan Count: %d\n", Tryptophan.freq);   
        printf("Valine Count: %d\n", Valine.freq);
}








    








    

















