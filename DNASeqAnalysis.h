#ifndef __DNASEQANALYSIS_H
#define __DNASEQANALYSIS_H

#include <limits.h>

#define EMP_QUEUE_ERR INT_MIN
#define MAX_LENGTH 1024
#define AMNT_NUCLEOTIDES 4
#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

typedef enum StrandType { DNA, COMPLEMENTARY, RNA } StrandType;
typedef enum ChargeType { POSITIVE, NEGATIVE, NEUTRAL } ChargeType;
typedef enum PolarityType { POLAR, NONPOLAR } PolarityType;

const char StopCodon[4] = "UAG\0";
const char Methionine[4] = "AUG\0";
const char Valine[4] = "GUC\0";
const char Histidine[4] = "CAC\0";
const char Leucine[4] = "CUG\0";
const char Threonine[4] = "ACA\0";
const char Proline[4] = "CCG\0";
const char Glutamine[4] = "GAG\0";
const char HealthyHemo[8] = "MVHLTPQ\0";

typedef struct node
{
  char data;
  struct node *next;
} node;

typedef struct LinkedList
{
  node *head;
  node *tail;
  int *nucleotideCount;
  StrandType strandType;
} LinkedList;

typedef struct AminoAcid
{
  char *codon;
  char code;
  ChargeType charge;
  PolarityType polarity;
  struct AminoAcid *next;
} AminoAcid;

typedef struct Queue
{
  char *array;
  int front;
  int back;
  int capacity;
  int size;
} Queue;

typedef struct Protein
{
  Queue *sequence;
  char *name;
} Protein;

#endif
