// Matthew Gomex
// COP 3223
// 8/20/2019
// Generates a strand of DNA based on pre-existing files

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "DNASeqAnalysis.h"

void printQueue(Queue *q);

node *createNode(char data)
{
  node *newNode = malloc(sizeof(node));

  newNode->data = data;
  newNode->next = NULL;

  return newNode;
}

LinkedList *createList(void)
{
  int i;

  LinkedList *newList = malloc(sizeof(LinkedList));

  newList->head = newList->tail = NULL;
  newList->strandType = DNA;
  newList->nucleotideCount = malloc(sizeof(int) * AMNT_NUCLEOTIDES);

  for (i = 0; i < AMNT_NUCLEOTIDES; i++)
    newList->nucleotideCount[i] = 0;

  return newList;
}

Queue *createQueue(int capacity)
{
  Queue *new = malloc(sizeof(Queue));
  new->array = malloc(sizeof(char) * capacity);
  new->capacity = capacity - 1;
  new->array[capacity] = '\0';
  new->front = 0;
  new->back = 0;
  new->size = 0;

  return new;
}

AminoAcid *createAminoAcid(char *codon, ChargeType charge, PolarityType polarity)
{
  AminoAcid *newAmino = malloc(sizeof(AminoAcid));

  newAmino->codon = codon;
  newAmino->charge = charge;
  newAmino->polarity = polarity;
  newAmino->next = NULL;

  return newAmino;
}

Protein *createProtein(char *name, int length)
{
  Protein *newProtein = malloc(sizeof(Protein));

  newProtein->sequence = createQueue(length);
  newProtein->name = name;

  return newProtein;
}

int isQueueFull(Queue *q)
{
  return (q == NULL || q->size == q->capacity);
}

int isQueueEmpty(Queue *q)
{
  return (q == NULL || q->size == 0);
}

void enqueue(Queue *q, char data)
{
  if (isQueueFull(q))
  {
    fprintf(stderr, "Error, Queue is Full in enqueue().\n");
    return;
  }

  else
  {
    q->array[q->back] = data;
    q->size++;
    q->back = q->size;
    return;
  }
}

int dequeue(Queue *q)
{
  char storeData;
  node *temp = NULL;

  if (isQueueEmpty(q))
    return EMP_QUEUE_ERR;

  storeData = q->array[q->front];
  q->front++;
  q->size--;

  return storeData;
}

void resetQueue(Queue *q)
{
  char throwOut;

  if (isQueueEmpty(q))
    return;

  else
  {
    while (!isQueueEmpty(q))
      throwOut = dequeue(q);

    q->front = q->back = 0;

    return;
  }
}

Queue *destroy_queue(Queue *q)
{
  if (q == NULL)
    return NULL;

  free(q->array);
  free(q);

  return NULL;
}

LinkedList *tailInsert(LinkedList *list, char data)
{
  if (list == NULL)
  {
    return NULL;
  }

  if (list->head == NULL)
  {
    list->head = list->tail = createNode(data);
    return list;
  }

  else
  {
    list->tail->next = createNode(data);
    list->tail = list->tail->next;
    return list;
  }
}

void insertNodes(LinkedList *list, char *str)
{
  int i, n = strlen(str);

  if (list == NULL)
  {
    fprintf(stderr, "Error, empty list in insertNodes().\n");
    return;
  }

  for (i = 0; i < n; i++)
  {
    if (str[i] == 'A')
      list->nucleotideCount[ADENINE]++;

    if (str[i] == 'C')
      list->nucleotideCount[CYTOSINE]++;

    if (str[i] == 'G')
      list->nucleotideCount[GUANINE]++;

    if (str[i] == 'T')
      list->nucleotideCount[THYMINE]++;

    list = tailInsert(list, str[i]);
  }

  return;
}

LinkedList *makeRNAStrand(LinkedList *list)
{
  node *temp = NULL;

  if (list == NULL || list->head == NULL)
  {
    fprintf(stderr, "Empty list in makeRNAstrand().\n");
    return NULL;
  }

  temp = list->head;

  while (temp != NULL)
  {
    if (temp->data == 'T')
    {
      temp->data = 'U';
      temp = temp->next;
    }
    else
      temp = temp->next;
  }

  list->strandType = RNA;

  return list;
}

LinkedList *makeComplementaryStrand(LinkedList *list)
{
  node *temp = NULL;

  if (list == NULL || list->head == NULL)
  {
    fprintf(stderr, "Empty list in makeComplementaryStrand().\n");
    return NULL;
  }

  temp = list->head;

  while (temp != NULL)
  {
    if (temp->data == 'A')
    {
      temp->data = 'T';
      temp = temp->next;
    }

    else if (temp->data == 'T')
    {
      temp->data = 'A';
      temp = temp->next;
    }

    else if (temp->data == 'C')
    {
      temp->data = 'G';
      temp = temp->next;
    }

    else if (temp->data == 'G')
    {
      temp->data = 'C';
      temp = temp->next;
    }
  }

  list->strandType = COMPLEMENTARY;

  return list;
}

LinkedList *copyStrand(char *fileName)
{
  LinkedList *list = NULL;
  FILE *ifp = NULL;
  char buffer[MAX_LENGTH];

  if ((ifp = fopen(fileName, "r")) == NULL)
  {
    fprintf(stderr, "Error, file: /%s/ cannot be opened\n", fileName);
    return NULL;
  }

  while (fscanf(ifp, "%s", buffer) != EOF)
  {
    list = createList();
    insertNodes(list, buffer);
  }

  return list;
}

void printStrand(LinkedList *list)
{
  node *temp = NULL;

  if (list == NULL || list->head == NULL)
  {
    printf("Empty list in printStrand().\n");
    return;
  }

  temp = list->head;

  while (temp != NULL)
  {
    printf("%c -> %c", temp->data, (temp->next == NULL) ? '\n' : ' ');
    temp = temp->next;
  }

  printf("%s Stats:\n\tAdenine Count: %d\n\tCytosine Count: %d\n\t%s Count: %d\n\tGuanine Count: %d\n\n",
         (list->strandType == COMPLEMENTARY) ? "Complementary DNA" : "DNA/RNA",
         list->nucleotideCount[ADENINE], list->nucleotideCount[CYTOSINE],
         (list->strandType == RNA) ? "Uracil" : "Thymine",
         list->nucleotideCount[THYMINE],
         list->nucleotideCount[GUANINE]);
}

Protein *createSequence(char *name, LinkedList *strand)
{
  Queue *q = NULL;
  node *temp = NULL;
  Protein *p = NULL;
  int i = 0, flag = 0, length = 0, counter = 0;

  if (strand == NULL || strand->head == NULL)
  {
    fprintf(stderr, "Empty strand in createSequence().\n");
    return NULL;
  }

  for (i = 0; i < 4; i++)
    length += strand->nucleotideCount[i];

  length /= 3;
  p = createProtein(name, length + 2);
  q = createQueue(4);

  temp = strand->head;

  while (temp != NULL)
  {
    enqueue(q, temp->data);
    counter++;
    temp = temp->next;

    if (counter % 3 == 0)
    {
      if (flag == 0 && (strcmp(q->array, Methionine) == 0))
      {
        printf("Methionine has been produced!\n");
        enqueue(p->sequence, 'M');
        resetQueue(q);
        flag = 1;
        continue;
      }

      if (flag && (strcmp(q->array, Valine) == 0))
      {
        printf("Valine has been produced!\n");
        enqueue(p->sequence, 'V');
        resetQueue(q);
        continue;
      }

      if (flag && (strcmp(q->array, Histidine) == 0))
      {
        printf("Histidine has been produced!\n");
        enqueue(p->sequence, 'H');
        resetQueue(q);
        continue;
      }

      if (flag && (strcmp(q->array, Leucine) == 0))
      {
        printf("Leucine has been produced!\n");
        enqueue(p->sequence, 'L');
        resetQueue(q);
        continue;
      }

      if (flag && (strcmp(q->array, Threonine) == 0))
      {
        printf("Threonine has been produced!\n");
        enqueue(p->sequence, 'T');
        resetQueue(q);
        continue;
      }

      if (flag && (strcmp(q->array, Proline) == 0))
      {
        printf("Proline has been produced!\n");
        enqueue(p->sequence, 'P');
        resetQueue(q);
        continue;
      }

      if (flag && (strcmp(q->array, Glutamine) == 0))
      {
        printf("Glutamine has been produced!\n");
        enqueue(p->sequence, 'Q');
        resetQueue(q);
        continue;
      }

      else
      {
        printf("Nothing was been produced!\n");
        resetQueue(q);
      }
    }
  }

  return p;
}

void checkForGeneticDisorders(Protein *p)
{
  if (p == NULL || p->sequence == NULL)
  {
    fprintf(stderr, "Error, empty protein in SickleCheck().\n");
    return;
  }

  printf("\nGenetic Disorders Results:\n\t\nProtein: %s\n", p->name);

  if (strcmp(p->sequence->array, HealthyHemo) == 0)
    printf("\n\tSickle Cell [N]\n");
  else
    printf("\n\tSickle Cell [Y]\n");

  return;

}

void printQueue(Queue *q)
{
  int i = 0;

  if (isQueueEmpty(q))
    return;

  else
  {
    printf("\n");

    for (i = 0; i < q->capacity; i++)
      printf("Q[%d]: %c ", i, q->array[i]);

    printf("\n");
    return;
  }
}

int main(void)
{
  LinkedList *userList = NULL;
  Protein *userProtein = NULL;
  char buffer[1024];
  int choice, loop = 1;

  printf("\nDNA Sequence Analysis Tool\n");
  printf("\n--------------------------------\n");
  printf("1 - Copy DNA Strand from file\n");
  printf("2 - Convert DNA Strand to Complementary DNA\n");
  printf("3 - Simulate Translation on DNA Strand\n");
  printf("4 - Simulate Transcription on RNA Strand\n");
  printf("5 - Print Genetic Disorders Report\n");
  printf("6 - Print current Genetic Strand\n");
  printf("7 - Print current protein\n");
  printf("8 - Free current DNA strand memory\n");
  printf("9 - Quit Program\n\n");

  scanf("%d", &choice);

  if (choice == 1)
  {
    printf("Please enter the name of the file you would like to copy\n");
    scanf("%s", buffer);
  }

  if (choice == 4)
  {
    printf("Please enter the name of the protein you would like to search for\n");
    scanf("%s", buffer);
  }

  while (choice != 9)
  {
    switch (choice) {

      case 1:
      {
        userList = copyStrand(buffer);
        loop = 1;
        break;
      }

      case 2:
      {
        userList = makeComplementaryStrand(userList);
        loop = 1;
        break;
      }

      case 3:
      {
        userList = makeRNAStrand(userList);
        loop = 1;
        break;
      }

      case 4:
      {
        userProtein = createSequence(buffer, userList);
        loop = 1;
        break;
      }

      case 5:
      {
        checkForGeneticDisorders(userProtein);
        loop = 1;
        break;
      }

      case 6:
      {
        printStrand(userList);
        loop = 1;
        break;
      }

      case 7:
      {
        printQueue(userProtein->sequence);
        loop = 1;
        break;
      }

      case 8:
      {
        printf("UNDER CONSTRUCTION\n");
        loop = 1;
        break;
      }

      case 9:
      {
        printf("Thank you for using DNA Sequence Analysis Tool!\n");
        exit(-1);
        break;
      }

      default:
      {
        printf("Invalid selection, please enter an integer between 1 and 9\n");
        loop = 1;
        break;
      }
    }

    printf("\n1 - Copy DNA Strand from file\n");
    printf("2 - Convert DNA Strand to Complementary DNA\n");
    printf("3 - Simulate Translation on DNA Strand\n");
    printf("4 - Simulate Transcription on RNA Strand\n");
    printf("5 - Print Genetic Disorders Report\n");
    printf("6 - Print current Genetic Strand\n");
    printf("7 - Print current protein\n");
    printf("8 - Free current DNA strand memory\n");
    printf("9 - Quit Program\n");
    printf("Note: Remember to always free memory before opening another file\n");

    scanf("%d", &choice);
    printf("\n\n");
  }

  return 0;
}
