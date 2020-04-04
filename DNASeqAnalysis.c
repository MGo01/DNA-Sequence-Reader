// Matthew Gomex
// COP 3223
// 8/20/2019
// Generates a strand of DNA based on pre-existing files

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#define MAX_LENGTH 1024
#define AMNT_NUCLEOTIDES 4
#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

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
} LinkedList;

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
  newList->nucleotideCount = malloc(sizeof(int) * AMNT_NUCLEOTIDES);

  for (i = 0; i < AMNT_NUCLEOTIDES; i++)
    newList->nucleotideCount[i] = 0;

  return newList;
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

LinkedList *copyStrand(char *fileName)
{
  LinkedList *list = NULL;
  FILE *ifp = NULL;
  char buffer[MAX_LENGTH];

  if ((ifp = fopen(fileName, "r")) == NULL)
  {
    fprintf(stderr, "Error, file: %s cannot be opened\n", fileName);
    return NULL;
  }

  while (fscanf(ifp, "%s", buffer) != EOF)
  {
    list = createList();
    insertNodes(list, buffer);
  }

  return list;
}

void printDNA(LinkedList *list)
{
  node *temp = NULL;

  if (list == NULL || list->head == NULL)
  {
    printf("Empty List\n");
    return;
  }

  temp = list->head;

  while (temp != NULL)
  {
    printf("%c -> %c", temp->data, (temp->next == NULL) ? '\n' : ' ');
    temp = temp->next;
  }

  printf("Adenine Count: %d\nCytosine Count: %d\nThymine Count: %d\nGuanine Count: %d\n",
          list->nucleotideCount[ADENINE], list->nucleotideCount[CYTOSINE], list->nucleotideCount[THYMINE],
          list->nucleotideCount[GUANINE]);
}

int main(void)
{
  LinkedList *dnaStrand = copyStrand("DNA.txt");
  printDNA(dnaStrand);

  return 0;
}
