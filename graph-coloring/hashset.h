#ifndef _HASHSET_H_
#define _HASHSET_H_

#include <stddef.h>
#include <stdio.h>

#define CUSTOM_ERROR_BASE                    1000
#define CUSTOM_ERROR_PARAM_ELEMENT           ( CUSTOM_ERROR_BASE + 13 )

#define CUSTOM_FREE( ptr ) if( ptr ){ free( ptr ); ptr = 0; }

typedef struct CustomCollection
{
    char * magic;         /* The magic string of collections                 */
    int size;             /* The size of the collections                     */

    /* A user defined element compare function                               */
    int (*compare)( const void * prev, const void * next );

    unsigned long changeCounter; /* Number of changes on the collections     */

} CustomSet;

typedef CustomSet CustomCollection;


typedef struct CustomHashSet_s
{
    CustomSet genericSet;     /* The generic set definition of the hash set     */

    /* The array of all pointers hashed                                      */
    unsigned char ** pointerArray;

    int capacity;          /* The capacity of the pointer array              */
    int stepSize;          /* The step size used for collision resolution    */
    double loadFactor;     /* The load factor of the hash set                */

    /* A user defined element hash value function                            */
    int (*hashValue)( const void * element );

} CustomHashSet;

typedef struct CustomTreeNode_s
{
    void * element;                     /* The element the node points to  */

    struct CustomTreeNode_s * prev;        /* The left node                   */
    struct CustomTreeNode_s * next;        /* The right node                  */

    struct CustomTreeNode_s * parent;      /* The parent node                 */
    int balance;                        /* AVL balance information of node */


} CustomTreeNode;


typedef struct CustomTreeSet_s
{
    CustomSet genericSet;       /* The generic set definition of the tree set   */

    CustomTreeNode * rootNode;  /* The root node of the AVL tree                */

} CustomTreeSet;


#define CUSTOM_SET_IS_HASH_SET( SET )\
    (SET ? (((CustomSet*)SET)->magic == CustomHashSetMagic) : 0 )
	
#define CUSTOM_SET_IS_TREE_SET( SET )\
    (SET ? (((CustomSet*)SET)->magic == CustomTreeSetMagic) : 0 )

extern int customSetContains(
    CustomSet * set,             /** The set to use                  */
    void * element            /** Element to look for             */
    );
	
extern void customSetClear(
        CustomSet * set     /** The set to clear */
        );
		
extern int customSetAdd(
        CustomSet * set,   /** The set to use                                */
        void * element  /** Element to be appended to this set            */
        );

extern CustomSet * customSetNewHashSet( void );

extern int customHtHashValue( const unsigned char * key, size_t keylen );

extern int customSetDefaultHashValue(
    const void *element     /** Element to calculate hash value for */
);

extern int    custom_errno;

static int customTreeSetContains(
CustomTreeSet * set,         /** The set to use                  */
void * element            /** Element to look for             */
);

void * custom_malloc0(
char   * tag,        /** tag used for memory leak detection */
size_t   size        /** number of bytes to allocate        */
);

static int customHashSetAdd(
CustomHashSet * set,            /** The set to use                              */
void * element               /** Element to be appended to this set          */
);

int customCollectionDefaultCompare(
    const void *left,     /** The left element for compare  */
    const void *right     /** The right element for compare */
);

static void customHashSetClear(
CustomHashSet * set                /** The set to clear */
);

static int customTreeSetAdd(
CustomTreeSet * set,         /** The set to use                              */
void * element            /** Element to be appended to this set          */
);

static void customTreeSetClear(
CustomTreeSet * set                /** The set to clear */
);

static int customHashSetContains(
CustomHashSet * set,         /** The set to use                  */
void * element            /** Element to look for             */
);

int customHtHashValue( const unsigned char * key, size_t keylen );

static void customTreeNodeFree(
CustomTreeNode * node                /** The node to free */
);

static CustomTreeNode * customTreeNodeCreate(
CustomTreeSet * set,
void * element
);

static CustomTreeNode * customTreeNodeInsert(
CustomTreeSet * set,                      /** The AVL tree to insert into    */
CustomTreeNode * parentNode,              /** The parent node to insert to   */
void * element,                        /** The element to insert          */
int  * heightChanged                   /** Set if the tree height changed */
);



#endif