#include<stdlib.h>
#include<stdio.h>
#include<stddef.h>
#include<string.h>
#include "hashset.h"


char * CustomHashSetMagic = "CustomHashSetMagic";
char * CustomTreeSetMagic = "CustomTreeSetMagic";

#define _CUSTOM_STEP_SIZE   3

#define CUSTOM_AVL_TREE_SET_PREV( node, referenceNode )\
{\
    if( node->prev != referenceNode )\
        if(( node->prev = referenceNode )){ node->prev->parent = node; }\
}

#define CUSTOM_AVL_TREE_SET_NEXT( node, referenceNode )\
{\
    if( node->next != referenceNode )\
        if(( node->next = referenceNode )){ node->next->parent = node; }\
}


CustomSet * customSetNewHashSet( void )
{
    CustomHashSet * customSet = (CustomHashSet *)custom_malloc0( "customSetNewHashSet", sizeof(CustomHashSet) );
    if( !customSet )
    {
        return NULL;
    }

    customSet->genericSet.magic = CustomHashSetMagic;
    customSet->hashValue = customSetDefaultHashValue;
    customSet->loadFactor = 0.75;
    customSet->stepSize = _CUSTOM_STEP_SIZE;

    return (CustomSet *)customSet;
}

int customSetAdd(
CustomSet * set,  /** The set to add to               */
void * element /** Element to be added to this set */
)
{
    if( !element )
    {
        return -1;
    }

    if( CUSTOM_SET_IS_TREE_SET( set ) )
    {
        return customTreeSetAdd( (CustomTreeSet *)set, element );
    }
    return customHashSetAdd( (CustomHashSet *)set, element );
}

void customSetClear(
CustomSet * set               /** The set to clear */
)
{
    if( CUSTOM_SET_IS_TREE_SET( set ) )
    {
        customTreeSetClear( (CustomTreeSet *)set );
        return;
    }

    customHashSetClear( (CustomHashSet *)set );
}

int customSetContains(
CustomSet * set,             /** The set to use                  */
void * element            /** Element to look for             */
)
{
    if( CUSTOM_SET_IS_HASH_SET( set ))
    {
        return customHashSetContains( (CustomHashSet *)set, element );
    }

    return customTreeSetContains( (CustomTreeSet *)set, element );
}

int customSetDefaultHashValue(
const void *element     /** Element to calculate hash value for */
)
{
    return customHtHashValue( (unsigned char *)&element, sizeof(void*) );
}

int customHt_J_Zobel_Hash( const unsigned char * key, size_t keylen )
{
    int ret = 104729;

    for( ; keylen-- > 0; key++ )
    {
        ret ^= ( (ret << 5) + *key + (ret >> 2) );
    }

    return ( ret & 0x7fffffff );
}

int customHtHashValue( const unsigned char * key, size_t keylen )
{
    return ( customHt_J_Zobel_Hash( key, keylen ) & 0x7fffffff );
}


int customCollectionElementCompare(
CustomCollection * collection,    /** The collection to compare the elements for   */
void *left,
void *right
)
{
    if( left == right )
    {
        return 0;
    }

    if( !left )
    {
        if( !right )
        {
            return 0;
        }
        return -1;
    }
    if( !right )
    {
        if( !left )
        {
            return 0;
        }
        return 1;
    }

    if( collection->compare )
    {
        /*
         * There is a specific compare function for the collection
         */
        return ( *( collection->compare ) )( &left, &right );
    }

    /*
     * Use the pointers of the objects to compare as default
     */
    return customCollectionDefaultCompare( &left, &right );
}


int customCollectionDefaultCompare(
    const void *left,     /** The left element for compare  */
    const void *right     /** The right element for compare */
)
{
    char * leftPointer = *(char**)left;
    char * rightPointer = *(char**)right;

    /*
     * Use the pointers of the objects to compare
     */
    if( leftPointer < rightPointer )
    {
        return ( -1 );
    }

    if( leftPointer == rightPointer )
    {
        return 0;
    }

    return ( 1 );
}

static int customHashSetAdd(
CustomHashSet * set,            /** The set to use                              */
void * element               /** Element to be appended to this set          */
)
{
    int mask = set->capacity - 1;
    int minCapacity = set->genericSet.size + 1;
    int neededCapacity = (int)( ( (double)minCapacity ) / set->loadFactor );
    int elementIndex;
    int index;
    unsigned char ** elementPtr;

    if( !element )
    {
        return 0;
    }

    for( ;; )
    {
        if( neededCapacity > set->capacity )
        {
         
            mask = set->capacity - 1;
        }

        elementIndex = set->hashValue( element ) & mask;
        index = elementIndex;

        for( ;; )
        {
            elementPtr = &( set->pointerArray[ index ] );
            if( !*elementPtr )
            {
                /*
                 * Store the element in the table
                 */
                *elementPtr = (unsigned char *)element;

                /*
                 * The element was stored in the hash table
                 */
                set->genericSet.size++;
                set->genericSet.changeCounter++;

                return 1;
            }

            if( !customCollectionElementCompare( (CustomCollection*)set, *elementPtr,
                                              element ) )
            {
                /*
                 * Element already in hash table
                 */
                return 0;
            }

            index = ( index + set->stepSize ) & mask;
            if( elementIndex == index )
            {
                /*
                 * Linear probing brought us back to the original index,
                 * break the loop
                 */
                break;
            }
        }

        /*
         * More space is needed in the hash table
         */
        minCapacity = set->capacity + 1;
        neededCapacity = (int)( ( (double)minCapacity ) / set->loadFactor );
    }
}

static void customHashSetClear(
CustomHashSet * set                /** The set to clear */
)
{
    if( set->capacity > 0 && set->pointerArray )
    {
        memset( set->pointerArray, 0, sizeof(void*) * set->capacity );
    }
    set->genericSet.size = 0;
    set->genericSet.changeCounter++;
}

static int customTreeSetAdd(
CustomTreeSet * set,         /** The set to use                              */
void * element            /** Element to be appended to this set          */
)
{
    int size = set->genericSet.size;
    int h = 0;
    CustomTreeNode * insertResult;

    if( !element )
    {
        return 0;
    }

    insertResult = customTreeNodeInsert( set, set->rootNode, element, &h );
    if( insertResult == NULL )
    {
        return -1;
    }

    /*
     * Remember the tree after the insert
     */
    insertResult->parent = NULL;
    set->rootNode = insertResult;

    if( size == set->genericSet.size )
    {
        return 0;
    }

    return 1;
}

static void customTreeSetClear(
CustomTreeSet * set                /** The set to clear */
)
{
    if( set->rootNode )
    {
        customTreeNodeFree( set->rootNode );
        set->rootNode = NULL;
    }

    set->genericSet.size = 0;
    set->genericSet.changeCounter++;
}

static int customTreeSetContains(
CustomTreeSet * set,         /** The set to use                  */
void * element            /** Element to look for             */
)
{
    CustomTreeNode * node = set->rootNode;
    int compareResult;

    if( set->genericSet.size == 0 )
    {
        return 0;
    }

    while( node )
    {
        compareResult = customCollectionElementCompare( (CustomCollection*)set, element,
                                                     node->element );
        if( !compareResult )
        {
            return 1;
        }

        if( compareResult < 0 )
        {
            node = node->prev;
        }
        else
        {
            node = node->next;
        }
    }
    return 0;
}

static int customHashSetContains(
CustomHashSet * set,         /** The set to use                  */
void * element            /** Element to look for             */
)
{
    int mask = set->capacity - 1;
    int index;
    void * pointer;

    if( set->genericSet.size == 0 || set->capacity < 1 )
    {
        return 0;
    }
    index = set->hashValue( element ) & mask;

    for( ;; )
    {
        pointer = set->pointerArray[ index ];
        if( !pointer )
        {
            return 0;
        }

        if( !customCollectionElementCompare( (CustomCollection*)set, pointer, element ) )
        {
            /*
             * Element in hash table
             */
            return 1;
        }

        index = ( index + set->stepSize ) & mask;
    }
}

void * custom_malloc0(
char   * tag,        /** tag used for memory leak detection */
size_t   size        /** number of bytes to allocate        */
)
{
    void * ptr = malloc( size );
    if( !ptr )
    {
        return( 0 );
    }

    memset( ptr, 0, size );


    return( ptr );
}

static CustomTreeNode * customTreeNodeInsert(
CustomTreeSet * set,                      /** The AVL tree to insert into    */
CustomTreeNode * parentNode,              /** The parent node to insert to   */
void * element,                        /** The element to insert          */
int  * heightChanged                   /** Set if the tree height changed */
)
{
    CustomTreeNode * p1;
    CustomTreeNode * p2;
    int compareResult;

    *heightChanged = 0;

    if( parentNode == NULL )
    {
        /*
         * Element is not in the tree yet, insert it.
         */
        *heightChanged = 1;
        return customTreeNodeCreate( set, element );
    }

    compareResult = customCollectionElementCompare( (CustomCollection*)set, element,
                                                 parentNode->element );
    if( !compareResult )
    {
        /*
         * Element already in tree
         */
        return parentNode;
    }

    if( compareResult < 0 )
    {
        /*
         * Insert into left sub tree
         */
        p1 = customTreeNodeInsert( set, parentNode->prev, element, heightChanged );
        if( !p1 )
        {
            return p1;
        }
        CUSTOM_AVL_TREE_SET_PREV( parentNode, p1 );

        if( !*heightChanged )
        {
            return parentNode;
        }

        /*
         * Left sub tree increased in height
         */
        if( parentNode->balance == 1 )
        {
            parentNode->balance = 0;
            *heightChanged = 0;
            return parentNode;
        }

        if( parentNode->balance == 0 )
        {
            parentNode->balance = -1;
            return parentNode;
        }

        /*
         * Balancing needed
         */
        p1 = parentNode->prev;

        if( p1->balance == -1 )
        {
            /*
             * Simple LL rotation
             */
            CUSTOM_AVL_TREE_SET_PREV( parentNode, p1->next );

            CUSTOM_AVL_TREE_SET_NEXT( p1, parentNode );
            parentNode->balance = 0;

            parentNode = p1;
            parentNode->balance = 0;
            *heightChanged = 0;
            return parentNode;
        }

        /*
         * double LR rotation
         */
        p2 = p1->next;

        CUSTOM_AVL_TREE_SET_NEXT( p1, p2->prev );

        CUSTOM_AVL_TREE_SET_PREV( p2, p1 );

        CUSTOM_AVL_TREE_SET_PREV( parentNode, p2->next );

        CUSTOM_AVL_TREE_SET_NEXT( p2, parentNode );

        if( p2->balance == -1 )
        {
            parentNode->balance = 1;
        }
        else
        {
            parentNode->balance = 0;
        }

        if( p2->balance == 1 )
        {
            p1->balance = -1;
        }
        else
        {
            p1->balance = 0;
        }
        parentNode = p2;
        parentNode->balance = 0;
        *heightChanged = 0;
        return parentNode;
    }

    /*
     * Insert into right sub tree
     */
    p1 = customTreeNodeInsert( set, parentNode->next, element, heightChanged );
    if( !p1 )
    {
        return p1;
    }
    CUSTOM_AVL_TREE_SET_NEXT( parentNode, p1 );

    if( !*heightChanged )
    {
        return parentNode;
    }

    /*
     * Right sub tree increased in height
     */
    if( parentNode->balance == -1 )
    {
        parentNode->balance = 0;
        *heightChanged = 0;
        return parentNode;
    }

    if( parentNode->balance == 0 )
    {
        parentNode->balance = 1;
        return parentNode;
    }

    /*
     * Balancing needed
     */
    p1 = parentNode->next;

    if( p1->balance == 1 )
    {
        /*
         * Simple RR rotation
         */
        CUSTOM_AVL_TREE_SET_NEXT( parentNode, p1->prev );

        CUSTOM_AVL_TREE_SET_PREV( p1, parentNode );
        parentNode->balance = 0;

        parentNode = p1;
        parentNode->balance = 0;
        *heightChanged = 0;
        return parentNode;
    }

    /*
     * double RL rotation
     */
    p2 = p1->prev;

    CUSTOM_AVL_TREE_SET_PREV( p1, p2->next );

    CUSTOM_AVL_TREE_SET_NEXT( p2, p1 );

    CUSTOM_AVL_TREE_SET_NEXT( parentNode, p2->prev );

    CUSTOM_AVL_TREE_SET_PREV( p2, parentNode );

    if( p2->balance == 1 )
    {
        parentNode->balance = -1;
    }
    else
    {
        parentNode->balance = 0;
    }

    if( p2->balance == -1 )
    {
        p1->balance = 1;
    }
    else
    {
        p1->balance = 0;
    }
    parentNode = p2;
    parentNode->balance = 0;
    *heightChanged = 0;
    return parentNode;
}

static CustomTreeNode * customTreeNodeCreate(
CustomTreeSet * set,
void * element
)
{
    CustomTreeNode * newNode = (CustomTreeNode *)custom_malloc0( "customTreeNodeCreate", sizeof(CustomTreeNode) );
    if( !newNode )
    {
        return newNode;
    }

    newNode->element = element;

    set->genericSet.size++;
    set->genericSet.changeCounter++;

    return newNode;
}

static void customTreeNodeFree(
CustomTreeNode * node                /** The node to free */
)
{
    if( node->prev )
    {
        customTreeNodeFree( node->prev);
    }

    if( node->next )
    {
        customTreeNodeFree( node->next);
    }

    CUSTOM_FREE( node );
}

