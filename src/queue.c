/***************************************************************************
 *   Copyright (C) 2006 by Piotr J. Durka Dobieslaw Ircha, Rafal Kus, Marek Matysiak   *
 *   durka@fuw.edu.pl, rircha@fuw.edu.pl, rkus@fuw.edu.pl				     	*
 *   Department of Biomedical Physics at Warsaw University			     		*
 *   http://brain.fuw.edu.pl, http://eeg.pl						     		*
 *												     		*
 *   This program is free software; you can redistribute it and/or modify	     		*
 *   it under the terms of the GNU General Public License as published by	     		*
 *   the Free Software Foundation; either version 2 of the License, or 		     	*
 *   (at your option) any later version.							     		*
 *												     		*
 *   This program is distributed in the hope that it will be useful,		     		*
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of	     	*
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 		*
 *   GNU General Public License for more details.					     		*
 *												     		*
 *   You should have received a copy of the GNU General Public License		     	*
 *   along with this program; if not, write to the					     		*
 *   Free Software Foundation, Inc.,							     		*
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.			     	*
 ***************************************************************************/

#include<stdlib.h>
#include"include/queue.h"
#include"include/types.h"

QueueNode* createNode(void)
{
    QueueNode *queueNode;
    queueNode = (QueueNode *)malloc(sizeof(QueueNode));
    queueNode->nextNode = NULL;
    return queueNode;
}

void freeQueueNode(QueueNode *queueNode, void (*funPtr)(void *data))
{
    if(funPtr == NULL)
	free((void *)(queueNode->data));
    else
	(*funPtr)(queueNode->data);

    free(queueNode);
}
	
Queue* createQueue(void)
{
    Queue *queue = (Queue *)malloc(sizeof(Queue));
    queue->size = 0;
    queue->firstNode = queue->lastNode = NULL;

    return queue;
}
	
Queue* addNode(Queue *queue, void *data)
{
    if(!queue->firstNode)
    {
	queue->firstNode = queue->lastNode = createNode();
	queue->firstNode->data = data;
	queue->size = queue->size + 1;
    }
    else
    {
	queue->lastNode->nextNode = createNode();
	queue->lastNode = queue->lastNode->nextNode;
	queue->lastNode->data = data;
	queue->size = queue->size + 1;
    }
    return queue;
}

void* readFirstNode(Queue *queue)
{
    return queue->firstNode->data;
}

void* readLastNode(Queue *queue)
{
    return queue->lastNode->data;
}

/* return data at n-th node, nodes are numerated from 0 */

void* readNNode(Queue *queue, unsigned int nodeNumber)
{
    unsigned int nodeCounter;
    QueueNode *node = queue->firstNode;

    for(nodeCounter = 0;nodeCounter<nodeNumber;nodeCounter++)
	node = node->nextNode;

    return node->data;
}
	
void* substractNode(Queue *queue)
{
    QueueNode *node;
    void *data;
	
    node = queue->firstNode;
    queue->firstNode = queue->firstNode->nextNode;
    data = node->data;
    free(node);
    queue->size = queue->size - 1;
    return data;
}

void clearQueue(Queue *queue, void (*ptrFreeFun)(void *))
{
    QueueNode *node;

    if(queue!=NULL)
    {
	while(queue->firstNode!=NULL)
	{
	    node = queue->firstNode;
	    queue->firstNode = queue->firstNode->nextNode;
		    		
	    freeQueueNode(node,ptrFreeFun);
	    queue->size = queue->size - 1;
	}	
    }

    queue->lastNode = queue->firstNode = NULL;
}
	
void freeQueue(Queue *queue, void (*ptrFreeFun)(void* ))
{
    QueueNode *node;
    if(queue!=NULL)
    {
	while(queue->firstNode!=NULL)
	{
	    node = queue->firstNode;
	    queue->firstNode = queue->firstNode->nextNode;
		
	    freeQueueNode(node,(*ptrFreeFun));
	    queue->size = queue->size - 1;
	}	
	free(queue);
    }
}
