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


#ifndef _QUEUE_H_

	#define _QUEUE_H_

	#include"types.h"

	QueueNode* createNode(void);
	void       freeQueueNode(QueueNode *queueNode, void ( *funPtr)(void *data));
	Queue*     createQueue(void);
	Queue*     addNode(Queue *queue, void *data);
	void*      readFirstNode(Queue *queue);
	void*      readLastNode(Queue *queue);
	void*      readNNode(Queue *queue, unsigned int nodeNumber);
	void*      substractNode(Queue *queue);
	void       clearQueue(Queue *queue, void (*ptrFreeFun)(void* ));
	void       freeQueue(Queue *queue, void (*ptrFreeFun)(void* ));
	
#endif
