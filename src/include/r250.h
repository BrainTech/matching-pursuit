/***************************************************************************
 *   Copyright (C) 2006 by Piotr J. Durka Dobieslaw Ircha, Rafal Kus       *
 *   durka@fuw.edu.pl, rircha@fuw.edu.pl, rkus@fuw.edu.pl                  *
 *   Department of Biomedical Physics at Warsaw University                 *
 *   http://brain.fuw.edu.pl, http://eeg.pl                                *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef _RAND_H_

	#define _RAND_H_

	/**** Function: r250_init
        	Description: initializes r250 random number generator. ****/

	void r250_init(int seed);

	/**** Function: r250 Description: returns a random unsigned integer k
			uniformly distributed in the interval 0 <= k < 65536.  ****/

	unsigned int r250(void);

	/**** Function: r250n Description: returns a random unsigned integer k
			uniformly distributed in the interval 0 <= k < n ****/
	unsigned int r250n(unsigned n);

	/**** Function: dr250
                Description: returns a random double z in range 0 <= z < 1.  ****/

	double dr250(void);

#endif
