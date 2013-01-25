/**   
 *      Copyright 2010 - 2011 by M. Caccamo (mario.caccamo@tgac.ac.uk)  and R. Ramirez-Gonzalez(Ricardo.Ramirez-Gonzalez@tgac.ac.uk)
 *      
 *      PyroCleanis free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 3 of the License, or
 *      (at your option) any later version.
 *      
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *      
 *      You should have received a copy of the GNU General Public License
 *      along with this program; if not, write to the Free Software
 *      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *      MA 02110-1301, USA.
 *
 */
 
#ifndef FLAGS_H_
#define FLAGS_H_
#include <stdlib.h>
#include <stdio.h>
#include <global.h>

//Global flags
#define  ALL_OFF  		   	 0

/**
 * This will be used to mark from which paths the node has been traversed
 * The list of defines give the value for each mask.
 *
 */
#ifndef SHORT_FLAGS
typedef unsigned int Flags;

#else
typedef unsigned short Flags;

#endif

boolean flags_check_for_flag(Flags f, Flags * db);

boolean flags_check_for_any_flag(Flags f, Flags * db);

void flags_action_clear_flags(Flags * db);

void flags_action_set_flag(Flags f, Flags * db);

void flags_action_unset_flag(Flags f, Flags * db);

#endif /* FLAGS_H_ */
