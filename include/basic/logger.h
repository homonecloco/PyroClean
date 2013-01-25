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

/*----------------------------------------------------------------------*
 * File:    logger.c                                                    *
 * Purpose: Handle debugging/log files.                                 *
 * Author:  Richard Leggett                                             *
 *          The Sainsbury Laboratory, Norwich, UK                       *
 *          richard.leggett@tsl.ac.uk                                   *
 * History: 19-Oct-10: RML: Pulled together from ad hoc code!           *
 *----------------------------------------------------------------------*/

int log_start(char* filename);
void log_printf(char* fmt, ...);
void log_and_screen_printf(char* fmt, ...);
void log_newline(void);
void log_get_timestamp(char* buffer);
void log_write_timestamp(int newline);
void log_progress_bar(int percentage);
