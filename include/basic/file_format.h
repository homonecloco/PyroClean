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

#ifndef FILE_FORMAT_H_
#define FILE_FORMAT_H_

typedef enum
 {
   FASTA = 0,
   FASTQ = 1,
   CTX   = 2,
   UNSPECIFIED = 3,
   ROCHE   = 4,
   HASH   = 5
 } FileFormat ;

#endif //FILE_FORMAT_H_
