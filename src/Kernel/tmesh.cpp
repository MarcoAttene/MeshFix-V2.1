/****************************************************************************
* TMesh                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2012: IMATI-GE / CNR                                         *
* All rights reserved.                                                      *
*                                                                           *
* This program is dual-licensed as follows:                                 *
*                                                                           *
* (1) You may use TMesh as free software; you can redistribute it and/or *
* modify it under the terms of the GNU General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or      *
* (at your option) any later version.                                       *
* In this case the program is distributed in the hope that it will be       *
* useful, but WITHOUT ANY WARRANTY; without even the implied warranty of    *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
* (2) You may use TMesh as part of a commercial software. In this case a *
* proper agreement must be reached with the Authors and with IMATI-GE/CNR   *
* based on a proper licensing contract.                                     *
*                                                                           *
****************************************************************************/

#include "tmesh_kernel.h"
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

namespace T_MESH
{

extern "C" void initPredicates();

void (* TMesh::display_message)(const char*, int) = NULL;

char *TMesh::app_name = NULL;
char *TMesh::app_version = NULL;
char *TMesh::app_year = NULL;
char *TMesh::app_authors = NULL;
char *TMesh::app_url = NULL;
char *TMesh::app_maillist = NULL;
const char *TMesh::filename = NULL;
bool TMesh::quiet = false;

void TMesh::init(void (*dm)(const char *, int))
{
 display_message = dm;
 app_name = NULL;
 app_version = NULL;
 app_year = NULL;
 app_authors = NULL;
 app_url = NULL;
 app_maillist = NULL;
 filename = NULL;
 quiet = false;
 initPredicates();
}


///////////// Prints a fatal error message and exits /////////////

void TMesh::error(const char *msg, ...)
{
 static char fmt[2048], fms[4096];
 va_list ap;
 va_start(ap, msg);
 strcpy(fmt,"\nERROR- ");
 strcat(fmt,msg);
 vsprintf(fms,fmt,ap);

 if (display_message != NULL)
  display_message(fms, DISPMSG_ACTION_ERRORDIALOG);
 else
 {
  fprintf(stderr,fms);
  exit(-1);
 }
}

///////////// Prints a warning message /////////////

void TMesh::warning(const char *msg, ...)
{
 if (quiet) return;
 static char fmt[2048], fms[4096];
 va_list ap;
 va_start(ap, msg);
 strcpy(fmt,"WARNING- ");
 strcat(fmt,msg);
 vsprintf(fms,fmt,ap);

 if (display_message != NULL) 
  display_message(fms, DISPMSG_ACTION_PUTMESSAGE);
 else
  fputs(fms, stderr);

 va_end(ap);
}

///////////// Prints an information message /////////////

void TMesh::info(const char *msg, ...)
{
 if (quiet) return;
 static char fmt[2048], fms[4096];
 va_list ap;
 va_start(ap, msg);
 strcpy(fmt,"INFO- ");
 strcat(fmt,msg);
 vsprintf(fms,fmt,ap);

 if (display_message != NULL) 
  display_message(fms, DISPMSG_ACTION_PUTMESSAGE);
 else
  printf(fms);

 va_end(ap);
}

///////// Reports progress status for a process //////////

void TMesh::begin_progress()
{
 if (quiet) return;
 if (display_message != NULL) 
  display_message("\n", DISPMSG_ACTION_PUTNEWLINE);
 else
  printf("\n");
}

void TMesh::report_progress(const char *msg, ...)
{
 if (quiet) return;
 static char fmt[2048] = "\r";
 static char fms[4096];
 static char rotating_bar[5] = "-\\|/";
 static unsigned char wc=0;

 if (msg == NULL)
 {
  sprintf(fms,"%c",rotating_bar[wc++]); if (wc==4) wc=0;
  strcpy(fmt+1,fms);

  if (display_message != NULL) 
   display_message(fmt, DISPMSG_ACTION_PUTPROGRESS);
  else
  {
   printf("%s",fmt);
   fflush(stdout);
  }
 }
 else
 {
  va_list ap;
  va_start(ap, msg);
  strcpy(fmt+1,msg);
  vsprintf(fms,fmt,ap);

  if (display_message != NULL) 
   display_message(fms, DISPMSG_ACTION_PUTPROGRESS);
  else
  {
   printf("%s", fms);
   fflush(stdout);
  }
  va_end(ap);
 }
}

void TMesh::end_progress()
{
 if (quiet) return;
 if (display_message != NULL) 
  display_message("\n", DISPMSG_ACTION_PUTNEWLINE);
 else
  printf("\n");
}

void TMesh::useRationals(bool u) 
{
#ifdef USE_HYBRID_KERNEL
	coord::use_rationals = u;
#endif
}

bool TMesh::isUsingRationals()
{
#ifdef USE_HYBRID_KERNEL
	return coord::use_rationals;
#else
	return false;
#endif
}

void TMesh::addMessageToLogFile(const char *msg)
{
	FILE *fp = fopen("tmesh.log", "a");
	fprintf(fp, msg);
	fclose(fp);
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
char *currentDateTime()
{
	time_t     now = time(0);
	struct tm  tstruct;
	static char buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	return buf;
}

void TMesh::logToFileAndExit(const char *s)
{
	static char msg[2048];
	sprintf(msg, "%s\nFILE: %s\nRETURN VALUE: %s\n\n", currentDateTime(), (filename) ? (filename) : ("unknown"), s);
	addMessageToLogFile(msg);
	TMesh::error(msg);
}

void TMesh::exitOnTimeout(clock_t ts)
{
	static clock_t beginning_time, timeout_secs;
	if (ts != 0) { beginning_time = clock(); timeout_secs = ts; }
	else if (((clock() - beginning_time) / 1000) > timeout_secs) logToFileAndExit("Timeout reached");
}

void TMesh::printElapsedTime(bool reset)
{
	static clock_t beginning_time;
	if (reset) beginning_time = clock();
	else printf("\n\n********** PARTIAL ELAPSED: %d msecs\n\n", (clock() - beginning_time));
}

} //namespace T_MESH
