// All programs that use Son.H, need to include
// windows.h and windef.h
#include <windows.h>
#include <windef.h>

#include <malloc.h>

// Header files for CED libraries
#include "SON.H"
#include "SONINTL.H"

// Matlab header files
#include <mex.h>
#include <matrix.h>

#define NUMFXNS 8
#define STRING_MATCH 0


/**********************************************************************/
/* void son_error(int)												  */
/*		Prints out a SON library error message to the screen		  */
/**********************************************************************/
void son_error(int);

/************************************************************************************/
/* void son_getMarkData(int, mxArray**, int, const mxArray**)                       */
/*		Gets the marker data from an already opened Spike2 file						*/
/************************************************************************************/
void son_getMarkData(int, mxArray**, int, const mxArray**);

void son_openOldFile(int, mxArray**, int, const mxArray**);
void son_closeFile(int, mxArray**, int, const mxArray**);
void son_chanMaxTime(int, mxArray**, int, const mxArray**);
void son_chanDivide(int, mxArray**, int, const mxArray**);
void son_getADCData(int, mxArray**, int, const mxArray**);
void son_getNumADCPts(int, mxArray**, int, const mxArray**);

/*******************************************************************************/
/*	Returns the file basic time interval.									   */
/*																			   */
/*	Input:	1. Handle to SMR file.											   */
/*																			   */
/*	Output: 1. Number of microseconds in the file basic time interval.		   */
/*******************************************************************************/
void son_getusPerTime(int, mxArray**, int, const mxArray**);

void son_getAllADCData(int, mxArray**, int, const mxArray**);


/*******************************************************************************/
/*	void son_chanKind(int, mxArray**, int, const mxArray**)					   */
/*																			   */
/*******************************************************************************/
void son_chanKind(int, mxArray**, int, const mxArray**);


/*******************************************************************************/
/*	Reads event data between two times into a user defined buffer.			   */
/*																			   */
/*	Input:	1. Handle to SMR file											   */
/*			2. Channel number to read										   */
/*			3. Maximum number of event times to return						   */
/*			4. Start time													   */
/*			5. End time														   */
/*																			   */
/*	Output:	1. Number of event points found									   */
/*			2. Structure array of event times								   */
/*																			   */
/*******************************************************************************/
void son_getEventData(int, mxArray**, int, const mxArray**);