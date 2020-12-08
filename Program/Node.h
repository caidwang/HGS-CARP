/*  ---------------------------------------------------------------------- //
    Hybrid Genetic Search for Arc Routing Problems -- HGS-CARP
    Copyright (C) 2016 Thibaut VIDAL

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//  ---------------------------------------------------------------------- */

#ifndef NODE_H
#define NODE_H

#include <iostream>
using namespace std;
class Route ;
#include "Route.h"

class Node
{

public :

// Access to the data of the problem
Params * params ;

// is-it a depot
bool isDepot ;

// index of the depot or customer
int cour ;

// place in the route
int place ;

// index of the day in which this customer is inserted
int day ;

// is this customer used on this day 
// (all customer nodes are created for each day, but not necessarily inserted in the sequence)
bool isPresent ;

// next depot or client in the route
Node * nextNode ;

// previous depot or client in the route
Node * pred ;

// associated route
Route * route ;

// pointer towards the preprocessed SeqData data structures
// "i" is considered to be the current customer
vector <SeqData *> seqi_j ; // data for (i,j) with j > i
vector <SeqData *> seqj_i ; // data for (j,i) (for the same subsequence as i_j, but reversed)
SeqData * seq0_i ; // data for (0,i)
SeqData * seqi_n ; // data for (i,n), n is the end of the route
SeqData * seqi_0 ; // data for (i,0) (for the reversed route)
SeqData * seqn_i ; // data for (n,i) (for the reversed route)
// the same pointers as (i,j) for some values, but simpler to call
SeqData * seq1 ; // data for (i) 
SeqData * seq12 ; // data for (i,i+1)
SeqData * seq21 ; // data for (i+1,i)
SeqData * seq123 ; // data for (i,i+1,i+2)
SeqData * seq321 ; // data for (i+2,i+1,i)

// cost of insertion in this day, if the considered customer had to be inserted
// This had to be generalized to the PCARP, as the demand may change as a function of the pattern choice, the
// costInsertion can be evaluated for all possible pattern which contain this day.
vector < double > costInsertion ;

// place where it would be inserted
// This had to be generalized to the PCARP, as the demand may change as a function of the pattern choice, the
// placeInsertion can be evaluated for all possible pattern which contain this day.
vector < Node * > placeInsertion ;

// possible moves for this customer and this day (granular search)
// todo moves是在什么时候如何添加的？
vector < int > moves ;

// constructor 1
Node(void);
	
// constructor 2
Node(bool isDepot, int cour, int day, bool isPresent, Node * nextNode , Node * pred, Route * route, Params * params);

// destructor
~Node(void);

// Copy constructor
Node(Node const& copy) ;

// Assignment operator in terms of the copy constructor
Node& operator=(Node const& copy);

// little function to correctly initialize the pointers
void setRemaining();

};

#endif
