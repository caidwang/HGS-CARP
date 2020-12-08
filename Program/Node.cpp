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

#include "Node.h"

Node::Node(void) {}

Node::Node(bool isDepot, int cour, int day, bool isPresent, Node * nextNode , Node * pred, Route * route, Params * params)
: isDepot(isDepot),cour(cour),day(day), isPresent(isPresent),nextNode(nextNode), pred(pred), route(route),params(params)
{
	int ccour = cour ;
	if (isDepot) ccour = 0 ;

	// Initialization of the costInsertion structure
	for (int i=0 ; i < params->cli[ccour].visits.size() ; i++)
	{
	    // 对于每种访问pattern，增加一个，如果是CARP，则只有一个pattern
		costInsertion.push_back(1.e30) ;
		placeInsertion.push_back(NULL) ;
	}
	place = -1 ;
}

Node::Node(Node const& copy)
{
	// Copy constructor
	isDepot = copy.isDepot ;
	cour = copy.cour ;
	place = copy.place ;
	day = copy.day ;
	isPresent = copy.isPresent ;
	nextNode = copy.nextNode ;
	pred = copy.pred ;
	route = copy.route ;
	params = copy.params ;
    costInsertion = copy.costInsertion ;
	placeInsertion = copy.placeInsertion ;
	moves = copy.moves ;
}

Node& Node::operator=(Node const& copy)
{
	// Copy constructor
	isDepot = copy.isDepot ;
	cour = copy.cour ;
	place = copy.place ;
	day = copy.day ;
	isPresent = copy.isPresent ;
	nextNode = copy.nextNode ;
	pred = copy.pred ;
	route = copy.route ;
	params = copy.params ;
    costInsertion = copy.costInsertion ;
	placeInsertion = copy.placeInsertion ;
	moves = copy.moves ;
	return *this;
}

Node::~Node(void){}

void Node::setRemaining()
{
	// seq1 has exactly the same meaning than seqi_j[0], but its more convenient to use and read
	seq1 = seqi_j[0] ;
	seq12 = seqi_j[1] ;
	seq123 = seqi_j[2] ;
	seq21 = seqj_i[1] ;
	seq321 = seqj_i[2] ;
}