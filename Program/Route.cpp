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

#include "Route.h"

Route::Route(void){}

Route::Route(int cour, Node * depot, Vehicle * vehicle, Params * params, Individual * indiv, int day) : cour(cour), depot(depot), vehicle(vehicle), params(params), individu(indiv), day(day)
{
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i++ )
	{
		coutInsertionClient.push_back(vector <double> ()) ;
		placeInsertionClient.push_back(vector <Node *> ());
		for (int p=0 ; p < params->cli[i].visits.size() ; p++)
		{
			coutInsertionClient[i].push_back(1.e30);
			placeInsertionClient[i].push_back(NULL);
		}
		nodeAndRouteTested.push_back(false);
	}
}

Route::~Route(void){}

void Route::reverse ()
{
	// Reversing the order of the nodes in the sequence
	if (!depot->nextNode->isDepot)
	{
		Node * temp ;
		Node * myDepot = depot ;
		Node * myDepotFin = depot->pred ;
		Node * cour = depot->nextNode ;

		while ( !cour->isDepot )
		{
			temp = cour->nextNode ;
			cour->nextNode = cour->pred ;
			cour->pred = temp ;
			cour = temp ;
		}

		temp = myDepot->nextNode ;
		myDepot->nextNode = myDepotFin->pred ;
		myDepotFin->pred = temp ;
		myDepot->nextNode->pred = myDepot ;
		myDepotFin->pred->nextNode = myDepotFin ;
	}
}

// 初始化了每个节点的预处理信息，包括从0->i, i->0, i->n, n->i, 还有i->j
// CARP的子问题的求解在这部分完成
void Route::updateRouteData (bool isForPrinting)
{
    /*
     * 1。 从0仓点出发，对0-i，i-0的子串进行预处理
     * 2。 从n仓点出发，对i->n,n->i的子串进行预处理
     * 3。 根据sizeSD的限定，从0仓点开始处理i->j, j->i结构的预处理
     * 4。 判断解的可行性
     */
    bool firstIt ;
    int place = 0 ;
    double Xvalue = 0 ;
    double Yvalue = 0 ;
    double nbNodes = 1 ;

    // Computing the auxiliary data on any subsequence (0..i), using forward recursion
    Node * node = depot ;
    node->place = place ;
    node->seq0_i->initialisation(node->cour, params, individu, day, isForPrinting);
    node->seqi_0->initialisation(node->cour, params, individu, day, false);
    Xvalue += params->cli[node->cour].coord.x ;
    Yvalue += params->cli[node->cour].coord.y ;

    firstIt = true ;
    while (!node->isDepot || firstIt) {
        firstIt = false ;
        node = node->nextNode ;
        Xvalue += params->cli[node->cour].coord.x ;
        Yvalue += params->cli[node->cour].coord.y ;
        nbNodes ++ ;
        place ++ ;
        node->place = place ;
        if (isForPrinting)
            node->seq0_i->concatOneAfterWithPathTracking(node->pred->seq0_i, node->cour, individu, day);
        else
            node->seq0_i->concatOneAfter(node->pred->seq0_i, node->cour, individu, day);
        node->seqi_0->concatOneBefore(node->pred->seqi_0, node->cour, individu, day);
    }

    // Computing the auxiliary data on any subsequence (i..n), using backward recursion
    node = depot->pred ;
    node->seqi_n->initialisation(node->cour, params, individu, day, false);
    node->seqn_i->initialisation(node->cour, params, individu, day, false);

    firstIt = true ;
    while (!node->isDepot || firstIt ) {
        firstIt = false ;
        node = node->pred ;
        node->seqi_n->concatOneBefore(node->nextNode->seqi_n, node->cour, individu, day);
        node->seqn_i->concatOneAfter(node->nextNode->seqn_i, node->cour, individu, day);
    }

    // Computing the auxiliary data on any subsequence (i..j), using forward recursion
    // To gain a bit of time, we limit this preprocessing to subsequences such that i..j does not contain more than "sizeSD" elements
    // (More intelligent strategies could be used (e.g., the hierarchical approach of Irnich))
    Node * noeudi ;
    Node * noeudj ;
    noeudi = depot ;
    for (int i=0 ; i <= depot->pred->place ; i++) {
        noeudi->seqi_j[0]->initialisation(noeudi->cour,params,individu,day,false);
        noeudj = noeudi->nextNode ;
        for (int j=1 ; j <= depot->pred->place - i && j < params->sizeSD ; j++) {
            noeudi->seqi_j[j]->concatOneAfter(noeudi->seqi_j[j-1],noeudj->cour,individu,day);
            noeudj = noeudj->nextNode ;
        }
        noeudi = noeudi->nextNode ;
    }

    noeudi = depot->pred ;
    for (int i=0 ; i <= depot->pred->place ; i++) {
        noeudi->seqj_i[0]->initialisation(noeudi->cour,params,individu,day,false);
        noeudj = noeudi->pred ;
        for (int j=1 ; j <= depot->pred->place - i && j < params->sizeSD ; j++) {
            noeudj->seqj_i[j]->concatOneAfter(noeudj->nextNode->seqj_i[j - 1], noeudj->cour, individu, day);
            noeudj = noeudj->pred ;
        }
        noeudi = noeudi->pred ;
    }

    // Checking the route feasibility
    double dist ;
    double violLoad ;
    double violDuration ;
    currentRouteCost = depot->pred->seq0_i->evaluation(depot->pred->seq0_i,depot->seq1,this->vehicle,dist,violDuration,violLoad);

    if (violDuration < 0.001 && violLoad < 0.001)
        isFeasible = true ;
    else
        isFeasible = false ;
}

// no insertion are computed
void Route::initiateInsertions()
{
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i ++ )
		for (int p=0 ; p < (int)params->cli[i].visits.size() ; p++)
			coutInsertionClient[i][p] = 1.e30 ;
}

